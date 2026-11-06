"""Per-sample SQLite index of every small variant call, annotated at build time.

A PacBio genome has 4-5 million small-variant calls, so the review screens cannot hold them
as Python objects or re-scan a VCF per request. The VCF is read once, each call is
annotated (gene body, ClinVar, gnomAD) and written to a SQLite file under the cache
directory; the API then filters, counts and pages with SQL. The index is derived data from a
real sample, so it lives outside the repo and is rebuilt whenever the VCF or either
annotation resource changes.

Annotation is by exact allele on trimmed alleles (see normalize.py). gnomAD frequencies are
joined inside SQLite from a reduced frequency database (see gnomad.py).
"""

from __future__ import annotations

import json
import os
import sqlite3
import time
from collections import Counter
from pathlib import Path
from typing import Dict, Optional

from . import vcf as vcf_mod
from .clinvar import ClinvarTable
from .genes import GeneIndex
from .gnotate import Gnotate
from .normalize import CHROM_ORDER, chrom_order

BUILDER_VERSION = 2
BATCH = 50_000

# GRCh38 pseudo-autosomal regions on chrX (1-based, inclusive)
PAR_X = ((10_001, 2_781_479), (155_701_383, 156_030_895))

SCHEMA = """
CREATE TABLE variants (
    id INTEGER PRIMARY KEY,
    chrom TEXT NOT NULL, chrom_order INTEGER NOT NULL, pos INTEGER NOT NULL,
    ref TEXT NOT NULL, alt TEXT NOT NULL, kind TEXT NOT NULL,
    npos INTEGER NOT NULL, nref TEXT NOT NULL, nalt TEXT NOT NULL,
    copies INTEGER NOT NULL, gt TEXT NOT NULL, phased INTEGER NOT NULL, ps TEXT, hap TEXT,
    qual REAL, filter TEXT, passed INTEGER NOT NULL, dp INTEGER, gq INTEGER, vaf REAL,
    genes TEXT,
    cv_group TEXT, cv_class TEXT, cv_stars INTEGER, cv_review TEXT, cv_id TEXT, cv_name TEXT,
    cv_cond TEXT, cv_evaluated TEXT, cv_nsub INTEGER, cv_gene TEXT,
    gn_af REAL, gn_ac INTEGER, gn_hom INTEGER,
    hp_af REAL, hp_ac INTEGER, hp_hom INTEGER
);
CREATE TABLE variant_genes (variant_id INTEGER NOT NULL, gene TEXT NOT NULL COLLATE NOCASE);
CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);
"""

INDEXES = """
CREATE INDEX idx_pos ON variants (chrom_order, pos);
CREATE INDEX idx_norm ON variants (chrom, npos, nref, nalt);
CREATE INDEX idx_cv ON variants (cv_group);
CREATE INDEX idx_af ON variants (gn_af);
CREATE INDEX idx_gene ON variant_genes (gene);
"""

_INSERT_COLUMNS = (
    "chrom, chrom_order, pos, ref, alt, kind, npos, nref, nalt, copies, gt, phased, ps, hap, qual, filter, "
    "passed, dp, gq, vaf, genes, cv_group, cv_class, cv_stars, cv_review, cv_id, cv_name, cv_cond, "
    "cv_evaluated, cv_nsub, cv_gene"
)
_PLACEHOLDERS = ",".join("?" * 31)


class Progress:
    """Read by the API while a build runs on another thread."""

    def __init__(self) -> None:
        self.phase = "starting"
        self.records = 0
        self.started = time.time()

    def as_dict(self) -> Dict[str, object]:
        return {"phase": self.phase, "records": self.records, "seconds": round(time.time() - self.started, 1)}


class QcAccumulator:
    """Whole-genome call metrics, accumulated as calls stream past."""

    INDEL_BINS = ((1, "1"), (2, "2"), (3, "3"), (4, "4"), (5, "5"), (10, "6-10"), (20, "11-20"), (49, "21-49"))

    def __init__(self) -> None:
        self.kind = Counter()
        self.zyg = Counter()
        self.filtered = 0
        self.ti = self.tv = 0
        self.het = self.phased_het = 0
        self.indel_hist = Counter()
        self.x_het = self.x_hom = self.auto_het = self.auto_hom = 0
        self.dp_sum = self.dp_n = self.gq_sum = self.gq_n = 0
        self.per_chrom = Counter()

    def add(self, c: vcf_mod.Call, passed: bool) -> None:
        if not passed:
            self.filtered += 1
            return
        self.kind[c.kind] += 1
        het = c.copies == 1
        self.zyg["het" if het else "hom"] += 1
        self.per_chrom[c.chrom] += 1
        if c.kind == "snv":
            pair = (c.norm_ref, c.norm_alt)
            if pair in (("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")):
                self.ti += 1
            else:
                self.tv += 1
        elif c.kind in ("ins", "del"):
            size = max(len(c.norm_ref), len(c.norm_alt))
            for limit, label in self.INDEL_BINS:
                if size <= limit:
                    self.indel_hist[label] += 1
                    break
        if het:
            self.het += 1
            if c.phased:
                self.phased_het += 1
        if c.dp is not None:
            self.dp_sum += c.dp
            self.dp_n += 1
        if c.gq is not None:
            self.gq_sum += c.gq
            self.gq_n += 1
        if c.chrom == "chrX" and not any(lo <= c.pos <= hi for lo, hi in PAR_X):
            if het:
                self.x_het += 1
            else:
                self.x_hom += 1
        elif c.chrom not in ("chrX", "chrY", "chrM"):
            if het:
                self.auto_het += 1
            else:
                self.auto_hom += 1

    def finish(self, stats: vcf_mod.ReadStats) -> Dict[str, object]:
        passed = sum(self.kind.values())
        auto_total = self.auto_het + self.auto_hom
        x_total = self.x_het + self.x_hom
        sex, basis = infer_sex(self.x_het, self.x_hom, self.auto_het, self.auto_hom)
        return {
            "records": stats.records,
            "calls_passing": passed,
            "calls_filtered": self.filtered,
            "snv": self.kind["snv"],
            "mnv": self.kind["mnv"],
            "insertions": self.kind["ins"],
            "deletions": self.kind["del"],
            "complex": self.kind["complex"],
            "het": self.zyg["het"],
            "hom": self.zyg["hom"],
            "het_hom_ratio": round(self.zyg["het"] / self.zyg["hom"], 3) if self.zyg["hom"] else None,
            "ti_tv": round(self.ti / self.tv, 3) if self.tv else None,
            "phased_het_fraction": round(self.phased_het / self.het, 4) if self.het else None,
            "mean_dp": round(self.dp_sum / self.dp_n, 1) if self.dp_n else None,
            "mean_gq": round(self.gq_sum / self.gq_n, 1) if self.gq_n else None,
            "indel_sizes": [{"bin": label, "count": self.indel_hist[label]} for _, label in self.INDEL_BINS],
            "per_chrom": [{"chrom": c, "count": self.per_chrom[c]} for c in sorted(self.per_chrom, key=chrom_order)],
            "symbolic_or_star_alleles": stats.symbolic_or_star,
            "sample_columns": stats.samples,
            "sample_name": stats.sample_name,
            "x_nonpar_calls": x_total,
            "autosomal_calls": auto_total,
            "x_nonpar_het_fraction": round(self.x_het / x_total, 4) if x_total else None,
            "autosomal_het_fraction": round(self.auto_het / auto_total, 4) if auto_total else None,
            "inferred_sex": sex,
            "inferred_sex_basis": basis,
        }


def infer_sex(x_het: int, x_hom: int, auto_het: int, auto_hom: int, min_x: int = 500):
    """Male / female / undetermined from the call set itself.

    A male carries one X outside the pseudo-autosomal regions, so almost none of his calls
    there are heterozygous; a female's heterozygous fraction there is close to her
    autosomal one. The ratio of the two fractions separates them when there are enough X
    calls (a truth set limited to chr1-22, or a small slice, will not have them)."""
    x_total, auto_total = x_het + x_hom, auto_het + auto_hom
    if x_total < min_x or auto_total == 0 or auto_het == 0:
        return "undetermined", f"too few chrX non-PAR calls ({x_total}) to infer sex"
    ratio = (x_het / x_total) / (auto_het / auto_total)
    basis = f"chrX non-PAR het fraction is {ratio:.2f}x the autosomal fraction ({x_total:,} X calls)"
    if ratio < 0.25:
        return "male", basis
    if ratio > 0.6:
        return "female", basis
    return "undetermined", basis


def _describe(g: Optional[Gnotate]) -> Dict[str, object]:
    if g is None:
        return {}
    return {"file": os.path.basename(g.path), "created": g.created, "fields": g.fields, "args": g.args}


def _field(g: Gnotate, suffix: str) -> str:
    for f in g.fields:
        if f.endswith("_" + suffix):
            return f
    raise ValueError(f"{g.path} has no *_{suffix} field (has {g.fields})")


def _join_frequencies(con: sqlite3.Connection, g: Gnotate, prefix: str, label: str, progress: Progress) -> int:
    """Fill ``<prefix>_af/_ac/_hom`` for every call found in a gnotate frequency zip.

    Matching is on the alleles as written in the VCF (position, ref, alt), the way gnomAD
    and slivar write them, one chromosome at a time."""
    af, ac, hom = _field(g, "af"), _field(g, "ac"), _field(g, "nhomalt")
    n = 0
    chroms = [r[0] for r in con.execute("SELECT DISTINCT chrom FROM variants ORDER BY chrom_order")]
    for chrom in chroms:
        if not g.has_chrom(chrom):
            continue
        progress.phase = f"joining {label} frequencies ({chrom})"
        rows = con.execute("SELECT id, pos, upper(ref), upper(alt) FROM variants WHERE chrom = ? ORDER BY id", (chrom,)).fetchall()
        found, values = g.lookup(chrom, [r[1] for r in rows], [r[2] for r in rows], [r[3] for r in rows])
        hits = [i for i, ok in enumerate(found) if ok]
        con.executemany(
            f"UPDATE variants SET {prefix}_af = ?, {prefix}_ac = ?, {prefix}_hom = ? WHERE id = ?",
            [(values[af][i], int(round(values[ac][i])), int(round(values[hom][i])), rows[i][0]) for i in hits],
        )
        n += len(hits)
    con.commit()
    return n


def build_index(
    vcf_path: str,
    db_path: str,
    genes: Optional[GeneIndex],
    clinvar: Optional[ClinvarTable],
    gnomad: Optional[Gnotate] = None,
    hprc: Optional[Gnotate] = None,
    progress: Optional[Progress] = None,
) -> Dict[str, object]:
    progress = progress or Progress()
    tmp = f"{db_path}.building"
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)
    if os.path.exists(tmp):
        os.remove(tmp)
    con = sqlite3.connect(tmp)
    con.executescript("PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF;" + SCHEMA)

    stats = vcf_mod.ReadStats()
    qc = QcAccumulator()
    progress.phase = "reading calls"
    batch, gene_rows, next_id = [], [], 1
    n_clinvar = 0
    for c in vcf_mod.iter_calls(vcf_path, stats):
        passed = c.filter in vcf_mod.PASS_FILTERS
        qc.add(c, passed)
        symbols = ()
        if genes is not None:
            start0 = c.norm_pos - 1
            symbols = genes.genes_at(c.chrom, start0, start0 + max(len(c.norm_ref), 1))
        cv = clinvar.lookup(c.chrom, c.norm_pos, c.norm_ref, c.norm_alt) if clinvar is not None else None
        if cv is not None:
            n_clinvar += 1
        batch.append(
            (
                c.chrom, chrom_order(c.chrom), c.pos, c.ref, c.alt, c.kind, c.norm_pos, c.norm_ref, c.norm_alt,
                c.copies, c.gt, int(c.phased), c.ps, c.haplotypes, c.qual, c.filter, int(passed), c.dp, c.gq, c.vaf,
                ",".join(symbols) or None,
                *(
                    (cv.group, cv.classification, cv.stars, cv.review, cv.variation_id, cv.name, cv.conditions,
                     cv.last_evaluated, cv.n_submitters, cv.gene)
                    if cv is not None
                    else (None,) * 10
                ),
            )
        )
        for g in symbols:
            gene_rows.append((next_id, g))
        next_id += 1
        progress.records = stats.records
        if len(batch) >= BATCH:
            con.executemany(f"INSERT INTO variants ({_INSERT_COLUMNS}) VALUES ({_PLACEHOLDERS})", batch)
            con.executemany("INSERT INTO variant_genes VALUES (?, ?)", gene_rows)
            batch, gene_rows = [], []
    if batch:
        con.executemany(f"INSERT INTO variants ({_INSERT_COLUMNS}) VALUES ({_PLACEHOLDERS})", batch)
        con.executemany("INSERT INTO variant_genes VALUES (?, ?)", gene_rows)

    progress.phase = "indexing"
    con.executescript(INDEXES)

    n_gnomad = _join_frequencies(con, gnomad, "gn", "gnomAD", progress) if gnomad is not None else 0
    n_hprc = _join_frequencies(con, hprc, "hp", "HPRC", progress) if hprc is not None else 0
    if n_gnomad or n_hprc:
        progress.phase = "indexing"
        con.execute("CREATE INDEX idx_hp ON variants (hp_af)")

    qc_result = qc.finish(stats)
    total = con.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
    st = os.stat(vcf_path)
    meta = {
        "builder_version": str(BUILDER_VERSION),
        "vcf_path": str(vcf_path),
        "vcf_stamp": f"{st.st_size}:{st.st_mtime_ns}",
        "clinvar_stamp": clinvar.stamp if clinvar is not None else "",
        "clinvar_release": clinvar.release if clinvar is not None else "",
        "clinvar_rows": str(len(clinvar)) if clinvar is not None else "0",
        "gnomad_stamp": gnomad.stamp if gnomad is not None else "",
        "gnomad_info": json.dumps(_describe(gnomad)),
        "hprc_stamp": hprc.stamp if hprc is not None else "",
        "hprc_info": json.dumps(_describe(hprc)),
        "genes_loaded": str(genes.n_genes) if genes is not None else "0",
        "n_variants": str(total),
        "n_clinvar_hits": str(n_clinvar),
        "n_gnomad_hits": str(n_gnomad),
        "n_hprc_hits": str(n_hprc),
        "qc": json.dumps(qc_result),
        "built": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    con.executemany("INSERT INTO meta VALUES (?, ?)", list(meta.items()))
    con.commit()
    con.close()
    os.replace(tmp, db_path)
    progress.phase = "ready"
    return meta


def read_meta(db_path: str) -> Optional[Dict[str, str]]:
    if not os.path.exists(db_path):
        return None
    try:
        con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
        try:
            return dict(con.execute("SELECT key, value FROM meta").fetchall())
        finally:
            con.close()
    except sqlite3.Error:
        return None


def is_current(
    meta: Optional[Dict[str, str]], vcf_path: str, clinvar_stamp: str, gnomad_stamp: str, hprc_stamp: str
) -> bool:
    if not meta or meta.get("builder_version") != str(BUILDER_VERSION):
        return False
    st = os.stat(vcf_path)
    return (
        meta.get("vcf_stamp") == f"{st.st_size}:{st.st_mtime_ns}"
        and meta.get("clinvar_stamp", "") == clinvar_stamp
        and meta.get("gnomad_stamp", "") == gnomad_stamp
        and meta.get("hprc_stamp", "") == hprc_stamp
    )
