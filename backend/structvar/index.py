"""Per-callset SQLite index of every SV call, annotated at build time.

A call set is one SV VCF, usually a joint call across a family (sawfish, pbsv or sniffles on a trio). The
VCF is read once; each call gets its gene impact (against the gene models), population matches (against
the SV catalogs), segmental-duplication / repeat / odd-region flags, and, for a trio, an inheritance class.
The API then filters, counts and pages with SQL. The index is derived from real participants' data, so it
lives outside the repo and is rebuilt when the VCF or any annotation resource changes.

Who is the child: a ``sample.json`` can name the proband and parents. Otherwise, for exactly three samples,
the child is the one whose genotypes are most consistent with the other two as its parents (the fewest
Mendelian violations on PASS calls made in all three), accepted only if it is clearly better than the
alternatives, and the basis is stored and shown. With any other number of samples no roles are assumed.
"""

from __future__ import annotations

import json
import os
import sqlite3
import time
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from . import vcf as vcf_mod
from .catalogs import Catalogs, Match
from .genes import GeneModels, impacts
from .normalize import SIZE_BINS, LAST_BIN, SV_TYPES, chrom_order, size_bin
from .tracks import Tracks

BUILDER_VERSION = 4
BATCH = 20_000
MAX_GENES_JSON = 25
MIN_TRIO_SITES = 200
TRIO_MARGIN = 0.5  # the best candidate's violation rate must be under this fraction of the next best

SCHEMA = """
CREATE TABLE svs (
    id INTEGER PRIMARY KEY,
    vid TEXT, chrom TEXT NOT NULL, chrom_order INTEGER NOT NULL, pos INTEGER NOT NULL, end INTEGER NOT NULL,
    svtype TEXT NOT NULL, size INTEGER, size_bin TEXT, filter TEXT, passed INTEGER NOT NULL, qual REAL,
    imprecise INTEGER NOT NULL, symbolic INTEGER NOT NULL, mate_chrom TEXT, mate_pos INTEGER, event TEXT,
    f_gt TEXT, f_copies INTEGER, f_ploidy INTEGER, f_gq INTEGER, f_alt INTEGER, f_dp INTEGER, f_cn REAL,
    p1_copies INTEGER, p1_dp INTEGER, p2_copies INTEGER, p2_dp INTEGER, inh TEXT, p_min_dp INTEGER,
    family TEXT,
    n_genes INTEGER NOT NULL, genes TEXT, impact TEXT NOT NULL, impact_rank INTEGER NOT NULL, gene_json TEXT,
    pop_n INTEGER NOT NULL, pop_af_max REAL, gnomad_af REAL, gnomad_popmax REAL, hprc_af REAL, pop_json TEXT,
    segdup_frac REAL, odd TEXT, repeat TEXT, repeat_frac REAL
);
CREATE TABLE sv_genes (sv_id INTEGER NOT NULL, symbol TEXT NOT NULL COLLATE NOCASE, rank INTEGER NOT NULL, coding INTEGER NOT NULL);
CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);
"""

INDEXES = """
CREATE INDEX idx_pos ON svs (chrom, pos);
CREATE INDEX idx_order ON svs (chrom_order, pos);
CREATE INDEX idx_type ON svs (svtype);
CREATE INDEX idx_impact ON svs (impact_rank);
CREATE INDEX idx_inh ON svs (inh);
CREATE INDEX idx_gene ON sv_genes (symbol);
"""


class Progress:
    def __init__(self) -> None:
        self.phase = "starting"
        self.records = 0
        self.started = time.time()

    def as_dict(self) -> Dict[str, object]:
        return {"phase": self.phase, "records": self.records, "seconds": round(time.time() - self.started, 1)}


# ---------------------------------------------------------------- family structure
def _copies(g: vcf_mod.Genotype) -> int:
    return g.copies


def _violation(c: int, p1: int, p2: int) -> bool:
    """A child genotype that its parents' genotypes cannot produce (copies of the ALT allele)."""
    return (c > 0 and p1 + p2 == 0) or (c == 0 and p1 == 2 and p2 == 2) or (c == 2 and (p1 == 0 or p2 == 0))


def infer_trio(calls: List[vcf_mod.SvCall], names: List[str]) -> Dict[str, Any]:
    """Pick the child of a three-sample callset from Mendelian consistency; see the module docstring."""
    out: Dict[str, Any] = {"proband": None, "parents": [], "source": "none", "basis": "no family structure assumed"}
    if len(names) != 3:
        return out
    viol, total = [0, 0, 0], 0
    for c in calls:
        if not c.passed or c.svtype == "BND":
            continue
        cp = [g.copies for g in c.samples]
        if min(cp) < 0:
            continue
        total += 1
        for child in range(3):
            p1, p2 = (cp[i] for i in range(3) if i != child)
            viol[child] += _violation(cp[child], p1, p2)
    if total < MIN_TRIO_SITES:
        out["basis"] = f"too few PASS calls made in all three samples ({total}) to infer who the child is"
        return out
    rates = [v / total for v in viol]
    order = sorted(range(3), key=lambda i: rates[i])
    best, second = order[0], order[1]
    detail = "; ".join(f"{names[i]} {100 * rates[i]:.1f}%" for i in range(3))
    if rates[second] > 0 and rates[best] <= TRIO_MARGIN * rates[second]:  # a tie (e.g. all zero) is not evidence
        out.update(
            proband=best, parents=[i for i in range(3) if i != best], source="inferred",
            basis=f"{names[best]} is the child: Mendelian violation rate if each sample were the child is {detail} ({total:,} PASS calls made in all three)",
        )
    else:
        out["basis"] = f"no sample is clearly the child; violation rate if each were the child: {detail} ({total:,} PASS calls)"
    return out


def resolve_roles(calls: List[vcf_mod.SvCall], names: List[str], declared: Dict[str, Any]) -> Dict[str, Any]:
    """Roles from sample.json when it names them, otherwise inferred."""
    proband = declared.get("proband")
    if proband in names:
        parents = [names.index(p) for p in declared.get("parents", []) if p in names]
        return {"proband": names.index(proband), "parents": parents[:2], "source": "sample.json",
                "basis": f"stated in sample.json: {proband} is the proband" + (f", parents {', '.join(names[i] for i in parents)}" if parents else "")}
    roles = infer_trio(calls, names)
    if roles["proband"] is None:
        roles["proband"] = 0 if names else None  # the first sample is shown; no family structure is assumed
    return roles


def inheritance(child: int, p1: Optional[int], p2: Optional[int]) -> Optional[str]:
    """de_novo | inherited | unknown for a call the child carries; None if the child does not carry it."""
    if child <= 0:
        return None
    ps = [p for p in (p1, p2) if p is not None]
    if len(ps) < 2:
        return None
    if any(p > 0 for p in ps):
        return "inherited"
    return "de_novo" if all(p == 0 for p in ps) else "unknown"


# ---------------------------------------------------------------- build
def _pop_summary(matches: Dict[str, Match]) -> Tuple[Optional[float], Optional[float], Optional[float], Optional[float]]:
    g, h = matches.get("gnomad_sv"), matches.get("hprc")
    afs = [x for x in (g.af if g else None, g.popmax if g else None, h.af if h else None) if x is not None]
    return (max(afs) if afs else None), (g.af if g else None), (g.popmax if g else None), (h.af if h else None)


def build_index(
    vcf_path: str, db_path: str, genes: Optional[GeneModels], catalogs: Optional[Catalogs], tracks: Optional[Tracks],
    declared: Optional[Dict[str, Any]] = None, progress: Optional[Progress] = None,
) -> Dict[str, str]:
    progress = progress or Progress()
    declared = declared or {}
    progress.phase = "reading calls"
    stats = vcf_mod.ReadStats()
    calls: List[vcf_mod.SvCall] = []
    for c in vcf_mod.iter_calls(vcf_path, stats):
        calls.append(c)
        progress.records = stats.records
    names = stats.sample_names
    roles = resolve_roles(calls, names, declared)
    focus: Optional[int] = roles["proband"]
    parents: List[int] = roles["parents"]
    role_of = {i: ("proband" if i == focus else "parent" if i in parents else "member") for i in range(len(names))}

    tmp = f"{db_path}.building"
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)
    if os.path.exists(tmp):
        os.remove(tmp)
    con = sqlite3.connect(tmp)
    con.executescript("PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF;" + SCHEMA)

    progress.phase = "annotating"
    rows: List[tuple] = []
    gene_rows: List[tuple] = []
    repeat_spans: List[Tuple[int, str, int, int]] = []
    qc_pass_types: Counter = Counter()
    qc_all_types: Counter = Counter()
    qc_filters: Counter = Counter()
    qc_sizes: Dict[str, Counter] = {}
    per_sample = [Counter() for _ in names]
    n_pass_nonbnd = n_denovo = n_imprecise = 0
    n_pop = n_gene = 0

    for sv_id, c in enumerate(calls, start=1):
        qc_all_types[c.svtype] += 1
        for f in c.filter.split(";"):
            qc_filters[f] += 1
        if c.imprecise:
            n_imprecise += 1
        counted = c.passed and c.svtype != "BND"
        if c.passed:
            qc_pass_types[c.svtype] += 1
        if counted:
            n_pass_nonbnd += 1
            qc_sizes.setdefault(size_bin(c.size) or "n/a", Counter())[c.svtype] += 1
            for i, g in enumerate(c.samples):
                per_sample[i]["carrier" if g.copies > 0 else "ref" if g.copies == 0 else "nocall"] += 1

        f = c.samples[focus] if (focus is not None and c.samples) else None
        p = [c.samples[i] for i in parents] if c.samples else []
        inh = inheritance(f.copies, p[0].copies, p[1].copies) if (f and len(p) == 2) else None
        p_dp = [g.dp for g in p if g.dp is not None]
        p_min_dp = min(p_dp) if len(p) == 2 and len(p_dp) == 2 else None
        if inh == "de_novo" and counted:
            n_denovo += 1
        family = json.dumps([
            {"sample": names[i], "role": role_of[i], "gt": g.gt, "copies": g.copies, "ploidy": g.ploidy, "partial": g.partial, "gq": g.gq,
             "alt": g.ad_alt, "ref": g.ad_ref, "dp": g.dp, "cn": g.cn}
            for i, g in enumerate(c.samples)
        ])

        hits: List[Dict[str, object]] = []
        if genes is not None:
            mate = (c.mate_chrom, c.mate_pos) if c.mate_chrom and c.mate_pos else None
            hits = impacts(genes, c.svtype, c.chrom, c.pos, c.end, mate)
        best = hits[0] if hits else None
        n_gene += bool(hits)
        for h in hits:
            gene_rows.append((sv_id, h["symbol"], h["rank"], int(h["biotype"] == "protein_coding")))

        matches: Dict[str, Match] = {}
        if catalogs is not None and c.svtype != "BND":
            qtype = c.svtype
            if qtype == "CNV" and f is not None and f.cn is not None:
                qtype = "DEL" if f.cn < 2 else "DUP" if f.cn > 2 else "CNV"
            matches = catalogs.match(qtype, c.chrom, c.pos, c.end, c.size)
        n_pop += bool(matches)
        pop_max, g_af, g_pm, h_af = _pop_summary(matches)
        pop_json = json.dumps({k: m._asdict() for k, m in matches.items()})

        segdup = odd = None
        if tracks is not None and c.svtype != "BND":
            end = max(c.end, c.pos)
            if tracks.segdups is not None:
                segdup = round(tracks.segdups.fraction(c.chrom, c.pos, end), 3)
            odd_names = tracks.odd_regions(c.chrom, c.pos, end)
            odd = ",".join(odd_names) or None
            repeat_spans.append((sv_id, c.chrom, c.pos, c.pos + 1 if c.svtype == "INS" else end))

        rows.append((
            sv_id, c.vid, c.chrom, chrom_order(c.chrom), c.pos, c.end, c.svtype, c.size, size_bin(c.size), c.filter, int(c.passed), c.qual,
            int(c.imprecise), int(c.symbolic), c.mate_chrom, c.mate_pos, c.event,
            f.gt if f else None, f.copies if f else None, f.ploidy if f else None, f.gq if f else None, f.ad_alt if f else None, f.dp if f else None, f.cn if f else None,
            p[0].copies if len(p) == 2 else None, p[0].dp if len(p) == 2 else None, p[1].copies if len(p) == 2 else None, p[1].dp if len(p) == 2 else None,
            inh, p_min_dp, family,
            len(hits), ",".join(str(h["symbol"]) for h in hits[:8]) or None, str(best["relation"]) if best else "none", int(best["rank"]) if best else 0,
            json.dumps(hits[:MAX_GENES_JSON]),
            len(matches), pop_max, g_af, g_pm, h_af, pop_json, segdup, odd, None, None,
        ))
        if len(rows) >= BATCH:
            con.executemany(f"INSERT INTO svs VALUES ({','.join('?' * 46)})", rows)
            con.executemany("INSERT INTO sv_genes VALUES (?, ?, ?, ?)", gene_rows)
            rows, gene_rows = [], []
    if rows:
        con.executemany(f"INSERT INTO svs VALUES ({','.join('?' * 46)})", rows)
        con.executemany("INSERT INTO sv_genes VALUES (?, ?, ?, ?)", gene_rows)

    if tracks is not None and tracks.has_repeats and repeat_spans:
        progress.phase = "labelling repeats"
        labels = tracks.repeat_labels(repeat_spans)
        con.executemany("UPDATE svs SET repeat = ?, repeat_frac = ? WHERE id = ?", [(lab, frac, k) for k, (lab, frac) in labels.items()])
    progress.phase = "indexing"
    con.executescript(INDEXES)

    bins = [b for _, b in SIZE_BINS] + [LAST_BIN]
    if not any(qc_sizes.get("< 50", Counter()).values()):
        bins.remove("< 50")  # most callers stop at 50 bp: no empty bar for them
    qc = {
        "records": stats.records, "skipped": stats.skipped, "invalid_records": stats.invalid, "source": stats.source, "samples": names,
        "by_type_all": {t: qc_all_types[t] for t in SV_TYPES}, "by_type_pass": {t: qc_pass_types[t] for t in SV_TYPES},
        "filters": [{"filter": k, "count": v} for k, v in qc_filters.most_common()],
        "pass_non_bnd": n_pass_nonbnd, "imprecise": n_imprecise,
        "size_bins": [{"bin": b, **{t: qc_sizes.get(b, Counter())[t] for t in SV_TYPES if t != "BND"}} for b in bins],
        "per_sample": [{"sample": names[i], "role": role_of[i], **{k: per_sample[i][k] for k in ("carrier", "ref", "nocall")}} for i in range(len(names))],
        "roles": {"proband": names[focus] if focus is not None else None, "parents": [names[i] for i in parents], "source": roles["source"], "basis": roles["basis"]},
        "de_novo_candidates": n_denovo if len(parents) == 2 else None,
        "with_gene": n_gene, "with_population_match": n_pop,
    }
    st = os.stat(vcf_path)
    meta = {
        "builder_version": str(BUILDER_VERSION), "vcf_path": str(vcf_path), "vcf_stamp": f"{st.st_size}:{st.st_mtime_ns}",
        "genes_stamp": genes.stamp if genes is not None else "", "genes_n": str(genes.n_genes) if genes is not None else "0",
        "catalogs_stamp": catalogs.stamp if catalogs is not None else "",
        "catalog_counts": json.dumps(catalogs.counts if catalogs is not None else {}),
        "tracks_stamp": tracks.stamp if tracks is not None else "", "repeats": "1" if tracks is not None and tracks.has_repeats else "0",
        "declared": json.dumps(declared, sort_keys=True), "n_svs": str(len(calls)), "qc": json.dumps(qc), "built": time.strftime("%Y-%m-%d %H:%M:%S"),
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


def is_current(meta: Optional[Dict[str, str]], vcf_path: str, stamps: Tuple[str, str, str], declared: Dict[str, Any]) -> bool:
    if not meta or meta.get("builder_version") != str(BUILDER_VERSION):
        return False
    st = os.stat(vcf_path)
    return (
        meta.get("vcf_stamp") == f"{st.st_size}:{st.st_mtime_ns}"
        and (meta.get("genes_stamp", ""), meta.get("catalogs_stamp", ""), meta.get("tracks_stamp", "")) == stamps
        and meta.get("declared") == json.dumps(declared, sort_keys=True)
    )
