"""Sample registry, index lifecycle and queries for the small-variant suite.

Samples are single-sample small-variant VCFs (DeepVariant output, or a truth set). They are
found under ``$SMALLVAR_RESULTS_ROOT`` (default ``~/smallvar_data/samples``) as
``<root>/<sample>/<anything>.vcf[.gz]`` with an optional ``sample.json`` (label, sex,
sex_source, caller). A small public HG002 slice ships in the repo so the suite works on a
fresh checkout; real samples stay outside it.

Each VCF is read once into a SQLite index (index.py) under ``$SMALLVAR_CACHE``, on a
background thread because a whole genome takes minutes; the screens poll ``status`` and
mount when it is ready. Annotation resources (ClinVar, gnomAD, HPRC, gene bodies) are
optional: when one is missing its columns stay empty and the API says so, so a missing
resource reads as "not evaluated", never as "nothing found".
"""

from __future__ import annotations

import csv
import json
import os
import sqlite3
import threading
import traceback
from pathlib import Path
from typing import Any, Dict, List, NamedTuple, Optional, Tuple

from . import index as index_mod
from .clinvar import ClinvarTable
from .genes import GeneIndex
from .gnotate import Gnotate
from .index import PAR_X, Progress

REPO_ROOT = Path(__file__).resolve().parents[2]
FIXTURE_DIR = REPO_ROOT / "backend" / "fixtures" / "hg002_smallvar"
REFERENCE_DIR = Path(__file__).resolve().parent / "reference"

DATA_ROOT = Path(os.environ.get("SMALLVAR_DATA", Path.home() / "smallvar_data"))
RESULTS_ROOT = Path(os.environ.get("SMALLVAR_RESULTS_ROOT", DATA_ROOT / "samples"))
CACHE_DIR = Path(os.environ.get("SMALLVAR_CACHE", DATA_ROOT / "cache"))
CLINVAR_PATH = os.environ.get("SMALLVAR_CLINVAR", str(DATA_ROOT / "clinvar_plp.grch38.tsv.gz"))
_GNOTATE_DIR = DATA_ROOT / "pacbio_bundle" / "static_resources" / "GRCh38" / "slivar_gnotate"
GNOMAD_PATH = os.environ.get("SMALLVAR_GNOMAD", str(_GNOTATE_DIR / "gnomad.hg38.v3.custom.v1.zip"))
HPRC_PATH = os.environ.get("SMALLVAR_HPRC", str(_GNOTATE_DIR / "hprc.deepvariant.glnexus.hg38.v1.zip"))
GENES_BED = os.environ.get("SMALLVAR_GENES_BED", str(REPO_ROOT / "project" / "uploads" / "genes.hg38.bed"))

GNOMAD_LABEL = "gnomAD v3.0 genomes"
HPRC_LABEL = "HPRC release 1 (DeepVariant + GLnexus)"
MAX_PAGE = 200


class NotReady(Exception):
    def __init__(self, status: Dict[str, Any]):
        super().__init__(status.get("state", "not ready"))
        self.status = status


class Sample(NamedTuple):
    key: str
    vcf: str
    info: Dict[str, Any]


# ---------------------------------------------------------------- discovery
def _read_info(directory: Path) -> Dict[str, Any]:
    p = directory / "sample.json"
    try:
        return json.loads(p.read_text()) if p.is_file() else {}
    except (OSError, ValueError):
        return {}


def _find_vcf(directory: Path) -> Optional[Path]:
    candidates = sorted(directory.glob("*.vcf.gz")) + sorted(directory.glob("*.vcf"))
    return candidates[0] if candidates else None


def discover() -> List[Sample]:
    found: List[Sample] = []
    if FIXTURE_DIR.is_dir() and _find_vcf(FIXTURE_DIR):
        found.append(Sample("HG002_chr22", str(_find_vcf(FIXTURE_DIR)), _read_info(FIXTURE_DIR)))
    if RESULTS_ROOT.is_dir():
        for d in sorted(RESULTS_ROOT.iterdir()):
            vcf = _find_vcf(d) if d.is_dir() else None
            if vcf:
                found.append(Sample(d.name, str(vcf), _read_info(d)))
    return found


def sample_keys() -> List[str]:
    return [s.key for s in discover()]


def default_sample() -> str:
    keys = sample_keys()
    return keys[0] if keys else ""


def resolve(key: Optional[str]) -> Sample:
    """Validate a sample key (None -> default). Raises KeyError if unknown."""
    samples = discover()
    wanted = key or (samples[0].key if samples else "")
    for s in samples:
        if s.key == wanted:
            return s
    raise KeyError(wanted)


# ---------------------------------------------------------------- annotation resources
def file_stamp(path: str) -> str:
    """Identifies a resource file by name, size and mtime; a change triggers a rebuild."""
    p = Path(path)
    if not p.is_file():
        return ""
    st = p.stat()
    return f"{p.name}:{st.st_size}:{st.st_mtime_ns}"


_resource_lock = threading.Lock()
_resources: Dict[str, Tuple[str, Any]] = {}


def _cached(name: str, path: str, loader) -> Optional[Any]:
    stamp = file_stamp(path)
    if not stamp:
        return None
    with _resource_lock:
        hit = _resources.get(name)
        if hit is None or hit[0] != stamp:
            _resources[name] = (stamp, loader(path))
        return _resources[name][1]


def clinvar_table() -> Optional[ClinvarTable]:
    return _cached("clinvar", CLINVAR_PATH, ClinvarTable.load)


def gnomad_gnotate() -> Optional[Gnotate]:
    return _cached("gnomad", GNOMAD_PATH, Gnotate)


def hprc_gnotate() -> Optional[Gnotate]:
    return _cached("hprc", HPRC_PATH, Gnotate)


def gene_index() -> Optional[GeneIndex]:
    return _cached("genes", GENES_BED, GeneIndex.from_bed)


def _stamps() -> Tuple[str, str, str]:
    return file_stamp(CLINVAR_PATH), file_stamp(GNOMAD_PATH), file_stamp(HPRC_PATH)


# ---------------------------------------------------------------- index lifecycle
_state_lock = threading.Lock()
_threads: Dict[str, threading.Thread] = {}
_progress: Dict[str, Progress] = {}
_errors: Dict[str, str] = {}


def db_path(sample: Sample) -> str:
    return str(CACHE_DIR / f"{sample.key}.sqlite")


def _current_meta(sample: Sample) -> Optional[Dict[str, str]]:
    meta = index_mod.read_meta(db_path(sample))
    if meta and index_mod.is_current(meta, sample.vcf, *_stamps()):
        return meta
    return None


def status(sample: Sample) -> Dict[str, Any]:
    if _current_meta(sample) is not None:
        return {"state": "ready"}
    with _state_lock:
        t = _threads.get(sample.key)
        if t is not None and t.is_alive():
            return {"state": "building", "progress": _progress[sample.key].as_dict()}
        if sample.key in _errors:
            return {"state": "error", "error": _errors[sample.key]}
    return {"state": "stale"}


def _build(sample: Sample, progress: Progress) -> None:
    try:
        progress.phase = "loading annotation resources"
        index_mod.build_index(
            sample.vcf, db_path(sample), gene_index(), clinvar_table(), gnomad_gnotate(), hprc_gnotate(), progress
        )
        with _state_lock:
            _errors.pop(sample.key, None)
        _overview_cache.pop(sample.key, None)
    except Exception as e:  # reported through status(); the screens show it
        traceback.print_exc()
        with _state_lock:
            _errors[sample.key] = f"{type(e).__name__}: {e}"


def ensure_index(sample: Sample) -> Dict[str, Any]:
    """Start a background build if the index is missing or stale; return the status."""
    st = status(sample)
    if st["state"] in ("stale", "error"):
        with _state_lock:
            _errors.pop(sample.key, None)
            progress = Progress()
            _progress[sample.key] = progress
            t = threading.Thread(target=_build, args=(sample, progress), daemon=True, name=f"smallvar-{sample.key}")
            _threads[sample.key] = t
            t.start()
        return {"state": "building", "progress": progress.as_dict()}
    return st


def _require_ready(sample: Sample) -> Dict[str, str]:
    meta = _current_meta(sample)
    if meta is None:
        raise NotReady(ensure_index(sample))
    return meta


def _connect(sample: Sample) -> sqlite3.Connection:
    con = sqlite3.connect(f"file:{db_path(sample)}?mode=ro", uri=True, check_same_thread=False)
    con.row_factory = sqlite3.Row
    return con


# ---------------------------------------------------------------- metadata and resources
def _qc(meta: Dict[str, str]) -> Dict[str, Any]:
    return json.loads(meta["qc"])


def _sex(sample: Sample, qc: Dict[str, Any]) -> Dict[str, Optional[str]]:
    declared = sample.info.get("sex")
    if declared in ("male", "female"):
        return {"sex": declared, "source": sample.info.get("sex_source") or "stated in sample.json"}
    return {"sex": qc.get("inferred_sex"), "source": qc.get("inferred_sex_basis")}


def resources(meta: Optional[Dict[str, str]]) -> List[Dict[str, Any]]:
    """Which annotation resources fed this sample's index, and which could not."""
    def installed(path: str) -> bool:
        return Path(path).is_file()

    out: List[Dict[str, Any]] = []
    cv_ok = bool(meta and meta.get("clinvar_stamp"))
    out.append({
        "key": "clinvar", "label": "ClinVar", "installed": cv_ok,
        "detail": (f"release {meta['clinvar_release']}, {int(meta['clinvar_rows']):,} GRCh38 small variants classified "
                   "pathogenic, likely pathogenic or conflicting; a variant outside this subset is not shown to be absent "
                   "from ClinVar") if cv_ok else "not installed: run scripts/real_data/build_smallvar_clinvar.py",
    })
    gn_ok = bool(meta and meta.get("gnomad_stamp"))
    gn_info = json.loads(meta["gnomad_info"]) if gn_ok and meta.get("gnomad_info") else {}
    out.append({
        "key": "gnomad", "label": GNOMAD_LABEL, "installed": gn_ok,
        "detail": (f"PASS sites only, overall allele frequency, from PacBio's WDL resource bundle ({gn_info.get('file')}, "
                   f"built {gn_info.get('created', '').replace('created on:', '').strip()}); absent means not seen at a PASS site "
                   "in that release") if gn_ok else "not installed: point SMALLVAR_GNOMAD at a slivar gnotate zip",
    })
    hp_ok = bool(meta and meta.get("hprc_stamp"))
    out.append({
        "key": "hprc", "label": HPRC_LABEL, "installed": hp_ok,
        "detail": "long-read cohort frequencies from PacBio's resource bundle: about 100 samples, so useful for spotting recurrent "
                  "long-read calls, not for ruling variants out. The cohort includes GIAB samples (HG001, HG003, HG004), so the "
                  "variants of related samples such as HG002 are over-represented" if hp_ok else "not installed",
    })
    genes_n = int(meta["genes_loaded"]) if meta else 0
    out.append({
        "key": "genes", "label": "Gene bodies", "installed": genes_n > 0,
        "detail": f"{genes_n:,} RefSeq gene loci; overlap with a gene body only, not coding consequence" if genes_n else "not installed",
    })
    for key, label, detail in (
        ("consequence", "Coding consequence", "needs a transcript-aware annotator such as SnpEff or VEP"),
        ("insilico", "In-silico scores", "CADD, REVEL, SpliceAI and AlphaMissense are not installed"),
        ("hpo", "Phenotype (HPO)", "no phenotype terms are loaded for this sample"),
        ("parents", "Parental genotypes", "single sample; no trio or duo"),
    ):
        out.append({"key": key, "label": label, "installed": False, "detail": f"not evaluated: {detail}"})
    return out


# The deck's map from ACMG/AMP criteria to the evidence behind them; shown so a reader can
# see which criteria this tool has evidence for.
def evidence_map(res: List[Dict[str, Any]], phased_fraction: Optional[float]) -> List[Dict[str, Any]]:
    have = {r["key"]: r["installed"] for r in res}
    return [
        {"criteria": "BA1 · BS1 · PM2", "evidence": "Population frequency", "available": have["gnomad"],
         "note": "gnomAD v3.0 overall frequency" if have["gnomad"] else "not evaluated"},
        {"criteria": "PP3 · BP4 · BP7", "evidence": "In-silico predictors", "available": have["insilico"], "note": "not evaluated"},
        {"criteria": "PVS1 · PM4 · BP3", "evidence": "Coding consequence", "available": have["consequence"], "note": "not evaluated"},
        {"criteria": "PS1 · PM5", "evidence": "ClinVar", "available": have["clinvar"],
         "note": "exact-allele match only; no same-residue (PS1/PM5) comparison" if have["clinvar"] else "not evaluated"},
        {"criteria": "PM3 · PS2 · PM6", "evidence": "Phasing and parental genotypes", "available": False,
         "note": (f"{100 * phased_fraction:.0f}% of heterozygous calls are phased; no parental genotypes"
                  if phased_fraction else "not evaluated: calls are unphased and there are no parental genotypes")},
        {"criteria": "PP4", "evidence": "Phenotype match (HPO)", "available": have["hpo"], "note": "not evaluated"},
    ]


# ---------------------------------------------------------------- queries
SORTS = {
    "position": "chrom_order {d}, pos {d}",
    "af": "gn_af IS NULL, gn_af {d}, chrom_order, pos",
    "stars": "cv_stars IS NULL, cv_stars {d}, chrom_order, pos",
    "gq": "gq IS NULL, gq {d}, chrom_order, pos",
}


class Filters(NamedTuple):
    pass_only: bool = True
    kinds: Tuple[str, ...] = ()
    zygosity: str = "any"  # any | het | hom
    chrom: str = ""
    in_gene: bool = False  # overlaps a gene body
    genes: Tuple[str, ...] = ()
    clinvar: str = "any"  # any | plp | plp_conflicting | has_match | no_match
    min_stars: int = 0
    max_af: Optional[float] = None  # gnomAD overall AF ceiling
    af_missing: str = "include"  # include | exclude: variants absent from gnomAD
    max_hprc_af: Optional[float] = None
    min_gq: Optional[int] = None
    min_dp: Optional[int] = None


def _clauses(f: Filters, have_gnomad: bool, have_hprc: bool, have_clinvar: bool) -> List[Tuple[str, str, List[Any]]]:
    """(label, SQL, params) per active filter, in the order the funnel reports them.

    A filter whose resource is not installed is left out, not applied: with no ClinVar column
    filled in, 'ClinVar pathogenic' would match nothing and read as a negative finding."""
    c: List[Tuple[str, str, List[Any]]] = []
    if f.pass_only:
        c.append(("passes the VCF FILTER", "passed = 1", []))
    if f.kinds:
        c.append((f"type is {'/'.join(f.kinds)}", f"kind IN ({','.join('?' * len(f.kinds))})", list(f.kinds)))
    if f.zygosity == "het":
        c.append(("heterozygous", "copies = 1", []))
    elif f.zygosity == "hom":
        c.append(("homozygous", "copies >= 2", []))
    if f.min_gq is not None:
        c.append((f"GQ ≥ {f.min_gq}", "(gq IS NULL OR gq >= ?)", [f.min_gq]))
    if f.min_dp is not None:
        c.append((f"depth ≥ {f.min_dp}", "(dp IS NULL OR dp >= ?)", [f.min_dp]))
    if f.chrom:
        c.append((f"on {f.chrom}", "chrom = ?", [f.chrom]))
    if f.in_gene:
        c.append(("in a gene body", "id IN (SELECT variant_id FROM variant_genes)", []))
    if f.genes:
        c.append((
            f"in {', '.join(f.genes[:4])}{'…' if len(f.genes) > 4 else ''}",
            f"id IN (SELECT variant_id FROM variant_genes WHERE gene IN ({','.join('?' * len(f.genes))}))",
            list(f.genes),
        ))
    cv = {
        "plp": ("ClinVar pathogenic / likely pathogenic", "cv_group = 'plp'"),
        "plp_conflicting": ("ClinVar P/LP or conflicting", "cv_group IN ('plp','conflicting')"),
        "has_match": ("any ClinVar match", "cv_group IS NOT NULL"),
        "no_match": ("no ClinVar match", "cv_group IS NULL"),
    }.get(f.clinvar)
    if cv and have_clinvar:
        c.append((cv[0], cv[1], []))
    if f.min_stars and have_clinvar:
        c.append((f"ClinVar review ≥ {f.min_stars} star{'s' if f.min_stars > 1 else ''}", "cv_stars >= ?", [f.min_stars]))
    if f.max_af is not None and have_gnomad:
        absent = "gn_af IS NULL OR " if f.af_missing == "include" else ""
        c.append((f"gnomAD AF ≤ {f.max_af:g}", f"({absent}gn_af <= ?)", [f.max_af]))
    if f.max_hprc_af is not None and have_hprc:
        c.append((f"HPRC AF ≤ {f.max_hprc_af:g}", "(hp_af IS NULL OR hp_af <= ?)", [f.max_hprc_af]))
    return c


def _where(clauses) -> Tuple[str, List[Any]]:
    if not clauses:
        return "", []
    return " WHERE " + " AND ".join(c[1] for c in clauses), [p for c in clauses for p in c[2]]


def _zygosity_label(row: sqlite3.Row, sex: Optional[str]) -> str:
    if row["chrom"] == "chrX" and sex == "male" and not any(lo <= row["pos"] <= hi for lo, hi in PAR_X):
        return "hemizygous"
    return "het" if row["copies"] == 1 else "hom"


def _row(row: sqlite3.Row, sex: Optional[str]) -> Dict[str, Any]:
    def has(prefix: str) -> bool:
        return row[f"{prefix}_af"] is not None

    return {
        "id": row["id"], "chrom": row["chrom"], "pos": row["pos"], "ref": row["ref"], "alt": row["alt"], "kind": row["kind"],
        "genes": row["genes"].split(",") if row["genes"] else [],
        "genotype": {
            "gt": row["gt"], "copies": row["copies"], "zygosity": _zygosity_label(row, sex), "phased": bool(row["phased"]),
            "phase_set": row["ps"], "haplotypes": row["hap"],
        },
        "quality": {"qual": row["qual"], "filter": row["filter"], "passed": bool(row["passed"]), "dp": row["dp"], "gq": row["gq"], "vaf": row["vaf"]},
        "clinvar": None if row["cv_group"] is None else {
            "classification": row["cv_class"], "group": row["cv_group"], "stars": row["cv_stars"], "review": row["cv_review"],
            "variation_id": row["cv_id"], "name": row["cv_name"], "conditions": row["cv_cond"],
            "last_evaluated": row["cv_evaluated"], "n_submitters": row["cv_nsub"], "genes": row["cv_gene"],
        },
        "gnomad": {"af": row["gn_af"], "ac": row["gn_ac"], "hom": row["gn_hom"]} if has("gn") else None,
        "hprc": {"af": row["hp_af"], "ac": row["hp_ac"], "hom": row["hp_hom"]} if has("hp") else None,
    }


def query_variants(
    sample: Sample, f: Filters, sort: str = "position", order: str = "asc", page: int = 1, page_size: int = 50,
    explain: bool = False,
) -> Dict[str, Any]:
    meta = _require_ready(sample)
    have_gnomad, have_hprc = bool(meta.get("gnomad_stamp")), bool(meta.get("hprc_stamp"))
    have_clinvar = bool(meta.get("clinvar_stamp"))
    sex = _sex(sample, _qc(meta))["sex"]
    clauses = _clauses(f, have_gnomad, have_hprc, have_clinvar)
    where, params = _where(clauses)
    page_size = max(1, min(page_size, MAX_PAGE))
    page = max(1, page)
    sort_sql = SORTS.get(sort, SORTS["position"]).format(d="DESC" if order == "desc" else "ASC")
    ignored = []
    if f.max_af is not None and not have_gnomad:
        ignored.append("max_af")
    if f.max_hprc_af is not None and not have_hprc:
        ignored.append("max_hprc_af")
    if f.clinvar != "any" and not have_clinvar:
        ignored.append("clinvar")
    if f.min_stars and not have_clinvar:
        ignored.append("min_stars")
    con = _connect(sample)
    try:
        total_all = int(meta["n_variants"])
        matching = con.execute(f"SELECT COUNT(*) FROM variants{where}", params).fetchone()[0]
        rows = con.execute(
            f"SELECT * FROM variants{where} ORDER BY {sort_sql} LIMIT ? OFFSET ?", params + [page_size, (page - 1) * page_size]
        ).fetchall()
        funnel: Optional[List[Dict[str, Any]]] = None
        if explain or matching == 0:
            funnel = [{"label": "all calls in the VCF", "count": total_all}]
            for i in range(len(clauses)):
                w, p = _where(clauses[: i + 1])
                funnel.append({"label": clauses[i][0], "count": con.execute(f"SELECT COUNT(*) FROM variants{w}", p).fetchone()[0]})
    finally:
        con.close()
    return {
        "total_all": total_all, "total_matching": matching, "page": page, "page_size": page_size,
        "rows": [_row(r, sex) for r in rows], "funnel": funnel, "ignored_filters": ignored,
        "gnomad_installed": have_gnomad, "hprc_installed": have_hprc, "clinvar_installed": have_clinvar,
    }


def get_variant(sample: Sample, variant_id: int) -> Optional[Dict[str, Any]]:
    meta = _require_ready(sample)
    con = _connect(sample)
    try:
        row = con.execute("SELECT * FROM variants WHERE id = ?", (variant_id,)).fetchone()
    finally:
        con.close()
    return _row(row, _sex(sample, _qc(meta))["sex"]) if row is not None else None


_overview_cache: Dict[str, Dict[str, Any]] = {}


def overview(sample: Sample) -> Dict[str, Any]:
    meta = _require_ready(sample)
    cached = _overview_cache.get(sample.key)
    if cached is not None and cached["_built"] == meta["built"]:
        return cached
    qc = _qc(meta)
    have_gnomad, have_hprc = bool(meta.get("gnomad_stamp")), bool(meta.get("hprc_stamp"))
    con = _connect(sample)
    try:
        def count(where: str = "", params: Tuple = ()) -> int:
            return con.execute(f"SELECT COUNT(*) FROM variants{where}", params).fetchone()[0]

        funnel = [{"label": "All calls in the VCF", "count": int(meta["n_variants"])}, {"label": "Pass the VCF filter", "count": count(" WHERE passed = 1")}]
        steps = ["passed = 1"]
        if have_gnomad:
            steps.append("(gn_af IS NULL OR gn_af <= 0.01)")
            funnel.append({"label": "Rare: gnomAD AF ≤ 1% or absent", "count": count(" WHERE " + " AND ".join(steps))})
        steps.append("id IN (SELECT variant_id FROM variant_genes)")
        funnel.append({"label": "In a gene body", "count": count(" WHERE " + " AND ".join(steps))})
        if meta.get("clinvar_stamp"):
            plp = " AND ".join(steps + ["cv_group = 'plp'"])
            funnel.append({"label": "ClinVar pathogenic / likely pathogenic", "count": count(" WHERE " + plp)})
        top = []
        if meta.get("clinvar_stamp"):
            top = [
                {"gene": r[0], "count": r[1]}
                for r in con.execute(
                    "SELECT g.gene, COUNT(*) FROM variants v JOIN variant_genes g ON g.variant_id = v.id "
                    "WHERE v.passed = 1 AND v.cv_group = 'plp' GROUP BY g.gene ORDER BY COUNT(*) DESC, g.gene LIMIT 12"
                )
            ]
        clinvar = {
            "plp": count(" WHERE passed = 1 AND cv_group = 'plp'"),
            "conflicting": count(" WHERE passed = 1 AND cv_group = 'conflicting'"),
        } if meta.get("clinvar_stamp") else None
        plp_rare = None
        if meta.get("clinvar_stamp") and have_gnomad:
            plp_rare = count(" WHERE passed = 1 AND cv_group = 'plp' AND (gn_af IS NULL OR gn_af <= 0.01)")
        gnomad_cov = None
        if have_gnomad:
            snv_like = count(" WHERE passed = 1")
            in_g = count(" WHERE passed = 1 AND gn_af IS NOT NULL")
            common = count(" WHERE passed = 1 AND gn_af > 0.01")
            gnomad_cov = {"found": in_g, "of": snv_like, "common_over_1pct": common}
    finally:
        con.close()
    res = resources(meta)
    out = {
        "_built": meta["built"], "sample": sample.key, "qc": qc, "funnel": funnel, "clinvar": clinvar,
        "clinvar_plp_rare": plp_rare, "top_plp_genes": top, "gnomad_coverage": gnomad_cov,
        "resources": res, "evidence": evidence_map(res, qc.get("phased_het_fraction")),
        "sex": _sex(sample, qc),
    }
    _overview_cache[sample.key] = out
    return out


def meta_for(sample: Sample) -> Dict[str, Any]:
    st = status(sample)
    base = {
        "sample": sample.key, "label": sample.info.get("label") or sample.key, "reference": "GRCh38",
        "caller": sample.info.get("caller"), "vcf": os.path.basename(sample.vcf), "status": st,
    }
    if st["state"] == "ready":
        meta = _current_meta(sample)
        qc = _qc(meta) if meta else {}
        base.update({"built": meta["built"] if meta else None, "n_variants": int(meta["n_variants"]) if meta else None,
                     "sex": _sex(sample, qc) if qc else None})
    return base


def gene_presets() -> List[Dict[str, Any]]:
    path = REFERENCE_DIR / "acmg_sf_v3_3.tsv"
    with open(path, encoding="utf-8", newline="") as f:
        genes = [r["gene"] for r in csv.DictReader(f, delimiter="\t")]
    return [{"key": "acmg_sf_v3_3", "label": "ACMG SF v3.3 genes", "genes": genes,
             "note": "A convenience list for filtering; this suite does not apply the ACMG SF reporting rules."}]
