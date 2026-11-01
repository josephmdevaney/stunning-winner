"""Call-set registry, index lifecycle and queries for the structural-variant suite.

A call set is one SV VCF (a joint call across a family, or a single sample) under
``$STRUCTVAR_RESULTS_ROOT`` (default ``~/sv_data/samples``) as ``<root>/<callset>/<name>.vcf[.gz]``, with
an optional ``sample.json``: ``{"label", "caller", "proband", "parents": [..]}``. Real families stay outside
the repo; the suite ships no sample data.

Each VCF is read once into a SQLite index (index.py) on a background thread; the screens poll ``status``
and mount when it is ready. The annotation resources (gene models, population catalogs, genome tracks)
are optional: a missing one leaves its columns empty, the API says so, and a filter that depends on it is
left out rather than applied, so a missing resource reads as "not evaluated", never as "nothing found".
"""

from __future__ import annotations

import json
import os
import sqlite3
import threading
import traceback
import zlib
from pathlib import Path
from typing import Any, Dict, List, NamedTuple, Optional, Tuple

from . import index as index_mod
from .catalogs import SOURCE_LABEL, SOURCES, Catalogs
from .genes import IMPACT_LABEL, GeneModels
from .index import Progress
from .normalize import SV_TYPES, canon_chrom
from .tracks import Tracks

DATA_ROOT = Path(os.environ.get("STRUCTVAR_DATA", Path.home() / "sv_data"))
RESULTS_ROOT = Path(os.environ.get("STRUCTVAR_RESULTS_ROOT", DATA_ROOT / "samples"))
CACHE_DIR = Path(os.environ.get("STRUCTVAR_CACHE", DATA_ROOT / "cache"))
GENES_PATH = os.environ.get("STRUCTVAR_GENES", str(DATA_ROOT / "sv_genes.grch38.tsv.gz"))
CATALOGS_PATH = os.environ.get("STRUCTVAR_CATALOGS", str(DATA_ROOT / "sv_catalogs.grch38.tsv.gz"))
BUNDLE_PATH = os.environ.get("STRUCTVAR_BUNDLE", str(DATA_ROOT / "pacbio_bundle" / "static_resources"))
MAX_PAGE = 200
IMPACT_FILTER_RANK = {"any": 0, "gene": 2, "exon": 3, "coding": 4, "whole_gene": 5}


class NotReady(Exception):
    def __init__(self, status: Dict[str, Any]):
        super().__init__(status.get("state", "not ready"))
        self.status = status


class Callset(NamedTuple):
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
    found = sorted(directory.glob("*.vcf.gz")) + sorted(directory.glob("*.vcf"))
    return found[0] if found else None


def discover() -> List[Callset]:
    out: List[Callset] = []
    if RESULTS_ROOT.is_dir():
        for d in sorted(RESULTS_ROOT.iterdir()):
            vcf = _find_vcf(d) if d.is_dir() else None
            if vcf:
                out.append(Callset(d.name, str(vcf), _read_info(d)))
    return out


def default_sample() -> str:
    found = discover()
    return found[0].key if found else ""


def resolve(key: Optional[str]) -> Callset:
    found = discover()
    wanted = key or (found[0].key if found else "")
    for c in found:
        if c.key == wanted:
            return c
    raise KeyError(wanted)


def _declared(c: Callset) -> Dict[str, Any]:
    return {k: c.info[k] for k in ("proband", "parents") if k in c.info}


# ---------------------------------------------------------------- annotation resources
def file_stamp(path: str) -> str:
    p = Path(path)
    if not p.is_file():
        return ""
    st = p.stat()
    return f"{p.name}:{st.st_size}:{st.st_mtime_ns}"


_resource_lock = threading.Lock()
_resources: Dict[str, Tuple[str, Any]] = {}


def _cached(name: str, stamp: str, loader) -> Optional[Any]:
    if not stamp:
        return None
    with _resource_lock:
        hit = _resources.get(name)
        if hit is None or hit[0] != stamp:
            _resources[name] = (stamp, loader())
        return _resources[name][1]


def gene_models() -> Optional[GeneModels]:
    return _cached("genes", file_stamp(GENES_PATH), lambda: GeneModels.load(GENES_PATH))


def catalogs() -> Optional[Catalogs]:
    return _cached("catalogs", file_stamp(CATALOGS_PATH), lambda: Catalogs.load(CATALOGS_PATH))


def tracks() -> Optional[Tracks]:
    stamp = Tracks.stamp_for(BUNDLE_PATH)
    return _cached("tracks", stamp, lambda: Tracks(BUNDLE_PATH))


def _stamps() -> Tuple[str, str, str]:
    return file_stamp(GENES_PATH), file_stamp(CATALOGS_PATH), Tracks.stamp_for(BUNDLE_PATH)


# ---------------------------------------------------------------- index lifecycle
_state_lock = threading.Lock()
_threads: Dict[str, threading.Thread] = {}
_progress: Dict[str, Progress] = {}
_errors: Dict[str, str] = {}


def db_path(c: Callset) -> str:
    return str(CACHE_DIR / f"{c.key}.sqlite")


def _current_meta(c: Callset) -> Optional[Dict[str, str]]:
    meta = index_mod.read_meta(db_path(c))
    return meta if meta and index_mod.is_current(meta, c.vcf, _stamps(), _declared(c)) else None


def status(c: Callset) -> Dict[str, Any]:
    if _current_meta(c) is not None:
        return {"state": "ready"}
    with _state_lock:
        t = _threads.get(c.key)
        if t is not None and t.is_alive():
            return {"state": "building", "progress": _progress[c.key].as_dict()}
        if c.key in _errors:
            return {"state": "error", "error": _errors[c.key]}
    return {"state": "stale"}


def _explain(e: Exception) -> str:
    if isinstance(e, (EOFError, zlib.error)) or type(e).__name__ == "BadGzipFile":
        return f"The VCF is truncated or corrupt ({e}). Copy it again and check it with `gzip -t`."
    return f"{type(e).__name__}: {e}"


def _build(c: Callset, progress: Progress) -> None:
    try:
        progress.phase = "loading annotation resources"
        index_mod.build_index(c.vcf, db_path(c), gene_models(), catalogs(), tracks(), _declared(c), progress)
        with _state_lock:
            _errors.pop(c.key, None)
        _overview_cache.pop(c.key, None)
    except Exception as e:  # reported through status(); the screens show it
        traceback.print_exc()
        with _state_lock:
            _errors[c.key] = _explain(e)


def ensure_index(c: Callset) -> Dict[str, Any]:
    st = status(c)
    if st["state"] in ("stale", "error"):
        with _state_lock:
            _errors.pop(c.key, None)
            progress = Progress()
            _progress[c.key] = progress
            t = threading.Thread(target=_build, args=(c, progress), daemon=True, name=f"structvar-{c.key}")
            _threads[c.key] = t
            t.start()
        return {"state": "building", "progress": progress.as_dict()}
    return st


def _require_ready(c: Callset) -> Dict[str, str]:
    meta = _current_meta(c)
    if meta is None:
        raise NotReady(ensure_index(c))
    return meta


def _connect(c: Callset) -> sqlite3.Connection:
    con = sqlite3.connect(f"file:{db_path(c)}?mode=ro", uri=True, check_same_thread=False)
    con.row_factory = sqlite3.Row
    return con


# ---------------------------------------------------------------- resources and evidence, for the overview
def resources(meta: Dict[str, str]) -> List[Dict[str, Any]]:
    counts = json.loads(meta.get("catalog_counts") or "{}")
    genes_ok, cat_ok, tr_ok = bool(meta.get("genes_stamp")), bool(meta.get("catalogs_stamp")), bool(meta.get("tracks_stamp"))
    out = [{
        "key": "genes", "label": "Gene models", "installed": genes_ok,
        "detail": (f"{int(meta['genes_n']):,} genes; MANE Select transcripts (UCSC mane.bb, NCBI/EMBL-EBI) where a gene has one, "
                   "otherwise the Ensembl 101 transcript with the longest coding sequence. Exon overlap with the transcript, not a "
                   "predicted consequence") if genes_ok else "not installed: run scripts/real_data/build_sv_genes.py",
    }]
    for s in SOURCES:
        what = {"gnomad_sv": "allele frequency overall and per ancestry (via dbVar nstd166)", "hprc": "allele frequency from about 100 long-read genomes",
                "decode": "presence only (no frequency in the file)", "eee": "number of discovery genomes carrying it (no frequency)"}[s]
        out.append({"key": s, "label": SOURCE_LABEL[s], "installed": cat_ok and s in counts,
                    "detail": f"{counts[s]:,} calls; {what}" if cat_ok and s in counts else "not installed: run scripts/real_data/build_sv_catalogs.py"})
    out.append({
        "key": "tracks", "label": "Segmental duplications, repeats, gaps", "installed": tr_ok,
        "detail": ("UCSC tracks from PacBio's resource bundle: segmental duplications, RepeatMasker and tandem repeats (SVs up to 100 kb"
                   + ("" if meta.get("repeats") == "1" else "; repeats not evaluated: tabix or the repeat file is missing")
                   + "), and gap / telomere / centromere / PAR / V(D)J regions") if tr_ok else "not installed",
    })
    for key, label, detail in (
        ("clingen", "ClinGen dosage sensitivity", "needed for the ACMG/ClinGen CNV scoring; not loaded"),
        ("insilico", "SV deleteriousness (CADD-SV)", "not loaded"),
        ("needlr", "needLR population frequency", "not loaded"),
    ):
        out.append({"key": key, "label": label, "installed": False, "detail": f"not evaluated: {detail}"})
    return out


# ---------------------------------------------------------------- queries
SORTS = {
    "position": "chrom_order {d}, pos {d}",
    "size": "size IS NULL, size {d}, chrom_order, pos",
    "impact": "impact_rank {d}, chrom_order, pos",
    "af": "pop_af_max IS NULL, pop_af_max {d}, chrom_order, pos",
    "gq": "f_gq IS NULL, f_gq {d}, chrom_order, pos",
}


class Filters(NamedTuple):
    pass_only: bool = True
    include_bnd: bool = False
    types: Tuple[str, ...] = ()
    min_size: Optional[int] = None
    max_size: Optional[int] = None
    carried: bool = True  # the focus sample carries it (the VCF holds the whole family's calls)
    zygosity: str = "any"  # any | het | hom
    inheritance: str = "any"  # any | de_novo | inherited | unknown
    min_parent_dp: Optional[int] = None
    min_gq: Optional[int] = None
    min_alt: Optional[int] = None
    impact: str = "any"  # any | gene | exon | coding | whole_gene
    genes: Tuple[str, ...] = ()
    chrom: str = ""
    max_pop_af: Optional[float] = None
    novel: bool = False  # in no population catalog
    max_segdup: Optional[float] = None
    no_odd: bool = False
    no_repeat: bool = False
    precise_only: bool = False


def _clauses(f: Filters, have_pop: bool, have_genes: bool, have_tracks: bool, trio: bool) -> List[Tuple[str, str, List[Any]]]:
    """(label, SQL, params) per active filter, in the order the funnel reports them. A filter whose
    resource is not installed is left out: with nothing annotated, 'rare' would match everything and
    'in an exon' nothing, and either would read as a finding."""
    c: List[Tuple[str, str, List[Any]]] = []
    if f.pass_only:
        c.append(("passes the caller's FILTER", "passed = 1", []))
    if not f.include_bnd:
        c.append(("not a breakend", "svtype != 'BND'", []))
    if f.types:
        c.append((f"type is {'/'.join(f.types)}", f"svtype IN ({','.join('?' * len(f.types))})", list(f.types)))
    if f.min_size is not None:
        c.append((f"size ≥ {f.min_size:,} bp", "(size IS NULL OR size >= ?)", [f.min_size]))
    if f.max_size is not None:
        c.append((f"size ≤ {f.max_size:,} bp", "(size IS NULL OR size <= ?)", [f.max_size]))
    if f.carried:
        c.append(("carried by the proband", "f_copies > 0", []))
    if f.zygosity == "het":
        c.append(("heterozygous", "f_copies = 1", []))
    elif f.zygosity == "hom":
        c.append(("homozygous", "f_copies >= 2", []))
    if trio and f.inheritance in ("de_novo", "inherited", "unknown"):
        label = {"de_novo": "de novo candidate", "inherited": "inherited", "unknown": "inheritance unknown"}[f.inheritance]
        c.append((label, "inh = ?", [f.inheritance]))
    if trio and f.min_parent_dp is not None:
        c.append((f"parents' depth ≥ {f.min_parent_dp}", "p_min_dp >= ?", [f.min_parent_dp]))
    if f.min_gq is not None:
        c.append((f"GQ ≥ {f.min_gq}", "(f_gq IS NULL OR f_gq >= ?)", [f.min_gq]))
    if f.min_alt is not None:
        c.append((f"≥ {f.min_alt} supporting reads", "(f_alt IS NULL OR f_alt >= ?)", [f.min_alt]))
    if f.chrom:
        c.append((f"on {f.chrom}", "chrom = ?", [f.chrom]))
    if have_genes and f.impact in IMPACT_FILTER_RANK and f.impact != "any":
        c.append((f"impact: {f.impact.replace('_', ' ')}", "impact_rank >= ?", [IMPACT_FILTER_RANK[f.impact]]))
    if have_genes and f.genes:
        c.append((f"in {', '.join(f.genes[:4])}{'…' if len(f.genes) > 4 else ''}",
                  f"id IN (SELECT sv_id FROM sv_genes WHERE symbol IN ({','.join('?' * len(f.genes))}))", list(f.genes)))
    if have_pop and f.max_pop_af is not None:
        c.append((f"population AF ≤ {f.max_pop_af:g}", "(pop_af_max IS NULL OR pop_af_max <= ?)", [f.max_pop_af]))
    if have_pop and f.novel:
        c.append(("in no population catalog", "pop_n = 0", []))
    if have_tracks and f.max_segdup is not None:
        c.append((f"segmental duplication < {f.max_segdup:.0%}", "(segdup_frac IS NULL OR segdup_frac < ?)", [f.max_segdup]))
    if have_tracks and f.no_odd:
        c.append(("not in a gap / telomere / PAR / V(D)J region", "odd IS NULL", []))
    if have_tracks and f.no_repeat:
        c.append(("not mostly repeat", "(repeat IS NULL OR repeat_frac < 0.5)", []))
    if f.precise_only:
        c.append(("precise breakpoints", "imprecise = 0", []))
    return c


def _where(clauses) -> Tuple[str, List[Any]]:
    if not clauses:
        return "", []
    return " WHERE " + " AND ".join(c[1] for c in clauses), [p for c in clauses for p in c[2]]


def _row(r: sqlite3.Row) -> Dict[str, Any]:
    return {
        "id": r["id"], "vid": r["vid"], "chrom": r["chrom"], "pos": r["pos"], "end": r["end"], "svtype": r["svtype"], "size": r["size"],
        "filter": r["filter"], "passed": bool(r["passed"]), "qual": r["qual"], "imprecise": bool(r["imprecise"]), "symbolic": bool(r["symbolic"]),
        "mate": {"chrom": r["mate_chrom"], "pos": r["mate_pos"]} if r["mate_chrom"] else None, "event": r["event"],
        "focus": {"gt": r["f_gt"], "copies": r["f_copies"], "gq": r["f_gq"], "alt": r["f_alt"], "dp": r["f_dp"], "cn": r["f_cn"]},
        "inheritance": r["inh"], "parent_min_dp": r["p_min_dp"], "family": json.loads(r["family"]) if r["family"] else [],
        "n_genes": r["n_genes"], "genes": json.loads(r["gene_json"]) if r["gene_json"] else [], "impact": r["impact"],
        "population": json.loads(r["pop_json"]) if r["pop_json"] else {}, "pop_af_max": r["pop_af_max"], "pop_n": r["pop_n"],
        "segdup_frac": r["segdup_frac"], "odd": r["odd"].split(",") if r["odd"] else [], "repeat": r["repeat"], "repeat_frac": r["repeat_frac"],
    }


def query_svs(c: Callset, f: Filters, sort: str = "position", order: str = "asc", page: int = 1, page_size: int = 50,
              explain: bool = False) -> Dict[str, Any]:
    meta = _require_ready(c)
    have_pop, have_genes, have_tracks = bool(meta.get("catalogs_stamp")), bool(meta.get("genes_stamp")), bool(meta.get("tracks_stamp"))
    qc = json.loads(meta["qc"])
    trio = len(qc["roles"]["parents"]) == 2
    clauses = _clauses(f, have_pop, have_genes, have_tracks, trio)
    where, params = _where(clauses)
    page_size = max(1, min(page_size, MAX_PAGE))
    page = max(1, page)
    sort_sql = SORTS.get(sort, SORTS["position"]).format(d="DESC" if order == "desc" else "ASC")
    ignored: List[str] = []
    if f.max_pop_af is not None and not have_pop:
        ignored.append("max_pop_af")
    if f.novel and not have_pop:
        ignored.append("novel")
    if (f.impact != "any" or f.genes) and not have_genes:
        ignored.append("impact" if f.impact != "any" else "genes")
    if (f.max_segdup is not None or f.no_odd or f.no_repeat) and not have_tracks:
        ignored.append("tracks")
    if f.inheritance != "any" and not trio:
        ignored.append("inheritance")
    if f.min_parent_dp is not None and not trio:
        ignored.append("min_parent_dp")
    con = _connect(c)
    try:
        total_all = int(meta["n_svs"])
        matching = con.execute(f"SELECT COUNT(*) FROM svs{where}", params).fetchone()[0]
        rows = con.execute(f"SELECT * FROM svs{where} ORDER BY {sort_sql} LIMIT ? OFFSET ?", params + [page_size, (page - 1) * page_size]).fetchall()
        funnel = None
        if explain or matching == 0:
            funnel = [{"label": "all calls in the VCF", "count": total_all}]
            for i in range(len(clauses)):
                w, p = _where(clauses[: i + 1])
                funnel.append({"label": clauses[i][0], "count": con.execute(f"SELECT COUNT(*) FROM svs{w}", p).fetchone()[0]})
    finally:
        con.close()
    return {
        "total_all": total_all, "total_matching": matching, "page": page, "page_size": page_size, "rows": [_row(r) for r in rows],
        "funnel": funnel, "ignored_filters": ignored, "population_installed": have_pop, "genes_installed": have_genes,
        "tracks_installed": have_tracks, "trio": trio,
    }


def get_sv(c: Callset, sv_id: int) -> Optional[Dict[str, Any]]:
    _require_ready(c)
    con = _connect(c)
    try:
        r = con.execute("SELECT * FROM svs WHERE id = ?", (sv_id,)).fetchone()
    finally:
        con.close()
    return _row(r) if r is not None else None


def region(c: Callset, chrom: str, start: int, end: int, limit: int = 300) -> Dict[str, Any]:
    """Genes (MANE exons) and the proband's SVs in a window, for the locus view."""
    meta = _require_ready(c)
    chrom = canon_chrom(chrom)
    start, end = max(1, start), max(start, end)
    models = gene_models() if meta.get("genes_stamp") else None
    genes = []
    if models is not None:
        for g in models.overlapping(chrom, start, end):
            genes.append({"symbol": g.symbol, "gene_id": g.gene_id, "biotype": g.biotype, "start": g.start, "end": g.end, "strand": g.strand,
                          "transcript": g.transcript, "model_source": g.model_source, "exons": g.exons, "cds": [g.cds_start, g.cds_end] if g.cds_start else None})
    con = _connect(c)
    try:
        rows = con.execute(
            "SELECT id, pos, end, svtype, size, f_copies, impact, imprecise FROM svs "
            "WHERE chrom = ? AND pos <= ? AND end >= ? AND passed = 1 AND svtype != 'BND' AND f_copies > 0 ORDER BY pos LIMIT ?",
            (chrom, end, start, limit),
        ).fetchall()
    finally:
        con.close()
    return {
        "chrom": chrom, "start": start, "end": end, "genes": genes[:200], "n_genes": len(genes),
        "svs": [{"id": r["id"], "pos": r["pos"], "end": r["end"], "svtype": r["svtype"], "size": r["size"], "copies": r["f_copies"], "impact": r["impact"], "imprecise": bool(r["imprecise"])} for r in rows],
    }


def gene_info(symbol: str) -> Optional[Dict[str, Any]]:
    models = gene_models()
    g = models.by_symbol.get(symbol) if models is not None else None
    if g is None:
        return None
    return {
        "symbol": g.symbol, "gene_id": g.gene_id, "biotype": g.biotype, "chrom": g.chrom, "start": g.start, "end": g.end, "strand": g.strand,
        "transcript": g.transcript, "ensembl_transcript": g.ensembl_transcript, "model_source": g.model_source, "n_exons": len(g.exons),
        "pli": g.pli, "oe_lof": g.oe_lof, "inheritance": list(g.inheritance), "hpo_n": g.hpo_n, "hpo_terms": list(g.hpo_terms), "diseases": list(g.diseases),
    }


_overview_cache: Dict[str, Dict[str, Any]] = {}


def overview(c: Callset) -> Dict[str, Any]:
    meta = _require_ready(c)
    cached = _overview_cache.get(c.key)
    if cached is not None and cached["_built"] == meta["built"]:
        return cached
    qc = json.loads(meta["qc"])
    have_pop, have_genes, have_tracks = bool(meta.get("catalogs_stamp")), bool(meta.get("genes_stamp")), bool(meta.get("tracks_stamp"))
    trio = len(qc["roles"]["parents"]) == 2
    con = _connect(c)
    try:
        def count(where: str) -> int:
            return con.execute(f"SELECT COUNT(*) FROM svs WHERE {where}").fetchone()[0]

        steps = [("All calls in the VCF", "1=1"), ("Pass the caller's FILTER", "passed = 1"), ("Not a breakend", "svtype != 'BND'")]
        steps.append(("Carried by the proband", "f_copies > 0"))
        if have_pop:
            steps.append(("Rare: population AF ≤ 1% or not found", "(pop_af_max IS NULL OR pop_af_max <= 0.01)"))
        if have_genes:
            steps.append(("Hits a gene exon or a whole gene", "impact_rank >= 3"))
        if have_genes and have_tracks:
            steps.append(("Outside segmental duplications and repeats", "(segdup_frac IS NULL OR segdup_frac < 0.5) AND (repeat IS NULL OR repeat_frac < 0.5)"))
        if trio:
            steps.append(("De novo candidates", "inh = 'de_novo'"))
        clauses, funnel = [], []
        for label, clause in steps:
            clauses.append(clause)
            funnel.append({"label": label, "count": count(" AND ".join(clauses))})
        top = []
        if have_genes:
            top = [{"gene": r[0], "count": r[1]} for r in con.execute(
                "SELECT g.symbol, COUNT(*) FROM svs s JOIN sv_genes g ON g.sv_id = s.id WHERE s.passed = 1 AND s.svtype != 'BND' AND s.f_copies > 0 "
                "AND g.rank >= 3 AND g.coding = 1 AND (s.pop_af_max IS NULL OR s.pop_af_max <= 0.01) GROUP BY g.symbol ORDER BY COUNT(*) DESC, g.symbol LIMIT 12")]
    finally:
        con.close()
    out = {"_built": meta["built"], "sample": c.key, "qc": qc, "funnel": funnel, "top_genes": top, "resources": resources(meta), "trio": trio}
    _overview_cache[c.key] = out
    return out


def meta_for(c: Callset) -> Dict[str, Any]:
    st = status(c)
    base: Dict[str, Any] = {"sample": c.key, "label": c.info.get("label") or c.key, "reference": "GRCh38", "caller": c.info.get("caller"),
                            "vcf": os.path.basename(c.vcf), "status": st}
    if st["state"] == "ready":
        meta = _current_meta(c)
        if meta:
            qc = json.loads(meta["qc"])
            base.update({"built": meta["built"], "n_svs": int(meta["n_svs"]), "roles": qc["roles"], "source": qc.get("source"), "samples": qc["samples"]})
    return base
