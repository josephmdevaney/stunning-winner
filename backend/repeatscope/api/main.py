"""FastAPI backend serving annotated TRGT data to the RepeatScope frontend.

Supports multiple assay "modes" (Genome-wide, PureTarget V2, Carrier, ...), each with
its own directory of single-sample TRGT VCFs and its own locus panel scope, configured
via `config/modes.json` (or REPEATSCOPE_MODES_CONFIG to point elsewhere):

    {
      "genome_wide":  {"label": "Genome-wide",  "vcf_dir": "fixtures/hg002_real",         "panel_bed": null},
      "puretarget_v2":{"label": "PureTarget V2", "vcf_dir": "fixtures/gv_puretarget_panel", "panel_bed": "fixtures/panels/puretarget_v2.repeat_definitions.hg38.bed"}
    }

`vcf_dir` holds one single-sample .vcf/.vcf.gz per sample (sample name read from the
VCF's own header). `panel_bed` (optional) is a repeat-catalog-format BED (`ID=;MOTIFS=;
STRUC=`) whose IDs scope that mode's Overview/QC to just that assay's loci; omit it (or
set null) to show STRchive's full 80-locus catalog, e.g. for genome-wide sequencing.

Every screen mirrors what's actually derivable from a repeat-calling VCF. Metrics that
require aligned-read data (mapq, genome-wide coverage, GC bias, read-length
histograms, ...) are reported under `unavailable_metrics` with a reason, rather than
fabricated.
"""

from __future__ import annotations

import dataclasses
import json
import os
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware

from ..annotate import (
    AnnotatedLocus,
    annotate_strchive_panel,
    iter_annotated_catalog,
    resolve_panel_locus_ids,
    simulate_pileup_reads,
)
from ..trgt import sample_names

REPO_ROOT = Path(__file__).resolve().parents[3]
BACKEND_ROOT = REPO_ROOT / "backend"

MODES_CONFIG_PATH = os.environ.get("REPEATSCOPE_MODES_CONFIG", str(BACKEND_ROOT / "config" / "modes.json"))
STRCHIVE_JSON = os.environ.get(
    "REPEATSCOPE_STRCHIVE_JSON", str(REPO_ROOT / "project" / "uploads" / "STRchive-loci.json")
)
CATALOG_BED = os.environ.get(
    "REPEATSCOPE_CATALOG_BED", str(REPO_ROOT / "project" / "uploads" / "repeat_catalog.hg38.bed")
)

CALLER_LABEL = "TRGT v5.1.0"
VCF_EXTENSIONS = (".vcf", ".vcf.gz")

app = FastAPI(title="RepeatScope API")
app.add_middleware(
    CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"]
)

# MitoLR is a sibling product served by the same app under /api/mito. Its data is
# self-contained (no shared VCF/mode state with RepeatScope above).
from mitolr.router import router as mito_router  # noqa: E402

app.include_router(mito_router)

# Per-mode state: {mode_key: {label, sample_paths: {sample: vcf_path}, locus_ids: set|None}}
_modes: dict[str, dict] = {}
_default_mode: str = ""

# Caches keyed by (mode, sample) so switching between samples/modes in the UI is instant
# after the first load of each.
_panels: dict[tuple[str, str], list[AnnotatedLocus]] = {}
_panel_by_id: dict[tuple[str, str], dict[str, AnnotatedLocus]] = {}
_catalog_caches: dict[tuple[str, str], list[dict]] = {}


def _resolve_config_path(p: Optional[str]) -> Optional[str]:
    if not p:
        return None
    path = Path(p)
    if not path.is_absolute():
        path = BACKEND_ROOT / p
    return str(path)


def _discover_samples_in_dir(vcf_dir: str) -> dict[str, str]:
    paths = sorted(
        p for p in Path(vcf_dir).iterdir() if p.name.endswith(VCF_EXTENSIONS) and p.is_file()
    )
    found: dict[str, str] = {}
    for p in paths:
        names = sample_names(p)
        if len(names) != 1:
            raise RuntimeError(f"{p} has {len(names)} samples ({names}); expected exactly one")
        found[names[0]] = str(p)
    return found


@app.on_event("startup")
def _load_modes() -> None:
    global _modes, _default_mode
    with open(MODES_CONFIG_PATH, encoding="utf-8") as f:
        config = json.load(f)

    modes: dict[str, dict] = {}
    for key, cfg in config.items():
        vcf_dir = _resolve_config_path(cfg.get("vcf_dir"))
        sample_paths: dict[str, str] = {}
        if vcf_dir and Path(vcf_dir).is_dir():
            sample_paths = _discover_samples_in_dir(vcf_dir)

        panel_bed = _resolve_config_path(cfg.get("panel_bed"))
        locus_ids: Optional[set] = None
        if panel_bed and Path(panel_bed).exists():
            locus_ids = resolve_panel_locus_ids(panel_bed, STRCHIVE_JSON)

        non_repeat_targets: list[dict] = []
        nrt_path = _resolve_config_path(cfg.get("non_repeat_targets"))
        if nrt_path and Path(nrt_path).exists():
            with open(nrt_path, encoding="utf-8") as f:
                non_repeat_targets = json.load(f)

        modes[key] = {
            "label": cfg.get("label", key),
            "sample_paths": sample_paths,
            "locus_ids": locus_ids,
            "non_repeat_targets": non_repeat_targets,
        }
    _modes = modes
    _default_mode = next(
        (k for k, m in modes.items() if m["sample_paths"]), next(iter(modes), "")
    )


def _resolve_mode(mode: Optional[str]) -> str:
    if mode is None:
        return _default_mode
    if mode not in _modes:
        raise HTTPException(status_code=404, detail=f"unknown mode {mode!r}")
    return mode


def _resolve_sample(mode: str, sample: Optional[str]) -> str:
    paths = _modes[mode]["sample_paths"]
    if not paths:
        raise HTTPException(status_code=404, detail=f"no samples available for mode {mode!r}")
    if sample is None:
        return sorted(paths)[0]
    if sample not in paths:
        raise HTTPException(status_code=404, detail=f"unknown sample {sample!r} for mode {mode!r}")
    return sample


def _get_panel(mode: str, sample: str) -> list[AnnotatedLocus]:
    key = (mode, sample)
    if key not in _panels:
        m = _modes[mode]
        vcf_path = m["sample_paths"][sample]
        panel = annotate_strchive_panel(vcf_path, STRCHIVE_JSON, sample=sample, locus_ids=m["locus_ids"])
        _panels[key] = panel
        _panel_by_id[key] = {l.id: l for l in panel}
    return _panels[key]


def _chrom_sort_key(c: str):
    order = {f"chr{i}": i for i in range(1, 23)}
    order["chrX"] = 23
    order["chrY"] = 24
    order["chrM"] = 25
    return order.get(c, 99)


def _get_catalog(mode: str, sample: str) -> list[dict]:
    key = (mode, sample)
    if key not in _catalog_caches:
        vcf_path = _modes[mode]["sample_paths"][sample]
        regions = iter_annotated_catalog(CATALOG_BED, vcf_path, sample=sample)
        _catalog_caches[key] = [
            {
                "chrom": r.chrom, "start": r.start, "stop": r.stop, "length": r.length,
                "motifs": r.motifs, "struc": r.struc, "genotyped": r.genotyped,
                "a1": r.a1, "a2": r.a2, "depth": r.depth,
            }
            for r in regions
        ]
    return _catalog_caches[key]


def _locus_summary(l: AnnotatedLocus) -> dict:
    return {
        "id": l.id, "gene": l.gene, "disease": l.disease, "region": l.region, "motif": l.motif,
        "a1": l.a1, "a2": l.a2, "zygosity": l.zygosity, "depth": l.depth,
        "status": l.status, "status_color": l.status_color, "genotyped": l.genotyped,
        "note": l.note, "clickable": True,
    }


def _non_repeat_summary(target: dict) -> dict:
    """A carrier-panel target whose variant type (CNV, gene conversion, structural
    inversion, ...) a repeat-genotyping VCF can't express or call -- shown honestly
    as "Not applicable" rather than omitted or misrepresented as a repeat result."""
    return {
        "id": f"__non_repeat__{target['gene']}",
        "gene": target["gene"],
        "disease": target["disease"],
        "region": f"{target['chrom']}:{target['start']}-{target['stop']}",
        "motif": "—",
        "a1": None, "a2": None, "zygosity": None, "depth": None,
        "status": "Not applicable",
        "status_color": "oklch(0.6 0.02 260)",
        "genotyped": False,
        "note": f"Requires {target['tool']} -- {target['mechanism']} (not a tandem repeat).",
        "clickable": False,
    }


@app.get("/api/modes")
def get_modes():
    return {
        "modes": [
            {"key": k, "label": m["label"], "sample_count": len(m["sample_paths"])}
            for k, m in _modes.items()
        ],
        "default": _default_mode,
    }


@app.get("/api/samples")
def get_samples(mode: Optional[str] = None):
    m = _resolve_mode(mode)
    samples = sorted(_modes[m]["sample_paths"])
    return {"samples": samples, "default": samples[0] if samples else "", "mode": m}


@app.get("/api/meta")
def get_meta(mode: Optional[str] = None, sample: Optional[str] = None):
    m = _resolve_mode(mode)
    s = _resolve_sample(m, sample)
    return {"sample": s, "mode": m, "source_vcf": _modes[m]["sample_paths"][s], "caller": CALLER_LABEL}


@app.get("/api/overview")
def get_overview(mode: Optional[str] = None, sample: Optional[str] = None):
    m = _resolve_mode(mode)
    s = _resolve_sample(m, sample)
    panel = _get_panel(m, s)
    genotyped = [l for l in panel if l.genotyped]
    pathogenic = [l for l in panel if l.status == "Full expansion"]
    depths = [l.depth for l in genotyped if l.depth is not None]
    mean_depth = round(sum(depths) / len(depths)) if depths else None
    under_covered = [l for l in genotyped if (l.depth or 0) < 10]

    non_repeat_targets = _modes[m]["non_repeat_targets"]
    loci = [_locus_summary(l) for l in panel] + [_non_repeat_summary(t) for t in non_repeat_targets]

    return {
        "sample": s,
        "mode": m,
        "caller": CALLER_LABEL,
        "stats": {
            "loci_analyzed": len(panel) + len(non_repeat_targets),
            "loci_genotyped": len(genotyped),
            "pathogenic_findings": len(pathogenic),
            "mean_read_depth": mean_depth,
            "qc_status": "Pass" if not under_covered else f"{len(under_covered)} loci below 10x",
        },
        "loci": loci,
    }


@app.get("/api/locus/{locus_id}")
def get_locus(locus_id: str, mode: Optional[str] = None, sample: Optional[str] = None):
    m = _resolve_mode(mode)
    s = _resolve_sample(m, sample)
    _get_panel(m, s)
    l = _panel_by_id[(m, s)].get(locus_id)
    if l is None:
        raise HTTPException(status_code=404, detail=f"unknown locus id {locus_id!r}")

    d = dataclasses.asdict(l)
    d["region"] = l.region
    d["pileup_reads"] = simulate_pileup_reads(l)
    d["pileup_source"] = (
        "Estimated: distributes each allele's real TRGT read-support count (SD) across its "
        "real reported confidence interval (ALLR); TRGT's VCF output does not include raw "
        "per-read spanning measurements."
        if l.genotyped
        else "unavailable"
    )
    d["methylation_source"] = (
        "Real per-allele mean 5mCpG methylation from TRGT's AM field (not per-read; a "
        "per-read CpG track requires the spanning-reads BAM, which this pipeline does not have)."
        if l.genotyped
        else "unavailable"
    )
    d["caller"] = CALLER_LABEL
    return d


@app.get("/api/qc/sample")
def get_sample_qc(mode: Optional[str] = None, sample: Optional[str] = None):
    m = _resolve_mode(mode)
    s = _resolve_sample(m, sample)
    panel = _get_panel(m, s)
    rows = []
    for l in panel:
        if not l.genotyped:
            flag = "No call"
        elif (l.depth or 0) >= 10:
            flag = "Pass"
        else:
            flag = "Low depth"
        rows.append({"gene": l.gene, "id": l.id, "depth": l.depth, "min_required": 10, "flag": flag})

    genotyped_depths = [l.depth for l in panel if l.depth is not None]
    passing = sum(1 for r in rows if r["flag"] == "Pass")
    return {
        "sample": s,
        "mode": m,
        "overall_status": "Pass" if passing == len(rows) else "Review",
        "loci_passing": passing,
        "loci_total": len(rows),
        "mean_locus_depth": round(sum(genotyped_depths) / len(genotyped_depths)) if genotyped_depths else None,
        "rows": rows,
        "unavailable_metrics": [
            {"label": "Mean HiFi read length", "reason": "Not present in a TRGT VCF; requires the aligned BAM/CRAM."},
            {"label": "Total reads mapped", "reason": "Not present in a TRGT VCF; requires the aligned BAM/CRAM."},
            {"label": "Flanking region MAPQ", "reason": "TRGT's VCF does not report per-locus flanking mapq."},
        ],
    }


@app.get("/api/qc/run")
def get_run_qc(mode: Optional[str] = None, sample: Optional[str] = None):
    m = _resolve_mode(mode)
    s = _resolve_sample(m, sample)
    catalog = _get_catalog(m, s)
    genotyped = [r for r in catalog if r["genotyped"]]
    depths = [r["depth"] for r in genotyped if r["depth"] is not None]

    by_chrom: dict[str, list[dict]] = {}
    for r in catalog:
        by_chrom.setdefault(r["chrom"], []).append(r)
    chrom_depth = []
    for chrom in sorted(by_chrom, key=_chrom_sort_key):
        regions = by_chrom[chrom]
        ds = [r["depth"] for r in regions if r["depth"] is not None]
        chrom_depth.append(
            {
                "chrom": chrom,
                "mean_repeat_region_depth": round(sum(ds) / len(ds), 1) if ds else None,
                "n_regions": len(regions),
                "n_genotyped": len(ds),
            }
        )

    return {
        "sample": s,
        "mode": m,
        "repeat_regions_total": len(catalog),
        "repeat_regions_genotyped": len(genotyped),
        "repeat_region_mean_depth": round(sum(depths) / len(depths), 1) if depths else None,
        "mean_depth_by_chromosome": chrom_depth,
        "note": (
            "These figures are repeat-region depth derived from the TRGT VCF's SD field, "
            "not whole-genome coverage -- a repeat-calling VCF only reports depth at "
            "cataloged repeat regions, not genome-wide."
        ),
        "unavailable_metrics": [
            {
                "label": "Genome-wide per-base coverage distribution",
                "reason": "Requires aligned BAM/CRAM + a coverage tool (e.g. mosdepth); not in a repeat-calling VCF.",
            },
            {
                "label": "Coverage vs. GC content",
                "reason": "Requires aligned reads genome-wide; not in a repeat-calling VCF.",
            },
            {
                "label": "Whole-chromosome mean coverage",
                "reason": "The chart above is repeat-region depth only, not whole-chromosome coverage.",
            },
            {
                "label": "Alignment & read metrics (% mapped, N50, duplicate rate, ...)",
                "reason": "Requires the aligner's own stats or a BAM/CRAM QC tool; not in a TRGT VCF.",
            },
        ],
    }


@app.get("/api/catalog")
def get_catalog(
    mode: Optional[str] = None, sample: Optional[str] = None,
    chrom: str = "all", motif: str = "", min_len: int = 0,
    page: int = 0, page_size: int = 50,
):
    m = _resolve_mode(mode)
    s = _resolve_sample(m, sample)
    all_regions = _get_catalog(m, s)
    filtered = all_regions
    if chrom != "all":
        filtered = [r for r in filtered if r["chrom"] == chrom]
    if motif:
        needle = motif.upper()
        filtered = [r for r in filtered if needle in r["motifs"].upper()]
    if min_len:
        filtered = [r for r in filtered if r["length"] >= min_len]

    total = len(filtered)
    start = page * page_size
    page_rows = filtered[start : start + page_size]
    chroms = sorted({r["chrom"] for r in all_regions}, key=_chrom_sort_key)

    return {
        "total_all": len(all_regions),
        "total_matching": total,
        "unique_motifs": len({r["motifs"] for r in all_regions}),
        "chromosomes": ["all"] + chroms,
        "page": page,
        "page_size": page_size,
        "rows": page_rows,
    }
