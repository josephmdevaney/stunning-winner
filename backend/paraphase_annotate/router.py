"""FastAPI routes for the Paraphase Annotation suite, mounted on the shared app
under /api/paraphase.

Runs PacBio's own interpret_paraphase.py (see runner.py) against a sample's
real Paraphase JSON + per-gene VCFs and serves its two output CSVs. All data
shaping lives in parser.py/runner.py; this module is a thin dispatch layer.
"""

from __future__ import annotations

from typing import Optional

from fastapi import APIRouter, HTTPException

from . import parser, regions, runner, samples

router = APIRouter(prefix="/api/paraphase", tags=["paraphase"])


def _resolve(sample: Optional[str]) -> str:
    if sample is None or sample not in samples.SAMPLES:
        raise HTTPException(status_code=404, detail=f"unknown sample {sample!r}")
    return sample


@router.get("/samples")
def get_samples():
    return {"samples": samples.sample_keys()}


@router.get("/status")
def get_status(sample: Optional[str] = None):
    return runner.status(_resolve(sample))


@router.post("/run")
def post_run(sample: Optional[str] = None):
    return runner.start_run(_resolve(sample))


@router.get("/run-info")
def get_run_info(sample: Optional[str] = None):
    """How this sample's annotation run was set up: the vendored script release,
    the Python it ran on, whether the JSON was adapted, and which Paraphase
    produced the data. `run` is null until a run has been started."""
    s = _resolve(sample)
    return {"sample": s, "script": runner.vendor_info(), "run": runner.run_meta(s)}


@router.get("/overview")
def get_overview(sample: Optional[str] = None):
    s = _resolve(sample)
    st = runner.status(s)
    if st["state"] != "done":
        return {"sample": s, **st, "regions": [], "run": runner.run_meta(s), "script": runner.vendor_info()}
    return {
        "sample": s,
        "state": "done",
        "regions": parser.parse_summary(runner.summary_csv(s)),
        "run": runner.run_meta(s),
        "script": runner.vendor_info(),
    }


@router.get("/region")
def get_region(sample: Optional[str] = None, region: Optional[str] = None):
    s = _resolve(sample)
    st = runner.status(s)
    if st["state"] != "done":
        return {"sample": s, "region": region, **st, "haplotypes": []}
    return {
        "sample": s,
        "region": region,
        "state": "done",
        "haplotypes": parser.parse_detailed(runner.detailed_csv(s), region),
    }


# --- Paraphase's own JSON, instantly (no annotation run needed) ---------------

@router.get("/json")
def get_json_table(sample: Optional[str] = None):
    """Every region in the sample's Paraphase JSON with its copy numbers, depth,
    haplotype counts and flags. Reads the file directly: no network, no wait."""
    s = _resolve(sample)
    if not samples.json_path(s).exists():
        raise HTTPException(status_code=404, detail=f"Paraphase JSON not found for {s!r}")
    try:
        return regions.table(s)
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))


@router.get("/json/region")
def get_json_region(sample: Optional[str] = None, region: Optional[str] = None):
    s = _resolve(sample)
    if not samples.json_path(s).exists():
        raise HTTPException(status_code=404, detail=f"Paraphase JSON not found for {s!r}")
    try:
        detail = regions.detail(s, region or "")
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))
    if detail is None:
        raise HTTPException(status_code=404, detail=f"region {region!r} not in {s!r}'s Paraphase JSON")
    return detail
