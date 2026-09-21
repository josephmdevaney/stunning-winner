"""FastAPI routes for the Paraphase Annotation suite, mounted on the shared app
under /api/paraphase.

Runs PacBio's own interpret_paraphase.py (see runner.py) against a sample's
real Paraphase JSON + per-gene VCFs and serves its two output CSVs. All data
shaping lives in parser.py/runner.py; this module is a thin dispatch layer.
"""

from __future__ import annotations

from typing import Optional

from fastapi import APIRouter, HTTPException

from . import parser, runner, samples

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


@router.get("/overview")
def get_overview(sample: Optional[str] = None):
    s = _resolve(sample)
    st = runner.status(s)
    if st["state"] != "done":
        return {"sample": s, **st, "regions": []}
    return {
        "sample": s,
        "state": "done",
        "regions": parser.parse_summary(runner.summary_csv(s)),
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
