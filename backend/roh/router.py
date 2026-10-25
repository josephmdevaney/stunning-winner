"""FastAPI routes for the ROH tool, mounted on the shared app under /api/roh.

Keyed by ``sample`` only -- ROH doesn't depend on a repeat panel or assay mode -- and only
samples in ``analysis.ROH_SAMPLES`` have a real ROH run.
"""

from __future__ import annotations

from typing import Optional

from fastapi import APIRouter, HTTPException

from . import analysis as roh

router = APIRouter(prefix="/api/roh", tags=["roh"])


def _resolve_roh_sample(sample: Optional[str]) -> str:
    if sample is None:
        keys = roh.sample_keys()
        if not keys:
            raise HTTPException(status_code=404, detail="no ROH samples available")
        return keys[0]
    if sample not in roh.ROH_SAMPLES:
        raise HTTPException(status_code=404, detail=f"no ROH run for sample {sample!r}")
    return sample


@router.get("/samples")
def get_roh_samples():
    keys = roh.sample_keys()
    return {"samples": keys, "default": keys[0] if keys else ""}


# Annotations use typing.Optional: FastAPI evaluates them at runtime and the repo supports
# Python 3.9 (see mitolr/router.py).
@router.get("/overview")
def get_roh_overview(sample: Optional[str] = None):
    s = _resolve_roh_sample(sample)
    return roh.overview_for_sample(s)


@router.get("/regions")
def get_roh_regions(
    sample: Optional[str] = None,
    chrom: str = "all",
    min_len_mb: float = 0,
    page: int = 0,
    page_size: int = 100,
):
    s = _resolve_roh_sample(sample)
    return roh.query_regions(s, chrom=chrom, min_len_mb=min_len_mb, page=page, page_size=page_size)
