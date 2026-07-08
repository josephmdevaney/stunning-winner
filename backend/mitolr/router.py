"""FastAPI routes for MitoLR, mounted on the shared app under /api/mito.

Sample-aware: the ``sample`` query param selects between the real HG002 MitoRSaw
calls (default) and the synthetic clinical demo. All data shaping lives in
``samples``; this module is a thin dispatch + validation layer.
"""

from __future__ import annotations

from fastapi import APIRouter, HTTPException

from . import samples

router = APIRouter(prefix="/api/mito", tags=["mitolr"])


def _resolve(sample: str | None) -> str:
    try:
        return samples.resolve(sample)
    except KeyError:
        raise HTTPException(status_code=404, detail=f"unknown sample {sample!r}")


@router.get("/samples")
def get_samples():
    return {
        "samples": [{"key": k, "label": samples.SAMPLES[k]["label"]} for k in samples.sample_keys()],
        "default": samples.DEFAULT_SAMPLE,
    }


@router.get("/meta")
def get_meta(sample: str | None = None):
    return samples.meta(_resolve(sample))


@router.get("/overview")
def get_overview(sample: str | None = None):
    return samples.overview(_resolve(sample))


@router.get("/variant/{variant_id}")
def get_variant(variant_id: str, sample: str | None = None):
    detail = samples.variant_detail(_resolve(sample), variant_id)
    if detail is None:
        raise HTTPException(status_code=404, detail=f"unknown variant id {variant_id!r}")
    return detail


@router.get("/qc")
def get_qc(sample: str | None = None):
    return samples.sample_qc(_resolve(sample))
