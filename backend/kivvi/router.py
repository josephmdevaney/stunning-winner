"""FastAPI routes for the Kivvi suite, mounted on the shared app under /api/kivvi.

Sample-aware via the ``sample`` query param (default: first discovered sample).
All data shaping lives in ``samples``/``parse``; this module is a thin dispatch +
validation layer.
"""

from __future__ import annotations

from typing import Optional

from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse

from . import samples

router = APIRouter(prefix="/api/kivvi", tags=["kivvi"])


def _resolve(sample: Optional[str]) -> str:
    try:
        return samples.resolve(sample)
    except KeyError:
        raise HTTPException(status_code=404, detail=f"unknown sample {sample!r}")


def _target(target: str) -> str:
    if target not in samples.TARGETS:
        raise HTTPException(status_code=404, detail=f"unknown target {target!r}")
    return target


@router.get("/samples")
def get_samples():
    return {
        "samples": [{"key": k, "label": k} for k in samples.sample_keys()],
        "default": samples.default_sample(),
    }


# Annotations use typing.Optional: FastAPI evaluates them at runtime and the repo
# supports Python 3.9 (see mitolr/router.py).
@router.get("/meta")
def get_meta(sample: Optional[str] = None):
    return samples.meta(_resolve(sample))


@router.get("/overview")
def get_overview(sample: Optional[str] = None):
    return samples.overview(_resolve(sample))


@router.get("/qc")
def get_qc(sample: Optional[str] = None):
    return samples.sample_qc(_resolve(sample))


@router.get("/target/{target}")
def get_target(target: str, sample: Optional[str] = None):
    key = _resolve(sample)
    detail = samples.target_detail(key, _target(target))
    if detail is None:
        raise HTTPException(status_code=404, detail=f"no {target} Kivvi output for sample {key!r}")
    return detail


@router.get("/target/{target}/svg")
def get_svg(target: str, sample: Optional[str] = None):
    key = _resolve(sample)
    path = samples.svg_path(key, _target(target))
    if path is None:
        raise HTTPException(status_code=404, detail=f"no {target} SVG for sample {key!r}")
    return FileResponse(path, media_type="image/svg+xml")
