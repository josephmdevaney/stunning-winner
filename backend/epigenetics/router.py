"""FastAPI routes for the Epigenetics suite, mounted on the shared app under /api/epigenetics.

Thin dispatch over service.py. Thresholds default to the documented research defaults and can be overridden
per request (the response always echoes the thresholds it used)."""

from __future__ import annotations

from typing import Optional

from fastapi import APIRouter, HTTPException

from . import samples, service
from .analysis import Thresholds
from .cpg import TabixError
from .episignature import service as episignature_service
from .episignature import vendor as episignature_vendor

router = APIRouter(prefix="/api/epigenetics", tags=["epigenetics"])


def _sample(key: Optional[str]) -> samples.Sample:
    registry = samples.load()
    if not registry:
        raise HTTPException(
            status_code=404,
            detail=f"no epigenetics samples are registered -- see {samples.config_path()} "
                   "(epigenetics.local.json.example shows the format)",
        )
    if key is None:
        return next(iter(registry.values()))
    if key not in registry:
        raise HTTPException(status_code=404, detail=f"unknown sample {key!r}")
    return registry[key]


def _window(window: str) -> str:
    if window not in ("core", "dmr"):
        raise HTTPException(status_code=422, detail="window must be 'core' or 'dmr'")
    return window


@router.get("/samples")
def get_samples():
    registry = samples.load()
    return {
        "samples": [{"key": s.key, "label": s.label, "note": s.note, "tracks": sorted(s.paths)} for s in registry.values()],
        "default": next(iter(registry), ""),
    }


@router.get("/imprinting")
def get_imprinting(
    sample: Optional[str] = None,
    window: str = "core",
    loss_below: Optional[float] = None,
    gain_above: Optional[float] = None,
    allelic_min_diff: Optional[float] = None,
):
    s = _sample(sample)
    base = Thresholds()
    th = Thresholds(
        loss_below=base.loss_below if loss_below is None else loss_below,
        gain_above=base.gain_above if gain_above is None else gain_above,
        allelic_min_diff=base.allelic_min_diff if allelic_min_diff is None else allelic_min_diff,
    )
    if not (0 <= th.loss_below < th.gain_above <= 100):
        raise HTTPException(status_code=422, detail="need 0 <= loss_below < gain_above <= 100")
    try:
        return service.overview(s, _window(window), th)
    except TabixError as e:
        raise HTTPException(status_code=503, detail=str(e))


@router.get("/imprinting/region")
def get_imprinting_region(sample: Optional[str] = None, id: Optional[str] = None, flank: int = 3000, window: str = "core"):
    s = _sample(sample)
    if not 0 <= flank <= 50000:
        raise HTTPException(status_code=422, detail="flank must be between 0 and 50000")
    try:
        out = service.region_track(s, id or "", flank, _window(window))
    except TabixError as e:
        raise HTTPException(status_code=503, detail=str(e))
    if out is None:
        raise HTTPException(status_code=404, detail=f"unknown imprinted DMR {id!r}")
    return out


@router.get("/episignatures")
def get_episignatures(sample: Optional[str] = None):
    """EpiSignaLR episignature scores for a genome. The models are fetched at setup, so until they are there this
    answers `installed: false` with the command to run instead of an error."""
    s = _sample(sample)
    if not episignature_vendor.installed():
        return {
            "installed": False, "setup_command": episignature_vendor.SETUP_COMMAND, "repo": episignature_vendor.REPO,
            "available": [{"id": k, "name": v[1], "gene": v[2]} for k, v in episignature_vendor.CATALOG.items()],
            "sample": s.key, "label": s.label, "models": [],
        }
    try:
        return {"installed": True, **episignature_service.overview(s)}
    except TabixError as e:
        raise HTTPException(status_code=503, detail=str(e))
    except ValueError as e:
        raise HTTPException(status_code=503, detail=f"the fetched EpiSignaLR files are inconsistent: {e}")
