"""FastAPI routes for the structural-variant suite, mounted on the shared app under /api/structvar.

Call-set-aware via the ``sample`` query param (default: the first discovered call set). Data routes answer
503 with the build status while an index is being built; the frontend polls ``/status`` and mounts its
screens once the call set is ready. All data shaping lives in ``samples``; this module only parses and
validates.
"""

from __future__ import annotations

from typing import Optional, Tuple

from fastapi import APIRouter, HTTPException

from . import samples
from .normalize import SV_TYPES, canon_chrom

router = APIRouter(prefix="/api/structvar", tags=["structvar"])


def _callset(key: Optional[str]) -> samples.Callset:
    try:
        return samples.resolve(key)
    except KeyError:
        raise HTTPException(status_code=404, detail=f"unknown call set {key!r}")


def _split(value: Optional[str]) -> Tuple[str, ...]:
    return tuple(v.strip() for v in (value or "").replace("\n", ",").split(",") if v.strip())


def _not_ready(e: samples.NotReady):
    raise HTTPException(status_code=503, detail=e.status)


@router.get("/samples")
def get_samples():
    return {"samples": [{"key": c.key, "label": c.info.get("label") or c.key} for c in samples.discover()], "default": samples.default_sample()}


@router.get("/status")
def get_status(sample: Optional[str] = None):
    """Starts the index build the first time a call set is asked about."""
    return samples.ensure_index(_callset(sample))


@router.get("/meta")
def get_meta(sample: Optional[str] = None):
    return samples.meta_for(_callset(sample))


@router.get("/overview")
def get_overview(sample: Optional[str] = None):
    try:
        return samples.overview(_callset(sample))
    except samples.NotReady as e:
        _not_ready(e)


@router.get("/svs")
def get_svs(
    sample: Optional[str] = None,
    pass_only: bool = True,
    include_bnd: bool = False,
    types: Optional[str] = None,
    min_size: Optional[int] = None,
    max_size: Optional[int] = None,
    carried: bool = True,
    zygosity: str = "any",
    inheritance: str = "any",
    min_parent_dp: Optional[int] = None,
    min_gq: Optional[int] = None,
    min_alt: Optional[int] = None,
    impact: str = "any",
    genes: Optional[str] = None,
    chrom: Optional[str] = None,
    max_pop_af: Optional[float] = None,
    novel: bool = False,
    max_segdup: Optional[float] = None,
    no_odd: bool = False,
    no_repeat: bool = False,
    precise_only: bool = False,
    sort: str = "position",
    order: str = "asc",
    page: int = 1,
    page_size: int = 50,
    explain: bool = False,
):
    if zygosity not in ("any", "het", "hom"):
        raise HTTPException(status_code=422, detail="zygosity must be any, het or hom")
    if inheritance not in ("any", "de_novo", "inherited", "unknown"):
        raise HTTPException(status_code=422, detail="unknown inheritance filter")
    if impact not in samples.IMPACT_FILTER_RANK:
        raise HTTPException(status_code=422, detail="unknown impact filter")
    f = samples.Filters(
        pass_only=pass_only, include_bnd=include_bnd, types=tuple(t for t in _split(types) if t in SV_TYPES), min_size=min_size, max_size=max_size,
        carried=carried, zygosity=zygosity, inheritance=inheritance, min_parent_dp=min_parent_dp, min_gq=min_gq, min_alt=min_alt, impact=impact,
        genes=_split(genes), chrom=canon_chrom(chrom) if chrom else "", max_pop_af=max_pop_af, novel=novel, max_segdup=max_segdup,
        no_odd=no_odd, no_repeat=no_repeat, precise_only=precise_only,
    )
    try:
        return samples.query_svs(_callset(sample), f, sort, order, page, page_size, explain)
    except samples.NotReady as e:
        _not_ready(e)


@router.get("/svs/{sv_id}")
def get_sv(sv_id: int, sample: Optional[str] = None):
    try:
        row = samples.get_sv(_callset(sample), sv_id)
    except samples.NotReady as e:
        _not_ready(e)
    if row is None:
        raise HTTPException(status_code=404, detail=f"no SV {sv_id}")
    return row


@router.get("/region")
def get_region(chrom: str, start: int, end: int, sample: Optional[str] = None):
    if end - start > 20_000_000:
        raise HTTPException(status_code=422, detail="region too large (limit 20 Mb)")
    try:
        return samples.region(_callset(sample), chrom, start, end)
    except samples.NotReady as e:
        _not_ready(e)


@router.get("/genes/{symbol}")
def get_gene(symbol: str):
    info = samples.gene_info(symbol)
    if info is None:
        raise HTTPException(status_code=404, detail=f"no gene {symbol!r}")
    return info
