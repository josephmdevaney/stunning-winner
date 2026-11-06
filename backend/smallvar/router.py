"""FastAPI routes for the small-variant suite, mounted on the shared app under /api/smallvar.

Sample-aware via the ``sample`` query param (default: first discovered sample). Data routes
answer 503 with the build status while a sample's index is being built; the frontend polls
``/status`` and only mounts its screens once the sample is ready. All data shaping lives in
``samples``; this module is a thin parsing and validation layer.
"""

from __future__ import annotations

from typing import Optional, Tuple

from fastapi import APIRouter, HTTPException

from . import samples
from .normalize import canon_chrom

router = APIRouter(prefix="/api/smallvar", tags=["smallvar"])

_KINDS = {"snv", "mnv", "ins", "del", "complex"}


def _sample(key: Optional[str]) -> samples.Sample:
    try:
        return samples.resolve(key)
    except KeyError:
        raise HTTPException(status_code=404, detail=f"unknown sample {key!r}")


def _split(value: Optional[str]) -> Tuple[str, ...]:
    return tuple(v.strip() for v in (value or "").replace("\n", ",").split(",") if v.strip())


@router.get("/samples")
def get_samples():
    return {
        "samples": [{"key": s.key, "label": s.info.get("label") or s.key} for s in samples.discover()],
        "default": samples.default_sample(),
    }


@router.get("/status")
def get_status(sample: Optional[str] = None):
    """Starts the index build the first time a sample is asked about."""
    return samples.ensure_index(_sample(sample))


@router.get("/meta")
def get_meta(sample: Optional[str] = None):
    return samples.meta_for(_sample(sample))


@router.get("/overview")
def get_overview(sample: Optional[str] = None):
    try:
        return samples.overview(_sample(sample))
    except samples.NotReady as e:
        raise HTTPException(status_code=503, detail=e.status)


@router.get("/variants")
def get_variants(
    sample: Optional[str] = None,
    pass_only: bool = True,
    kinds: Optional[str] = None,
    zygosity: str = "any",
    chrom: Optional[str] = None,
    in_gene: bool = False,
    genes: Optional[str] = None,
    clinvar: str = "any",
    min_stars: int = 0,
    max_af: Optional[float] = None,
    af_missing: str = "include",
    max_hprc_af: Optional[float] = None,
    min_gq: Optional[int] = None,
    min_dp: Optional[int] = None,
    sort: str = "position",
    order: str = "asc",
    page: int = 1,
    page_size: int = 50,
    explain: bool = False,
):
    kind_list = tuple(k for k in _split(kinds) if k in _KINDS)
    if zygosity not in ("any", "het", "hom"):
        raise HTTPException(status_code=422, detail="zygosity must be any, het or hom")
    if clinvar not in ("any", "plp", "plp_conflicting", "has_match", "no_match"):
        raise HTTPException(status_code=422, detail="unknown clinvar filter")
    if af_missing not in ("include", "exclude"):
        raise HTTPException(status_code=422, detail="af_missing must be include or exclude")
    f = samples.Filters(
        pass_only=pass_only, kinds=kind_list, zygosity=zygosity, chrom=canon_chrom(chrom) if chrom else "", in_gene=in_gene,
        genes=_split(genes), clinvar=clinvar, min_stars=max(0, min(min_stars, 4)), max_af=max_af,
        af_missing=af_missing, max_hprc_af=max_hprc_af, min_gq=min_gq, min_dp=min_dp,
    )
    try:
        return samples.query_variants(_sample(sample), f, sort, order, page, page_size, explain)
    except samples.NotReady as e:
        raise HTTPException(status_code=503, detail=e.status)


@router.get("/variants/{variant_id}")
def get_variant(variant_id: int, sample: Optional[str] = None):
    try:
        row = samples.get_variant(_sample(sample), variant_id)
    except samples.NotReady as e:
        raise HTTPException(status_code=503, detail=e.status)
    if row is None:
        raise HTTPException(status_code=404, detail=f"no variant {variant_id}")
    return row


@router.get("/gene-presets")
def get_gene_presets():
    return {"presets": samples.gene_presets()}
