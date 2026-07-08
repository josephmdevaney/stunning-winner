"""FastAPI routes for MitoLR, mounted on the shared app under /api/mito.

Serves the seeded mtDNA call set, per-variant detail (with an IGV-style read
pileup and circular-map inputs), and the sample QC dashboard. The variant data
is generated once at import and cached; a production build would swap the
``generate``/``qc`` modules for a mitochondrial variant caller's output joined
against a pathogenic-variant knowledge base.
"""

from __future__ import annotations

from functools import lru_cache

from fastapi import APIRouter, HTTPException

from . import generate, qc
from .reference import GENE_MAP, GENOME_LENGTH

router = APIRouter(prefix="/api/mito", tags=["mitolr"])

SAMPLE_ID = "MT-0417-P2"
INSTRUMENT_LINE = "PacBio HiFi · mtDNA amplicon"
REFERENCE = "rCRS (NC_012920.1)"
RUN_DATE = "2026-07-02"
# Prototype placeholder caller name -- not a real tool to integrate against.
CALLER_LABEL = "MitoCaller v1.4"

# Heteroplasmy interpretation bands shown on the detail screen's reference card.
REFERENCE_RANGES = [
    {"label": "Not detected", "range": "<5%"},
    {"label": "Low-level heteroplasmy", "range": "5–30%"},
    {"label": "Clinically significant", "range": ">30%"},
]


@lru_cache(maxsize=1)
def _variants() -> list[dict]:
    return generate.build_variants()


@lru_cache(maxsize=1)
def _variants_by_id() -> dict[str, dict]:
    return {v["id"]: v for v in _variants()}


def _meta() -> dict:
    return {
        "sample": SAMPLE_ID,
        "instrument_line": INSTRUMENT_LINE,
        "reference": REFERENCE,
        "run_date": RUN_DATE,
        "caller": CALLER_LABEL,
    }


@router.get("/meta")
def get_meta() -> dict:
    return _meta()


@router.get("/overview")
def get_overview() -> dict:
    variants = _variants()
    haplo = generate.haplogroup()
    pathogenic = [v for v in variants if v["status"] == "Pathogenic"]
    consequence_options = sorted({v["consequence"] for v in variants})

    return {
        **_meta(),
        "stats": {
            "variants_called": len(variants),
            "pathogenic_findings": len(pathogenic),
            "haplogroup": haplo["name"],
            "qc_status": "Pass",
        },
        "variants": variants,
        "haplogroup": haplo,
        "consequence_options": consequence_options,
    }


@router.get("/variant/{variant_id}")
def get_variant(variant_id: str) -> dict:
    variant = _variants_by_id().get(variant_id)
    if variant is None:
        raise HTTPException(status_code=404, detail=f"unknown variant id {variant_id!r}")

    return {
        **_meta(),
        "variant": variant,
        "alt_base": generate.alt_base(variant),
        "pileup_reads": generate.simulate_pileup(variant),
        "pileup_source": (
            "Estimated: distributes the variant's read-depth across mutant vs wild-type "
            "reads in proportion to the reported heteroplasmy. A true per-read view "
            "requires the spanning-reads BAM."
        ),
        "gene_map": GENE_MAP,
        "genome_length": GENOME_LENGTH,
        "reference_ranges": REFERENCE_RANGES,
    }


@router.get("/qc")
def get_qc() -> dict:
    return {**_meta(), **qc.build_qc()}
