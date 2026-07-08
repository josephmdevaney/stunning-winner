"""Sample registry + per-sample data builders for MitoLR.

Two samples, mirroring RepeatScope's multi-sample UX:

* ``HG002`` — real PacBio HiFi calls parsed from a MitoRSaw VCF (default).
* ``MT-0417-P2`` — the synthetic clinical-demo sample (seeded), retained so the
  pathogenic-finding UI states are still demonstrable.

The builders return the exact shapes the frontend consumes, so the router is a
thin dispatch layer over ``sample`` selection.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Optional

from . import generate, qc, vcf
from .reference import GENE_MAP, GENOME_LENGTH, HAPLOGROUP

BACKEND_ROOT = Path(__file__).resolve().parent

REFERENCE_RANGES = [
    {"label": "Not detected", "range": "<5%"},
    {"label": "Low-level heteroplasmy", "range": "5–30%"},
    {"label": "Clinically significant", "range": ">30%"},
]

# key -> descriptor. `source` selects the builder.
SAMPLES: dict[str, dict] = {
    "HG002": {
        "label": "HG002 (MitoRSaw)",
        "source": "vcf",
        "vcf": str(BACKEND_ROOT / "fixtures" / "HG002.GRCh38.mitorsaw.vcf.gz"),
        # Real chrM QC + per-read pileups precomputed from the haplotagged BAM
        # (scripts/real_data/build_mito_bam_qc.py). Optional: falls back to
        # VCF-derived QC if absent.
        "bam_qc": str(BACKEND_ROOT / "fixtures" / "HG002.chrM.bam_qc.json"),
        "caller": "mitorsaw 0.2.3-cba3800",
        "instrument_line": "PacBio HiFi · GRCh38 chrM",
        "reference": "GRCh38 chrM (rCRS · NC_012920.1)",
        "run_date": "GIAB HG002",
    },
    "MT-0417-P2": {
        "label": "MT-0417-P2 (demo)",
        "source": "seeded",
        "caller": "MitoCaller v1.4 (synthetic)",
        "instrument_line": "PacBio HiFi · mtDNA amplicon",
        "reference": "rCRS (NC_012920.1)",
        "run_date": "2026-07-02",
    },
}

DEFAULT_SAMPLE = "HG002"


def sample_keys() -> list[str]:
    return list(SAMPLES)


def resolve(sample: Optional[str]) -> str:
    if sample is None:
        return DEFAULT_SAMPLE
    if sample not in SAMPLES:
        raise KeyError(sample)
    return sample


def meta(sample: str) -> dict:
    s = SAMPLES[sample]
    return {
        "sample": sample,
        "instrument_line": s["instrument_line"],
        "reference": s["reference"],
        "run_date": s["run_date"],
        "caller": s["caller"],
    }


@lru_cache(maxsize=None)
def variants(sample: str) -> list[dict]:
    s = SAMPLES[sample]
    if s["source"] == "vcf":
        return vcf.build_variants(s["vcf"])
    return generate.build_variants()


@lru_cache(maxsize=None)
def _bam_qc(sample: str) -> Optional[dict]:
    """Precomputed aligned-BAM chrM stats for this sample, if available."""
    path = SAMPLES[sample].get("bam_qc")
    if path and Path(path).exists():
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    return None


def variant_by_id(sample: str, variant_id: str) -> Optional[dict]:
    return next((v for v in variants(sample) if v["id"] == variant_id), None)


def _haplogroup(sample: str) -> dict:
    if SAMPLES[sample]["source"] == "seeded":
        return dict(HAPLOGROUP)
    # Real sample: haplogroup requires a phylogenetic classifier (Haplogrep /
    # PhyloTree), which is out of scope for a variant-call parser. Report that
    # honestly and surface the actual SNVs vs rCRS as informative markers.
    snv_markers = [
        f"{v['change'].split('>')[0]}{v['pos']}{v['change'].split('>')[1]}"
        for v in variants(sample)
        if ">" in v["change"]
    ]
    return {
        "name": "Unclassified",
        "confidence": None,
        "region": "Requires phylogenetic classification (e.g. Haplogrep against PhyloTree)",
        "markers": snv_markers,
    }


def overview(sample: str) -> dict:
    vs = variants(sample)
    haplo = _haplogroup(sample)
    pathogenic = [v for v in vs if v["status"] == "Pathogenic"]
    return {
        **meta(sample),
        "stats": {
            "variants_called": len(vs),
            "pathogenic_findings": len(pathogenic),
            "haplogroup": haplo["name"],
            "qc_status": "Pass",
        },
        "variants": vs,
        "haplogroup": haplo,
        "consequence_options": sorted({v["consequence"] for v in vs}),
    }


def variant_detail(sample: str, variant_id: str) -> Optional[dict]:
    v = variant_by_id(sample, variant_id)
    if v is None:
        return None

    bam = _bam_qc(sample)
    is_snv = len(v["change"]) == 3 and ">" in v["change"]  # e.g. "A>G"
    bam_var = bam["variants"].get(v["id"]) if bam else None

    if bam_var and is_snv:
        # Real per-read alleles from the aligned BAM. Used only for SNVs: naive
        # pileup at an anchor base can't quantify homopolymer indels (MitoRSaw's
        # haplotype-aware VAF is authoritative for those).
        pileup = bam_var["pileup"]
        pileup_source = (
            f"Real per-read alleles at chrM:{v['pos']} from the aligned HiFi BAM "
            f"({bam_var['alt_reads']} alt / {bam_var['ref_reads']} ref reads) — independently "
            f"consistent with MitoRSaw's reported heteroplasmy."
        )
    elif SAMPLES[sample]["source"] == "vcf":
        pileup = vcf.simulate_pileup(v)
        pileup_source = (
            "Per-allele read support from MitoRSaw's AD field (ref vs alt read counts), "
            "subsampled to a legible number of rows. For homopolymer indels the caller's "
            "haplotype-aware VAF is used rather than a raw pileup, which cannot resolve them."
        )
    else:
        pileup = generate.simulate_pileup(v)
        pileup_source = (
            "Estimated: distributes the variant's read depth across mutant vs wild-type "
            "reads in proportion to the reported heteroplasmy."
        )
    return {
        **meta(sample),
        "variant": v,
        "alt_base": generate.alt_base(v),
        "pileup_reads": pileup,
        "pileup_source": pileup_source,
        "gene_map": GENE_MAP,
        "genome_length": GENOME_LENGTH,
        "reference_ranges": REFERENCE_RANGES,
    }


def sample_qc(sample: str) -> dict:
    bam = _bam_qc(sample)
    if bam is not None:
        return {**meta(sample), **qc.build_bam_qc(bam)}
    if SAMPLES[sample]["source"] == "vcf":
        return {**meta(sample), **qc.build_vcf_qc(variants(sample))}
    return {**meta(sample), **qc.build_qc()}
