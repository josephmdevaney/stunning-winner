"""Sample QC dashboard data for MitoLR.

Sequencing-run-level QC for a mtDNA amplicon / enrichment PacBio HiFi sample:
read-length distribution, per-position coverage across the 16,569 bp genome,
and per-region coverage adequacy. Values are seeded (see generate.seeded_rand)
so they are stable; in production they would come from the aligned BAM/CRAM and
a coverage tool. Metrics that this pipeline genuinely cannot produce are listed
under ``unavailable_metrics`` rather than fabricated -- mirroring RepeatScope.
"""

from __future__ import annotations

import math

from .generate import seeded_rand
from .reference import GENOME_LENGTH

# Per-region coverage adequacy (matches the design's curated QC table).
_QC_REGIONS = [
    {"region": "D-loop / control region", "depth": 5120, "strandBalance": "0.51", "flag": "Pass"},
    {"region": "MT-RNR1 (12S rRNA)", "depth": 4980, "strandBalance": "0.49", "flag": "Pass"},
    {"region": "MT-RNR2 (16S rRNA)", "depth": 4710, "strandBalance": "0.50", "flag": "Pass"},
    {"region": "MT-TL1 (tRNA-Leu)", "depth": 5340, "strandBalance": "0.52", "flag": "Pass"},
    {"region": "MT-ATP6/8", "depth": 4460, "strandBalance": "0.48", "flag": "Pass"},
    {"region": "MT-ND4", "depth": 3980, "strandBalance": "0.47", "flag": "Pass"},
    {"region": "MT-ND6", "depth": 2210, "strandBalance": "0.39", "flag": "Low strand balance"},
    {"region": "MT-CYB", "depth": 4890, "strandBalance": "0.51", "flag": "Pass"},
]

_INSTRUMENT = [
    {"label": "Instrument", "value": "Revio · SMRT Cell 25M"},
    {"label": "Mean HiFi read Q", "value": "Q32 (99.94% accuracy)"},
    {"label": "mtDNA copy number", "value": "312 copies/cell"},
    {"label": "Mean mtDNA coverage", "value": "4,820x"},
]


def _read_length_histogram(rng) -> list[dict]:
    """~18 bins, bell-shaped around the full-length ~16.5 kb amplicon."""
    bins = []
    for i in range(18):
        center = 2 + i * 1.1
        peak = math.exp(-((center - 16.5) ** 2) / 6) * 100
        bins.append(
            {"kb": round(center, 1), "count": max(2, math.floor(peak + rng() * 5 + 0.5))}
        )
    return bins


def _coverage_track(rng) -> list[dict]:
    """Linear depth vs. position, 60 sampled points across the genome."""
    n_pts = 60
    pts = []
    for i in range(n_pts):
        base = 4500 + math.sin(i / n_pts * math.pi * 2) * 800
        depth = max(500, math.floor(base + (rng() - 0.5) * 1200 + 0.5))
        pos = round(i / (n_pts - 1) * GENOME_LENGTH)
        pts.append({"pos": pos, "depth": depth})
    return pts


def build_qc() -> dict:
    rng = seeded_rand(42)
    hist = _read_length_histogram(rng)
    coverage = _coverage_track(rng)

    flagged = [r for r in _QC_REGIONS if r["flag"] != "Pass"]
    passing = len(_QC_REGIONS) - len(flagged)

    return {
        "overall_status": "Pass" if not flagged else "Review",
        "stats": [
            {
                "label": "Overall status",
                "value": "Pass" if not flagged else "Review",
                "sub": f"{passing}/{len(_QC_REGIONS)} regions sufficiently covered",
                "dotColor": "oklch(0.55 0.14 150)" if not flagged else "oklch(0.62 0.15 75)",
            },
            {"label": "Mean HiFi read length", "value": "15.9kb", "sub": "N50 16.4kb", "dotColor": "oklch(0.6 0.14 220)"},
            {"label": "Total reads mapped", "value": "184K", "sub": "mtDNA-mapped, MAPQ≥20", "dotColor": "oklch(0.6 0.14 220)"},
            {
                "label": "Flagged regions",
                "value": str(len(flagged)),
                "sub": (flagged[0]["region"].split(" (")[0] + " — low strand balance") if flagged else "None",
                "dotColor": "oklch(0.62 0.15 75)" if flagged else "oklch(0.55 0.14 150)",
            },
        ],
        "read_length_histogram": hist,
        "coverage_track": coverage,
        "genome_length": GENOME_LENGTH,
        "regions": _QC_REGIONS,
        "min_required": 1000,
        "instrument": _INSTRUMENT,
        "coverage_y_axis_max_default": 6000,
        "unavailable_metrics": [
            {
                "label": "Single-molecule CpG methylation track",
                "reason": "Per-read 5mCpG calls require the aligned kinetics/BAM; the variant summary alone does not carry them.",
            },
            {
                "label": "NUMT (nuclear-mitochondrial insertion) contamination estimate",
                "reason": "Distinguishing authentic mtDNA reads from NUMT-derived reads requires whole-genome alignment, not an mtDNA-only amplicon summary.",
            },
        ],
    }
