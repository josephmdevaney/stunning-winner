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
from .reference import GENOME_LENGTH, MT_FEATURES, features_for_pos

# Region lengths (bp) for gating the strand-balance metric: it is approximated
# from read midpoints, which is unreliable for loci too small to catch enough
# read centers (the tRNAs), so those are shown as "—".
_REGION_LEN = {f["name"]: f["stop"] - f["start"] + 1 for f in MT_FEATURES}
_REGION_LEN["D-loop"] = (GENOME_LENGTH - 16024 + 1) + 576
_STRAND_MIN_LEN = 200


def build_bam_qc(bam: dict) -> dict:
    """Rich Sample QC from an aligned-BAM stats JSON (see scripts/real_data/
    build_mito_bam_qc.py). Read-length distribution, per-base coverage, mapQ,
    read quality, and per-region depth are all real; only a few run-level
    metrics remain out of reach from an mtDNA-only BAM slice."""
    n_pass = sum(1 for r in bam["regions"] if r["flag"] == "Pass")
    n_reg = len(bam["regions"])
    cov_max = bam.get("coverage_max") or 1000
    mean_len_kb = (bam["mean_read_length_bp"] / 1000) if bam.get("mean_read_length_bp") else None

    regions = []
    for r in bam["regions"]:
        length = _REGION_LEN.get(r["region"], 999)
        # Suppress the midpoint-approximated strand balance for tiny tRNA loci.
        strand = r["strandBalance"] if length >= _STRAND_MIN_LEN else "—"
        regions.append({**r, "strandBalance": strand})

    return {
        "overall_status": "Pass" if n_pass == n_reg else "Review",
        "stats": [
            {
                "label": "Overall status",
                "value": "Pass" if n_pass == n_reg else "Review",
                "sub": f"{n_pass}/{n_reg} regions ≥ {bam['min_required']}x",
                "dotColor": "oklch(0.55 0.14 150)" if n_pass == n_reg else "oklch(0.62 0.15 75)",
            },
            {
                "label": "Mean coverage",
                "value": f"{bam['coverage_mean']:,}x",
                "sub": f"min {bam['coverage_min']:,}x · max {bam['coverage_max']:,}x",
                "dotColor": "oklch(0.6 0.14 220)",
            },
            {
                "label": "chrM reads mapped",
                "value": f"{bam['total_reads']:,}",
                "sub": f"mean length {mean_len_kb:.1f}kb" if mean_len_kb else "—",
                "dotColor": "oklch(0.6 0.14 220)",
            },
            {
                "label": "Read quality",
                "value": bam.get("mean_read_q") or "—",
                "sub": f"mapQ {round(bam['mean_mapq'])}" if bam.get("mean_mapq") is not None else "—",
                "dotColor": "oklch(0.6 0.14 220)",
            },
        ],
        "read_length_histogram": bam["read_length_histogram"],
        "coverage_track": bam["coverage_track"],
        "coverage_note": "Per-base read depth across chrM, from the aligned HiFi BAM",
        "genome_length": bam["length"],
        "regions": regions,
        "min_required": bam["min_required"],
        "instrument": [
            {"label": "Caller", "value": "mitorsaw 0.2.3-cba3800"},
            {"label": "Reference", "value": "GRCh38 chrM (rCRS)"},
            {"label": "Mean HiFi read length", "value": f"{bam['mean_read_length_bp']:,} bp" if bam.get("mean_read_length_bp") else "—"},
            {"label": "Mean HiFi read Q", "value": bam.get("mean_read_q") or "—"},
        ],
        "coverage_y_axis_max_default": int(math.ceil(cov_max * 1.15 / 250) * 250),
        "unavailable_metrics": [
            {"label": "mtDNA copy number (copies/cell)", "reason": "Requires nuclear-normalized depth (autosomal coverage); the chrM slice alone can't yield copy number."},
            {"label": "Per-region strand balance for tRNA loci", "reason": "Approximated from read midpoints; loci under 200 bp catch too few read centers to be reliable, so they are shown as —."},
            {"label": "NUMT-derived read fraction", "reason": "MitoRSaw filters NUMT reads during calling; a per-read NUMT flag is not carried in this BAM slice."},
        ],
    }

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


def build_vcf_qc(variants: list[dict]) -> dict:
    """Sample QC honestly derived from a variant-call VCF.

    A VCF carries per-call read depth (DP) but not read-length distributions,
    genome-wide coverage, strand balance, or instrument QC -- those need the
    aligned BAM (or MitoRSaw's --output-hap-stats JSON). So this reports the
    depth-at-called-sites it *can* compute and lists the rest as unavailable,
    rather than fabricating a coverage track."""
    MIN_REQUIRED = 500
    depths = [v["depth"] for v in variants if v.get("depth") is not None]
    mean_depth = round(sum(depths) / len(depths)) if depths else None
    low_sites = [v for v in variants if (v.get("depth") or 0) < MIN_REQUIRED]

    # Depth at each called site, placed by real genomic position (a sparse depth
    # profile, NOT interpolated genome-wide coverage).
    coverage = [
        {"pos": v["pos"], "depth": v["depth"]}
        for v in sorted(variants, key=lambda x: x["pos"])
        if v.get("depth") is not None
    ]
    cov_max = max((c["depth"] for c in coverage), default=1000)

    # Per-region adequacy: aggregate called-site depth by rCRS feature. A site
    # inside an overlapping reading frame (MT-ATP8/MT-ATP6, MT-ND4L/MT-ND4)
    # counts toward both genes' coverage adequacy, not just the first-listed one.
    by_region: dict[str, list[int]] = {}
    for v in variants:
        if v.get("depth") is None:
            continue
        feats = features_for_pos(v["pos"])
        names = ["D-loop" if f["name"] == "MT-CR" else f["name"] for f in feats] or ["Intergenic"]
        for name in names:
            by_region.setdefault(name, []).append(v["depth"])
    regions = []
    for name, ds in sorted(by_region.items(), key=lambda kv: -max(kv[1])):
        mean = round(sum(ds) / len(ds))
        regions.append(
            {
                "region": f"{name} ({len(ds)} site{'s' if len(ds) != 1 else ''})",
                "depth": mean,
                "strandBalance": "—",  # not derivable from the VCF
                "flag": "Pass" if mean >= MIN_REQUIRED else "Low depth",
            }
        )

    return {
        "overall_status": "Pass" if not low_sites else "Review",
        "stats": [
            {
                "label": "Overall status",
                "value": "Pass" if not low_sites else "Review",
                "sub": f"{len(variants) - len(low_sites)}/{len(variants)} called sites ≥ {MIN_REQUIRED}x",
                "dotColor": "oklch(0.55 0.14 150)" if not low_sites else "oklch(0.62 0.15 75)",
            },
            {"label": "Variants called", "value": str(len(variants)), "sub": "chrM, vs rCRS", "dotColor": "oklch(0.6 0.14 220)"},
            {
                "label": "Mean call depth",
                "value": f"{mean_depth}x" if mean_depth is not None else "—",
                "sub": f"across {len(depths)} called sites",
                "dotColor": "oklch(0.6 0.14 220)",
            },
            {
                "label": "Low-depth sites",
                "value": str(len(low_sites)),
                "sub": f"below {MIN_REQUIRED}x" if low_sites else "None",
                "dotColor": "oklch(0.62 0.15 75)" if low_sites else "oklch(0.55 0.14 150)",
            },
        ],
        # Read-length distribution is unavailable from a VCF -> empty; the UI
        # renders an inline "unavailable" note in that card.
        "read_length_histogram": [],
        "coverage_track": coverage,
        "coverage_note": "Read depth (DP) at called sites — a sparse profile, not interpolated genome-wide coverage",
        "genome_length": GENOME_LENGTH,
        "regions": regions,
        "min_required": MIN_REQUIRED,
        "instrument": [
            {"label": "Caller", "value": "mitorsaw 0.2.3-cba3800"},
            {"label": "Reference", "value": "GRCh38 chrM (rCRS)"},
            {"label": "Contig length", "value": f"{GENOME_LENGTH:,} bp"},
            {"label": "Mean call depth", "value": f"{mean_depth}x" if mean_depth is not None else "—"},
        ],
        "coverage_y_axis_max_default": int(math.ceil(cov_max * 1.15 / 250) * 250),
        "unavailable_metrics": [
            {"label": "HiFi read-length distribution", "reason": "Requires the aligned BAM; not present in a variant-call VCF."},
            {"label": "Genome-wide per-base coverage", "reason": "The track above is depth at called sites only; a continuous coverage profile needs the BAM (e.g. mosdepth) or MitoRSaw's hap-stats JSON."},
            {"label": "Per-strand balance", "reason": "Strand-specific read counts are not emitted in this VCF."},
            {"label": "mtDNA copy number & instrument Q", "reason": "Require the aligned BAM / run-level metrics, not the mtDNA variant calls."},
        ],
    }


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
        "coverage_note": "Linear depth track, position 1–16,569 (rCRS)",
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
