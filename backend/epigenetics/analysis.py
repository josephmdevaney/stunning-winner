"""Per-region methylation statistics and status calls for imprinted DMRs.

Everything here is a pure function of CpG lists, so it can be tested without any files.

Estimator. Methylation over a region is the *pooled* fraction -- sum of estimated methylated counts over
the sum of estimated counts at every CpG in the window -- which weights each CpG by its read support. (It is
also what the EpiSignaLR paper uses for DMR-level values.) pb-CpG-tools "model" mode reports model-based
estimates rather than raw read counts; those are what is pooled here.

Calls. There is no consensus threshold for an abnormal imprinted-DMR methylation level (the 2024 MLID joint
statement says so explicitly: labs should use ranges derived from their own controls). The bands below are
therefore *defaults*, conventional and not validated on a cohort, and every result carries the thresholds
that produced it. A control-derived range can replace them per region.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Optional, Sequence

from .cpg import CpG

STATUS_IMPRINTED = "imprinted"        # an intact-looking imprint
STATUS_LOSS = "loss"                  # loss of methylation: both alleles unmethylated
STATUS_GAIN = "gain"                  # gain of methylation: both alleles methylated
STATUS_BORDERLINE = "borderline"      # neither clearly intact nor clearly disturbed
STATUS_LOW_COVERAGE = "low_coverage"  # too few CpGs / too little depth to call
STATUS_NO_DATA = "no_data"

ABERRANT = (STATUS_LOSS, STATUS_GAIN)


@dataclass(frozen=True)
class Thresholds:
    min_cpgs: int = 5                  # CpGs needed in the window to call at all
    min_mean_coverage: float = 8.0     # combined-track mean depth needed to call
    min_hap_mean_coverage: float = 4.0  # per-haplotype mean depth needed to use a haplotype
    loss_below: float = 20.0           # percent; below this, the region is unmethylated
    gain_above: float = 80.0           # percent; above this, the region is methylated
    imprint_low: float = 30.0          # unphased: hemimethylation band that is consistent with an imprint
    imprint_high: float = 70.0
    allelic_min_diff: float = 40.0     # percentage points between haplotypes that counts as allele-specific

    def as_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class TrackStats:
    n_cpgs: int
    mean_coverage: float
    methylation: float  # pooled percent, 0-100

    def as_dict(self) -> dict:
        return {"n_cpgs": self.n_cpgs, "mean_coverage": round(self.mean_coverage, 1), "methylation": round(self.methylation, 1)}


def track_stats(cpgs: Sequence[CpG]) -> Optional[TrackStats]:
    """Pooled methylation over `cpgs`, or None when there is no estimated read support at all."""
    if not cpgs:
        return None
    total = sum(c.est_mod + c.est_unmod for c in cpgs)
    if total <= 0:
        return None
    return TrackStats(
        n_cpgs=len(cpgs),
        mean_coverage=sum(c.coverage for c in cpgs) / len(cpgs),
        methylation=100.0 * sum(c.est_mod for c in cpgs) / total,
    )


def _usable(stats: Optional[TrackStats], min_cpgs: int, min_cov: float) -> bool:
    return stats is not None and stats.n_cpgs >= min_cpgs and stats.mean_coverage >= min_cov


def classify(
    combined: Optional[TrackStats],
    hap1: Optional[TrackStats],
    hap2: Optional[TrackStats],
    th: Thresholds = Thresholds(),
) -> dict:
    """Status, allele pattern and notes for one region.

    With two usable haplotypes the call is allele-aware: a large methylation difference between them is the
    imprint itself (one allele methylated, the other not), whatever the pooled value. Without them it falls
    back to the pooled value alone, which is weaker evidence and is labelled "unphased"."""
    notes: list[str] = []
    if combined is None:
        return {"status": STATUS_NO_DATA, "phasing": "none", "pattern": None, "allelic_difference": None,
                "notes": ["No CpGs with read support in this window."]}
    if not _usable(combined, th.min_cpgs, th.min_mean_coverage):
        why = []
        if combined.n_cpgs < th.min_cpgs:
            why.append(f"{combined.n_cpgs} CpGs (need {th.min_cpgs})")
        if combined.mean_coverage < th.min_mean_coverage:
            why.append(f"mean depth {combined.mean_coverage:.1f}x (need {th.min_mean_coverage:g}x)")
        return {"status": STATUS_LOW_COVERAGE, "phasing": "none", "pattern": None, "allelic_difference": None,
                "notes": ["Not enough data to call: " + ", ".join(why) + "."]}

    phased = _usable(hap1, th.min_cpgs, th.min_hap_mean_coverage) and _usable(hap2, th.min_cpgs, th.min_hap_mean_coverage)
    m = combined.methylation
    if phased:
        a, b = hap1.methylation, hap2.methylation  # type: ignore[union-attr]
        diff = abs(a - b)
        if diff >= th.allelic_min_diff:
            status, pattern = STATUS_IMPRINTED, "allele-specific"
        elif max(a, b) < th.loss_below:
            status, pattern = STATUS_LOSS, "both alleles unmethylated"
        elif min(a, b) > th.gain_above:
            status, pattern = STATUS_GAIN, "both alleles methylated"
        else:
            status, pattern = STATUS_BORDERLINE, "alleles similar, intermediate"
        return {"status": status, "phasing": "phased", "pattern": pattern, "allelic_difference": round(diff, 1), "notes": notes}

    if (hap1 is not None) or (hap2 is not None):
        notes.append("Haplotype-resolved methylation was not usable for this region (too few phased CpGs or too little phased depth).")
    if m < th.loss_below:
        status = STATUS_LOSS
    elif m > th.gain_above:
        status = STATUS_GAIN
    elif th.imprint_low <= m <= th.imprint_high:
        status = STATUS_IMPRINTED
    else:
        status = STATUS_BORDERLINE
    return {"status": status, "phasing": "unphased", "pattern": "unphased", "allelic_difference": None, "notes": notes}


def region_result(
    region: dict,
    combined: Sequence[CpG],
    hap1: Sequence[CpG],
    hap2: Sequence[CpG],
    th: Thresholds = Thresholds(),
    window: str = "dmr",
) -> dict:
    """One region's stats and call, shaped for the API. `region` is an entry of imprinted_dmrs.json."""
    c, h1, h2 = track_stats(combined), track_stats(hap1), track_stats(hap2)
    call = classify(c, h1, h2, th)
    return {
        "id": region["id"], "name": region["name"], "chrom": region["chrom"],
        "category": region["category"], "primary_diagnostic": region["primary_diagnostic"],
        # The joint statement calls category-X DMRs (polymorphic or tissue-variable) not interpretable for
        # imprinting disorders: they are measured and shown, but never counted toward a disorder or MLID.
        "interpretable": region["category"] != "X",
        "methylated_allele": region["methylated_allele"], "disorders": region["disorders"],
        "disorder_notes": region["disorder_notes"],
        "window": window,
        "combined": c.as_dict() if c else None,
        "hap1": h1.as_dict() if h1 else None,
        "hap2": h2.as_dict() if h2 else None,
        **call,
    }
