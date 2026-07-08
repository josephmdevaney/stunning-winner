"""Load the STRchive disease-locus catalog and classify repeat calls against it.

STRchive (https://github.com/dashnowlab/STRchive) publishes, per disease locus, the
GRCh38 coordinates plus benign/intermediate/pathogenic repeat-count thresholds. This
module reproduces the same locus fields used by the original design prototype's
`buildLoci()` (see `project/Repeat Expansion Analyzer.dc.html`), sourced directly from
the raw upstream STRchive-loci.json rather than the prototype's pre-condensed copy.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

from .motifs import is_cyclic_rotation


@dataclass(frozen=True)
class StrchiveLocus:
    id: str
    gene: str
    disease: str
    chrom: str
    start: int
    stop: int
    motif: str
    inheritance: str
    evidence: str
    benign_min: int | None
    benign_max: int | None
    intermediate_min: int | None
    intermediate_max: int | None
    pathogenic_min: int | None
    pathogenic_max: int | None
    # Most STRchive loci are single-motif (`motif` above covers them). A handful (RFC1,
    # the FAME genes, ...) instead have several candidate motifs that can occupy the same
    # repeat tract, where *which* motif expanded determines pathogenicity, not just length.
    pathogenic_motifs: tuple[str, ...] = field(default_factory=tuple)
    benign_motifs: tuple[str, ...] = field(default_factory=tuple)
    unknown_motifs: tuple[str, ...] = field(default_factory=tuple)

    @property
    def region(self) -> str:
        return f"{self.chrom}:{self.start}-{self.stop}"

    @property
    def is_motif_dependent(self) -> bool:
        """True only when STRchive gives us an explicit benign-motif contrast list.

        Nearly every locus has a non-empty `pathogenic_motifs` (even simple single-motif
        ones like HTT, where it's just `("CAG",)`), so that alone isn't a useful signal.
        What actually matters is whether we have something to positively rule an allele
        *in* as benign by motif identity -- without that, flagging an unrecognized motif
        as "uncertain significance" would be overreach for loci we simply lack full
        motif-classification data on, rather than a real ambiguity.
        """
        return bool(self.benign_motifs)


def _first_motif(raw: dict) -> str:
    for key in (
        "pathogenic_motif_reference_orientation",
        "reference_motif_reference_orientation",
        "benign_motif_reference_orientation",
    ):
        values = raw.get(key) or []
        if values:
            return values[0]
    return ""


def load_strchive_loci(json_path: str | Path) -> list[StrchiveLocus]:
    """Parse STRchive's raw per-locus JSON export into StrchiveLocus records.

    `json_path` should point at STRchive's `data/*loci.json` export (a list of
    per-locus dicts with fields like `start_hg38`, `pathogenic_motif_reference_orientation`,
    etc.) -- see https://github.com/dashnowlab/STRchive/tree/main/data.
    """
    with open(json_path, encoding="utf-8") as f:
        raw_loci = json.load(f)

    loci = []
    for raw in raw_loci:
        loci.append(
            StrchiveLocus(
                id=raw["id"],
                gene=raw["gene"],
                disease=raw["disease"],
                chrom=raw["chrom"],
                start=raw["start_hg38"],
                stop=raw["stop_hg38"],
                motif=_first_motif(raw),
                inheritance="/".join(raw.get("inheritance") or []),
                evidence="/".join(raw.get("evidence") or []),
                benign_min=raw.get("benign_min"),
                benign_max=raw.get("benign_max"),
                intermediate_min=raw.get("intermediate_min"),
                intermediate_max=raw.get("intermediate_max"),
                pathogenic_min=raw.get("pathogenic_min"),
                pathogenic_max=raw.get("pathogenic_max"),
                pathogenic_motifs=tuple(raw.get("pathogenic_motif_reference_orientation") or ()),
                benign_motifs=tuple(raw.get("benign_motif_reference_orientation") or ()),
                unknown_motifs=tuple(raw.get("unknown_motif_reference_orientation") or ()),
            )
        )
    return loci


# Status labels and their display colors, matching the original design's badge palette.
STATUS_NORMAL = "Normal"
STATUS_INTERMEDIATE = "Intermediate"
STATUS_PATHOGENIC = "Full expansion"
STATUS_UNCERTAIN = "Uncertain significance"

STATUS_COLORS = {
    STATUS_PATHOGENIC: "oklch(0.55 0.20 25)",
    STATUS_UNCERTAIN: "oklch(0.5 0.15 300)",
    STATUS_INTERMEDIATE: "oklch(0.62 0.15 75)",
    STATUS_NORMAL: "oklch(0.55 0.14 150)",
}

# Severity ranking used to pick the worse of two alleles' classifications.
_SEVERITY = {STATUS_NORMAL: 0, STATUS_INTERMEDIATE: 1, STATUS_UNCERTAIN: 2, STATUS_PATHOGENIC: 3}


@dataclass(frozen=True)
class Classification:
    status: str
    status_color: str
    note: str


def _is_contraction_locus(locus: StrchiveLocus) -> bool:
    """True if this locus's pathogenic range sits *below* its benign range.

    Almost all STRchive loci are expansion diseases (more repeats = worse). A couple
    (e.g. CHNG3_MIR7-2, HMNR7_VWA1) are contraction diseases instead -- fewer repeats
    than the benign reference is what's pathogenic. STRchive's `locus_tags` includes a
    "contraction" tag, but it's overloaded (e.g. SBMA_AR carries it too, for an entirely
    normal *expansion* disease where contractions of the expanded allele have merely
    been observed as a side phenomenon) -- so direction is inferred from the numbers
    themselves: if the pathogenic range's floor is below the benign range's ceiling,
    pathogenic values must be the *lower* ones.

    HMNR7_VWA1 is a special case handled separately by `_is_deviation_locus` --
    see there.
    """
    return (
        locus.pathogenic_min is not None
        and locus.benign_max is not None
        and locus.pathogenic_min < locus.benign_max
    )


def _is_deviation_locus(locus: StrchiveLocus) -> bool:
    """True for the one STRchive locus (HMNR7_VWA1) whose pathogenic range *contains*
    its single benign value, rather than sitting cleanly above or below it: benign is
    exactly 2 repeats, pathogenic is 1-3 -- i.e. GeneReviews' own description is "any
    deviation from 2 motifs is thought to be pathogenic", not a length threshold in
    either direction. Treating this with `_classify_length_only`'s min/max comparison
    is actively wrong: contraction-direction logic (call pathogenic when the smaller
    allele is <= pathogenic_max=3) flags the completely normal 2/2 genotype as
    pathogenic, since 2 <= 3.
    """
    return (
        locus.benign_min is not None
        and locus.benign_min == locus.benign_max
        and locus.pathogenic_min is not None
        and locus.pathogenic_max is not None
        and locus.pathogenic_min <= locus.benign_min <= locus.pathogenic_max
    )


def _classify_deviation_locus(locus: StrchiveLocus, a1: int | None, a2: int | None) -> Classification:
    benign = locus.benign_min  # benign_min == benign_max, enforced by _is_deviation_locus
    values = [v for v in (a1, a2) if v is not None]
    deviant = [v for v in values if v != benign]
    if deviant:
        worst = max(deviant, key=lambda v: abs(v - benign))
        return Classification(
            status=STATUS_PATHOGENIC,
            status_color=STATUS_COLORS[STATUS_PATHOGENIC],
            note=f"{worst} repeats deviates from the benign value ({benign}); any deviation is pathogenic",
        )
    return Classification(
        status=STATUS_NORMAL,
        status_color=STATUS_COLORS[STATUS_NORMAL],
        note=f"Both alleles at the benign value ({benign} repeats)",
    )


def _classify_length_only(locus: StrchiveLocus, a1: int | None, a2: int | None) -> Classification:
    """Classify by repeat count alone: the larger of the two alleles compared against
    pathogenic, then intermediate, then benign thresholds (smaller allele instead, for
    the rare contraction-direction loci)."""
    values = [v for v in (a1, a2) if v is not None]
    contraction = _is_contraction_locus(locus)
    called = min(values) if contraction else max(values)
    cmp_patho = (lambda v, t: v <= t) if contraction else (lambda v, t: v >= t)
    cmp_inter = (lambda v, t: v <= t) if contraction else (lambda v, t: v >= t)
    patho_threshold = locus.pathogenic_max if contraction else locus.pathogenic_min
    inter_threshold = locus.intermediate_max if contraction else locus.intermediate_min

    if patho_threshold is not None and cmp_patho(called, patho_threshold):
        direction = "≤" if contraction else "≥"
        return Classification(
            status=STATUS_PATHOGENIC,
            status_color=STATUS_COLORS[STATUS_PATHOGENIC],
            note=f"{'Below' if contraction else 'Above'} pathogenic threshold ({direction}{patho_threshold} repeats)",
        )
    if inter_threshold is not None and cmp_inter(called, inter_threshold):
        lo, hi = locus.intermediate_min, locus.intermediate_max
        lo = lo if lo is not None else "?"
        hi = hi if hi is not None else "?"
        return Classification(
            status=STATUS_INTERMEDIATE,
            status_color=STATUS_COLORS[STATUS_INTERMEDIATE],
            note=f"Intermediate allele ({lo}–{hi}); consider follow-up",
        )
    benign_max = locus.benign_max if locus.benign_max is not None else "?"
    return Classification(
        status=STATUS_NORMAL,
        status_color=STATUS_COLORS[STATUS_NORMAL],
        note=f"Within benign range (≤{benign_max})",
    )


def _classify_allele_by_motif(
    locus: StrchiveLocus, count: int, motif: str | None
) -> Classification:
    """Classify one allele of a motif-dependent locus (RFC1, FAME genes, ...).

    The motif actually present matters more than raw length: a benign motif expanded
    to any length isn't disease-causing, while a pathogenic-consistent motif's length
    still needs to clear the usual thresholds. A motif matching neither the known
    benign nor pathogenic list is uncertain significance -- a genuinely unclassified
    finding, not an assumed-safe default.
    """
    if motif is None:
        # No usable per-segment motif data (e.g. MC missing) -- fall back to length alone.
        return _classify_length_only(locus, count, None)

    if any(is_cyclic_rotation(motif, m) for m in locus.benign_motifs):
        return Classification(
            status=STATUS_NORMAL, status_color=STATUS_COLORS[STATUS_NORMAL],
            note=f"{motif} is a benign motif at this locus regardless of length ({count} repeats)",
        )

    if any(is_cyclic_rotation(motif, m) for m in locus.pathogenic_motifs):
        if locus.pathogenic_min is not None and count >= locus.pathogenic_min:
            return Classification(
                status=STATUS_PATHOGENIC, status_color=STATUS_COLORS[STATUS_PATHOGENIC],
                note=f"{motif} (pathogenic motif) expanded to {count} repeats, "
                f"above pathogenic threshold (≥{locus.pathogenic_min})",
            )
        if locus.intermediate_min is not None and count >= locus.intermediate_min:
            return Classification(
                status=STATUS_INTERMEDIATE, status_color=STATUS_COLORS[STATUS_INTERMEDIATE],
                note=f"{motif} (pathogenic motif) at {count} repeats, in the intermediate range "
                f"({locus.intermediate_min}–{locus.intermediate_max})",
            )
        return Classification(
            status=STATUS_NORMAL, status_color=STATUS_COLORS[STATUS_NORMAL],
            note=f"{motif} (pathogenic-consistent motif) within benign range ({count} repeats)",
        )

    return Classification(
        status=STATUS_UNCERTAIN, status_color=STATUS_COLORS[STATUS_UNCERTAIN],
        note=f"{motif} ({count} repeats) is not a characterized benign or pathogenic motif "
        "at this locus -- clinical significance unknown",
    )


def _classify_motif_aware(
    locus: StrchiveLocus,
    a1: int | None, a2: int | None,
    a1_motif: str | None, a2_motif: str | None,
) -> Classification:
    candidates = []
    if a1 is not None:
        candidates.append(_classify_allele_by_motif(locus, a1, a1_motif))
    if a2 is not None:
        candidates.append(_classify_allele_by_motif(locus, a2, a2_motif))
    return max(candidates, key=lambda c: _SEVERITY[c.status])


def classify(
    locus: StrchiveLocus,
    a1: int | None, a2: int | None,
    a1_motif: str | None = None, a2_motif: str | None = None,
) -> Classification:
    """Classify a genotype call against a locus's benign/intermediate/pathogenic thresholds.

    For most STRchive loci this is a straightforward length comparison (see
    `_classify_length_only`). For the minority of loci where pathogenicity depends on
    *which* motif expanded (RFC1, the FAME genes, ...) -- `locus.is_motif_dependent` --
    each allele's actual motif (`a1_motif`/`a2_motif`, from `TrgtCall.dominant_motifs()`)
    is checked against STRchive's own benign/pathogenic motif lists, with an allele
    matching neither reported as uncertain significance rather than silently guessed.
    HMNR7_VWA1 is a one-off third case (`_is_deviation_locus`): its pathogenic range
    isn't a length threshold at all, it's "any deviation from exactly 2 repeats."
    """
    values = [v for v in (a1, a2) if v is not None]
    if not values:
        return Classification(
            status=STATUS_NORMAL, status_color=STATUS_COLORS[STATUS_NORMAL],
            note="No genotype call at this locus",
        )

    if _is_deviation_locus(locus):
        return _classify_deviation_locus(locus, a1, a2)
    if locus.is_motif_dependent:
        return _classify_motif_aware(locus, a1, a2, a1_motif, a2_motif)
    return _classify_length_only(locus, a1, a2)
