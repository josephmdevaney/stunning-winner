"""Tie the STRchive catalog and a TRGT VCF together into annotated, UI-ready records.

This is the core of "take a repeat VCF and do annotation": for every STRchive disease
locus, find the sample's TRGT genotype call at that position and classify it against
STRchive's benign/intermediate/pathogenic thresholds. Everything here is either a
value read directly from the VCF, or a value clearly derived from one (and labeled as
such) -- nothing is randomly simulated the way the original design prototype was.
"""

from __future__ import annotations

import random
import zlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from . import strchive
from .catalog import CatalogRegion, iter_catalog_regions
from .strchive import StrchiveLocus, load_strchive_loci
from .trgt import TrgtCall, TrgtIndex


def _seed_for(key: str) -> int:
    return zlib.crc32(key.encode("utf-8"))


def resolve_panel_locus_ids(panel_bed_path: str | Path, strchive_json_path: str | Path) -> set[str]:
    """Map a panel BED's own region IDs to STRchive locus ids.

    Assay vendors' own catalog IDs occasionally drift from STRchive's naming (e.g. a
    panel's "OPDM_ABCD3" vs STRchive's "OPDM5_ABCD3" for the same locus) even when the
    coordinates are identical -- so this matches by exact ID first, falling back to
    chrom + coordinate overlap rather than silently dropping an unmatched region from
    the panel's scope.
    """
    loci = load_strchive_loci(strchive_json_path)
    by_id = {l.id: l for l in loci}
    resolved: set[str] = set()
    for region in iter_catalog_regions(panel_bed_path):
        if region.id in by_id:
            resolved.add(region.id)
            continue
        for locus in loci:
            if locus.chrom == region.chrom and abs(locus.start - region.start) < 200 and abs(locus.stop - region.stop) < 200:
                resolved.add(locus.id)
                break
    return resolved


@dataclass
class AnnotatedLocus:
    id: str
    gene: str
    disease: str
    chrom: str
    start: int
    stop: int
    motif: str
    inheritance: str
    evidence: str
    benign_min: Optional[int]
    benign_max: Optional[int]
    intermediate_min: Optional[int]
    intermediate_max: Optional[int]
    pathogenic_min: Optional[int]
    pathogenic_max: Optional[int]

    genotyped: bool = False
    trid: Optional[str] = None
    phased: bool = False
    a1: Optional[int] = None
    a2: Optional[int] = None
    a1_bp: Optional[int] = None
    a2_bp: Optional[int] = None
    a1_range_bp: Optional[tuple[int, int]] = None
    a2_range_bp: Optional[tuple[int, int]] = None
    a1_depth: Optional[int] = None
    a2_depth: Optional[int] = None
    depth: Optional[int] = None
    a1_purity: Optional[float] = None
    a2_purity: Optional[float] = None
    a1_methylation: Optional[float] = None
    a2_methylation: Optional[float] = None
    a1_motif: Optional[str] = None
    a2_motif: Optional[str] = None
    zygosity: Optional[str] = None
    status: str = "No call"
    status_color: str = "oklch(0.6 0.02 260)"
    note: str = "No TRGT genotype call found at this locus"

    @property
    def region(self) -> str:
        return f"{self.chrom}:{self.start}-{self.stop}"

    def to_dict(self) -> dict:
        d = dict(self.__dict__)
        d["region"] = self.region
        return d


def _annotate_locus(locus: StrchiveLocus, call: Optional[TrgtCall]) -> AnnotatedLocus:
    out = AnnotatedLocus(
        id=locus.id,
        gene=locus.gene,
        disease=locus.disease,
        chrom=locus.chrom,
        start=locus.start,
        stop=locus.stop,
        motif=locus.motif,
        inheritance=locus.inheritance,
        evidence=locus.evidence,
        benign_min=locus.benign_min,
        benign_max=locus.benign_max,
        intermediate_min=locus.intermediate_min,
        intermediate_max=locus.intermediate_max,
        pathogenic_min=locus.pathogenic_min,
        pathogenic_max=locus.pathogenic_max,
    )

    if call is None or not call.genotyped:
        return out

    # Motif-dependent loci (RFC1, ...) have several candidate motifs at the same tract;
    # targeting `locus.motif` specifically would return that one motif's count -- often
    # zero, if the sample's actual expansion is a *different* candidate motif entirely.
    # Summing all segments instead gives the true count regardless of which motif is
    # present (only one segment is ever nonzero in practice for these single-tract loci).
    a1, a2 = call.repeat_units(None if locus.is_motif_dependent else locus.motif)
    a1_bp, a2_bp = call.allele_lengths_bp
    a1_depth, a2_depth = call.allele_depths
    a1_range, a2_range = call.allele_length_ranges
    a1_purity, a2_purity = call.allele_purity_frac
    a1_meth, a2_meth = call.allele_methylation_pct
    a1_motif, a2_motif = call.dominant_motifs()

    out.genotyped = True
    out.trid = call.trid
    out.phased = call.phased
    out.a1, out.a2 = a1, a2
    out.a1_bp, out.a2_bp = a1_bp, a2_bp
    out.a1_range_bp, out.a2_range_bp = a1_range, a2_range
    out.a1_depth, out.a2_depth = a1_depth, a2_depth
    out.depth = call.total_depth
    out.a1_purity, out.a2_purity = a1_purity, a2_purity
    out.a1_methylation, out.a2_methylation = a1_meth, a2_meth
    out.a1_motif, out.a2_motif = a1_motif, a2_motif
    if a1 is not None and a2 is None:
        out.zygosity = "Hemizygous"
    elif a1 is not None and a1 == a2:
        out.zygosity = "Homozygous"
    else:
        out.zygosity = "Heterozygous"

    cls = strchive.classify(locus, a1, a2, a1_motif, a2_motif)
    out.status, out.status_color, out.note = cls.status, cls.status_color, cls.note
    return out


def annotate_strchive_panel(
    vcf_path: str | Path,
    strchive_json_path: str | Path,
    sample: Optional[str] = None,
    locus_ids: Optional[set] = None,
) -> list[AnnotatedLocus]:
    """Annotate STRchive loci against a sample's TRGT VCF calls.

    `locus_ids`, when given, scopes this to a specific assay's panel (e.g. PureTarget
    V2's 38 loci) instead of STRchive's full 80 -- so a targeted-panel sample's Overview
    only lists what that assay actually attempts, rather than 42+ inevitable "no call"
    rows for loci it was never going to cover.
    """
    loci = load_strchive_loci(strchive_json_path)
    if locus_ids is not None:
        loci = [l for l in loci if l.id in locus_ids]
    index = TrgtIndex.from_vcf(vcf_path, sample=sample)

    results = []
    for locus in loci:
        call = index.by_trid.get(locus.id) or index.find_overlap(locus.chrom, locus.start, locus.stop)
        results.append(_annotate_locus(locus, call))
    return results


def simulate_pileup_reads(locus: AnnotatedLocus) -> list[dict]:
    """Approximate per-read repeat-unit values for the pileup chart.

    TRGT's VCF reports summary genotype calls (allele length + confidence interval +
    supporting-read count), not each spanning read's individual measurement -- that
    lives in TRGT's separate spanning-reads BAM output, which this pipeline doesn't
    have. To still render a pileup, we distribute each allele's real supporting-read
    count (SD) across its real reported confidence interval (ALLR). This is real
    depth and a real range, just not real individual read values -- callers should
    treat this as an estimate, not raw per-read data.
    """
    if not locus.genotyped:
        return []

    rng = random.Random(_seed_for(locus.id))
    motif_len = max(1, len(locus.motif)) if locus.motif else 1
    reads: list[dict] = []
    for cluster, (units, bp_range, depth) in enumerate(
        [
            (locus.a1, locus.a1_range_bp, locus.a1_depth),
            (locus.a2, locus.a2_range_bp, locus.a2_depth),
        ],
        start=1,
    ):
        if units is None:
            continue
        n = depth if depth and depth > 0 else 6
        if bp_range:
            lo = round(bp_range[0] / motif_len)
            hi = round(bp_range[1] / motif_len)
        else:
            lo = hi = units
        for _ in range(n):
            value = units if lo == hi else rng.randint(min(lo, hi), max(lo, hi))
            reads.append({"value": max(1, value), "cluster": cluster})
    reads.sort(key=lambda r: r["value"])
    return reads


@dataclass
class AnnotatedCatalogRegion:
    chrom: str
    start: int
    stop: int
    length: int
    motifs: str
    struc: str
    genotyped: bool = False
    a1: Optional[int] = None
    a2: Optional[int] = None
    depth: Optional[int] = None


def annotate_catalog_region(region: CatalogRegion, call: Optional[TrgtCall]) -> AnnotatedCatalogRegion:
    out = AnnotatedCatalogRegion(
        chrom=region.chrom, start=region.start, stop=region.stop, length=region.length,
        motifs=region.motifs, struc=region.struc,
    )
    if call is None or not call.genotyped:
        return out
    target_motif = region.motifs.split(",")[0] if region.motifs else None
    a1, a2 = call.repeat_units(target_motif)
    out.genotyped = True
    out.a1, out.a2 = a1, a2
    out.depth = call.total_depth
    return out


def iter_annotated_catalog(
    bed_path: str | Path, vcf_path: str | Path, sample: Optional[str] = None
):
    """Yield AnnotatedCatalogRegion for every region in the genome-wide repeat catalog BED."""
    index = TrgtIndex.from_vcf(vcf_path, sample=sample)
    for region in iter_catalog_regions(bed_path):
        call = index.by_trid.get(region.id) or index.find_overlap(region.chrom, region.start, region.stop)
        yield annotate_catalog_region(region, call)
