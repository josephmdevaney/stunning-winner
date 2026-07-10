"""gnomAD v3.1 mtDNA population frequency: an exact (pos, ref, alt) lookup
against gnomAD's real mitochondrial callset (56,434 genomes), from a reduced
TSV gnomAD itself publishes (release/3.1/vcf/genomes/gnomad.genomes.v3.1.
sites.chrM.reduced_annotations.tsv) alongside their full sites VCF.

gnomAD reports two allele frequencies per site because heteroplasmy makes
"the" allele frequency ambiguous for mtDNA:

    AF_hom  frequency restricted to calls with heteroplasmy >= 0.95
    AF_het  frequency restricted to calls with heteroplasmy in [0.10, 0.95)

`filters` is gnomAD's own site-level QC flag (PASS, or e.g.
artifact_prone_site / indel_stack / npg for known-troublesome mtDNA
positions such as homopolymer tracts) -- surfaced verbatim rather than
silently trusting every matched row.

A variant absent from this table was never observed in gnomAD's 56,434
genomes; that's reported as such, not silently treated the same as a very
low (but present) frequency.

This lookup only ever runs against real MitoRSaw-derived calls (see vcf.py).
The synthetic demo sample's positions are randomly generated within real
gene boundaries and must never be checked against a real population
database -- a coincidental match would misrepresent a fabricated call as
having a real population frequency.
"""

from __future__ import annotations

import csv
from functools import lru_cache
from pathlib import Path
from typing import Optional, TypedDict

FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"
FREQ_TSV = FIXTURES_DIR / "gnomad.genomes.v3.1.sites.chrM.reduced_annotations.tsv"

SOURCE_LABEL = "gnomAD v3.1 mtDNA callset"
COHORT_SIZE = 56434  # genomes in the gnomAD v3.1 mitochondrial release


class GnomadHit(TypedDict):
    af_hom: float
    af_het: float
    ac_hom: int
    ac_het: int
    an: int
    filters: list[str]
    pass_qc: bool
    max_observed_heteroplasmy: float


@lru_cache(maxsize=1)
def _table() -> dict[tuple[int, str, str], GnomadHit]:
    table: dict[tuple[int, str, str], GnomadHit] = {}
    with open(FREQ_TSV, encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            key = (int(r["position"]), r["ref"], r["alt"])
            filters = r["filters"].split(",")
            table[key] = {
                "af_hom": float(r["AF_hom"]),
                "af_het": float(r["AF_het"]),
                "ac_hom": int(r["AC_hom"]),
                "ac_het": int(r["AC_het"]),
                "an": int(r["AN"]),
                "filters": filters,
                "pass_qc": filters == ["PASS"],
                "max_observed_heteroplasmy": float(r["max_observed_heteroplasmy"]),
            }
    return table


def frequency_for(pos: int, ref: str, alt: str) -> Optional[GnomadHit]:
    """The gnomAD v3.1 record for this exact (pos, ref, alt), or None if the
    variant was never observed in the callset."""
    return _table().get((pos, ref, alt))


# ACMG/AMP frequency-based evidence thresholds, per McCormick et al. 2020 ("Specifications
# of the ACMG/AMP standards and guidelines for mitochondrial DNA variant interpretation",
# Hum Mutat 41(12):2028-2057) -- the ClinGen mtDNA Expert Panel's own specification of
# BA1/BS1/PM2 for mtDNA, derived from their analysis of real MITOMAP pathogenic-variant
# frequency distributions. The thresholds are deliberately non-contiguous: the 0.0099-0.01
# and 0.00002-0.005 gaps are real findings from that analysis, not omissions here.
ACMG_BA1_THRESHOLD = 0.01
ACMG_BS1_RANGE = (0.005, 0.0099)
ACMG_PM2_SUPPORTING_THRESHOLD = 0.00002


def acmg_frequency_evidence(hit: Optional[GnomadHit]) -> list[dict]:
    """BA1 (stand-alone benign) / BS1 (strong benign) / PM2_supporting (supporting
    pathogenic) against gnomAD's real mtDNA callset. Total population frequency is
    AF_hom + AF_het -- the variant's real occurrence at any heteroplasmy level
    >=10%, not just the homoplasmic slice -- since McCormick et al. specify a plain
    "allele frequency" without splitting by heteroplasmy band.

    Only one code (or none, in the gaps between thresholds) is ever returned, since
    the three ranges are mutually exclusive by construction."""
    if hit is None:
        return [
            {
                "code": "PM2_supporting",
                "strength": "Supporting (pathogenic)",
                "note": f"Absent from {SOURCE_LABEL} (n={COHORT_SIZE:,}); frequency below the 0.002% threshold.",
            }
        ]
    af = hit["af_hom"] + hit["af_het"]
    if af > ACMG_BA1_THRESHOLD:
        return [
            {
                "code": "BA1",
                "strength": "Stand-alone (benign)",
                "note": f"Frequency {af:.2%} in {SOURCE_LABEL} exceeds the 1% stand-alone benign threshold.",
            }
        ]
    lo, hi = ACMG_BS1_RANGE
    if lo <= af <= hi:
        return [
            {
                "code": "BS1",
                "strength": "Strong (benign)",
                "note": f"Frequency {af:.2%} in {SOURCE_LABEL} falls in the 0.5-0.99% strong-benign range.",
            }
        ]
    if af < ACMG_PM2_SUPPORTING_THRESHOLD:
        return [
            {
                "code": "PM2_supporting",
                "strength": "Supporting (pathogenic)",
                "note": f"Frequency {af:.4%} in {SOURCE_LABEL} is below the 0.002% supporting-pathogenic threshold.",
            }
        ]
    return []
