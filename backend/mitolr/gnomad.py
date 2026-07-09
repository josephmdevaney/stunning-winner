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
