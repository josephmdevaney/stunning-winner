"""Chromosome names and allele normalisation shared by every join in this suite.

ClinVar, gnomAD and the sample VCF can write the same variant differently (``chr7`` vs
``7``; a padding base kept or dropped). Every lookup goes through ``variant_key`` so that
both sides are compared in one trimmed form."""

from __future__ import annotations

from typing import Tuple

CHROM_ORDER = {f"chr{i}": i for i in range(1, 23)}
CHROM_ORDER.update({"chrX": 23, "chrY": 24, "chrM": 25})


def canon_chrom(chrom: str) -> str:
    """``7`` / ``chr7`` -> ``chr7``; ``MT`` / ``chrMT`` / ``M`` -> ``chrM``."""
    c = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if c in ("M", "MT", "m", "mt"):
        return "chrM"
    return "chr" + c


def chrom_order(chrom: str) -> int:
    return CHROM_ORDER.get(canon_chrom(chrom), 99)


def trim_alleles(pos: int, ref: str, alt: str) -> Tuple[int, str, str]:
    """Strip the bases shared by ``ref`` and ``alt``: common suffix first, then common
    prefix (moving ``pos`` along with it). Either allele can come out empty, which is the
    parsimonious form of a pure insertion or deletion."""
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while ref and alt and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def variant_key(chrom: str, pos: int, ref: str, alt: str) -> Tuple[str, int, str, str]:
    pos, ref, alt = trim_alleles(pos, ref.upper(), alt.upper())
    return canon_chrom(chrom), pos, ref, alt


def variant_kind(ref: str, alt: str) -> str:
    """snv | mnv | ins | del | complex, from the trimmed alleles."""
    if not ref and alt:
        return "ins"
    if ref and not alt:
        return "del"
    if len(ref) == 1 and len(alt) == 1:
        return "snv"
    if len(ref) == len(alt):
        return "mnv"
    return "complex"
