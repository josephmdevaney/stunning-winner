"""Chromosome names, SV types and size bins shared across the structural-variant suite."""

from __future__ import annotations

from typing import Optional

CHROM_ORDER = {f"chr{i}": i for i in range(1, 23)}
CHROM_ORDER.update({"chrX": 23, "chrY": 24, "chrM": 25})

SV_TYPES = ("DEL", "DUP", "INS", "INV", "CNV", "BND")

# (upper bound in bp, label); the last bin is open. Callers differ in their smallest SV (sawfish keeps deletions from 15 bp and
# insertions from 35 bp), so calls under 50 bp get their own bin instead of being counted as 50–99.
SIZE_BINS = ((50, "< 50"), (100, "50–99"), (300, "100–299"), (1_000, "300–999"), (10_000, "1–10 kb"), (100_000, "10–100 kb"), (1_000_000, "100 kb–1 Mb"))
LAST_BIN = "> 1 Mb"


def canon_chrom(chrom: str) -> str:
    c = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if c in ("M", "MT", "m", "mt"):
        return "chrM"
    return "chr" + c


def chrom_order(chrom: str) -> int:
    return CHROM_ORDER.get(canon_chrom(chrom), 99)


def size_bin(size: Optional[int]) -> Optional[str]:
    if size is None:
        return None
    for limit, label in SIZE_BINS:
        if size < limit:
            return label
    return LAST_BIN
