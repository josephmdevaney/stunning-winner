"""Parses interpret_paraphase.py's own two output CSVs (column names are that
script's, not ours -- see its SUMMARY_HEADER/DETAILED_HEADER) into the shapes
the API returns."""

from __future__ import annotations

import csv
from pathlib import Path

from . import enrich


def _read_csv(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _int_or_none(v: str) -> int | None:
    return int(v) if v not in ("", "NA", None) else None


def parse_summary(path: Path) -> list[dict]:
    rows = []
    for r in _read_csv(path):
        with_patho = _int_or_none(r["Num_haplotypes_with_pathogenic_variants"])
        without_patho = _int_or_none(r["Num_haplotypes_without_pathogenic_variants"])
        rows.append({
            "region": r["Region"],
            "region_median_depth": _int_or_none(r["Region_median_depth"]),
            "gene": r["Gene"] if r["Gene"] != "NA" else None,
            "total_haplotypes": _int_or_none(r["Total_haplotypes"]),
            "n_haplotypes_with_pathogenic": with_patho,
            "n_haplotypes_without_pathogenic": without_patho,
            # PacBio's own bar for clinical significance in a recessive-disease
            # model: zero haplotypes free of a pathogenic variant is the
            # long-read analog of biallelic mutations (see annotation/README.md).
            "zero_functional_haplotypes": without_patho == 0,
            "info": r["Useful_information_from_Paraphase_JSON"] or None,
        })
    return rows


def parse_detailed(path: Path, region: str | None = None) -> list[dict]:
    rows = []
    for r in _read_csv(path):
        if region is not None and r["Region"] != region:
            continue
        # Adds c./p. notation + a verified ClinVar classification via a
        # separate, fast re-query (see enrich.py's own docstring for why this
        # isn't done inside interpret_paraphase.py's own much slower run).
        variants = enrich.parse_pathogenic_variants(r["Pathogenic_variants_annotated"])
        rows.append({
            "region": r["Region"],
            "gene": r["Gene"],
            "transcript": r["Transcript"],
            "transcript_start": _int_or_none(r["Transcript_start"]),
            "transcript_end": _int_or_none(r["Transcript_end"]),
            "haplotype": r["Haplotype"],
            "haplotype_boundary": r["Haplotype_boundary"],
            "boundary_check": r["Boundary_check"],
            "known_nonfunctional": r["Pathogenic_variants_annotated"] == "known_nonfunctional",
            "pathogenic_variants": variants,
        })
    return rows
