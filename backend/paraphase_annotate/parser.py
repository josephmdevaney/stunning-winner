"""Parses interpret_paraphase.py's own two output CSVs (column names are that
script's, not ours -- see its SUMMARY_HEADER/DETAILED_HEADER) into the shapes
the API returns."""

from __future__ import annotations

import csv
from pathlib import Path

from . import enrich


# interpret_paraphase.py joins the values of a cell with commas but does not quote it, so the csv module reads every
# value after the first as an extra column. Those spill over under this key and `_last_cell` puts them back; without it
# a haplotype's pathogenic variants beyond the first are silently dropped.
_SPILL = "_spilled_cells"


def _read_csv(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, restkey=_SPILL))


def _last_cell(row: dict, column: str) -> str:
    """The value of a row's final column, with any comma-separated values that spilled past it joined back on."""
    spilled = row.get(_SPILL)
    return ",".join([row[column], *spilled]) if spilled else row[column]


def _int_or_none(v: str) -> int | None:
    return int(v) if v not in ("", "NA", None) else None


def _number_or_none(v: str) -> int | float | None:
    """A depth, which Paraphase writes as `50` for some samples and `39.0` or `39.5` for others: a whole number comes
    back as an int, anything else as a float."""
    if v in ("", "NA", None):
        return None
    number = float(v)
    return int(number) if number.is_integer() else number


def parse_summary(path: Path) -> list[dict]:
    rows = []
    for r in _read_csv(path):
        with_patho = _int_or_none(r["Num_haplotypes_with_pathogenic_variants"])
        without_patho = _int_or_none(r["Num_haplotypes_without_pathogenic_variants"])
        rows.append({
            "region": r["Region"],
            "region_median_depth": _number_or_none(r["Region_median_depth"]),
            "gene": r["Gene"] if r["Gene"] != "NA" else None,
            "total_haplotypes": _int_or_none(r["Total_haplotypes"]),
            "n_haplotypes_with_pathogenic": with_patho,
            "n_haplotypes_without_pathogenic": without_patho,
            # PacBio's own bar for clinical significance in a recessive-disease
            # model: zero haplotypes free of a pathogenic variant is the
            # long-read analog of biallelic mutations (see annotation/README.md).
            "zero_functional_haplotypes": without_patho == 0,
            "info": _last_cell(r, "Useful_information_from_Paraphase_JSON") or None,
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
        cell = _last_cell(r, "Pathogenic_variants_annotated")
        variants = enrich.parse_pathogenic_variants(cell)
        rows.append({
            "region": r["Region"],
            "gene": r["Gene"],
            "transcript": r["Transcript"],
            "transcript_start": _int_or_none(r["Transcript_start"]),
            "transcript_end": _int_or_none(r["Transcript_end"]),
            "haplotype": r["Haplotype"],
            "haplotype_boundary": r["Haplotype_boundary"],
            "boundary_check": r["Boundary_check"],
            "known_nonfunctional": cell == "known_nonfunctional",
            "pathogenic_variants": variants,
        })
    return rows
