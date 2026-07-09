"""gnomAD v4.1 mtDNA constraint metrics: gene- and region-level tolerance to
variation, from https://gnomad.broadinstitute.org (mitochondrial genome
constraint, Feb 2025). Two curated TSVs, committed as fixtures:

    gnomad_v4.1_mtdna_gene_constraint.tsv
        Per (gene, consequence class) observed vs expected sum maximum
        heteroplasmy across gnomAD v4.1's ~56K-genome mtDNA callset, with a
        90% CI on the obs:exp ratio. Consequence is one of synonymous /
        missense / stop_gain (protein genes) or RNA_variant (tRNA/rRNA genes).

    gnomad_v4.1_mtdna_regional_constraint.tsv
        Sub-gene regions with significantly lower obs:exp than their gene's
        average (regional missense constraint for protein genes; regional
        constraint for RNA variants in the two rRNA genes). Start/End
        position are genomic (chrM) coordinates, so matching a variant's
        position needs no codon-frame knowledge.

This is population-level intolerance to variation, not a per-variant
frequency or pathogenicity call -- a low obs:exp means the gene/region as a
whole is depleted of that class of variation in the population, which is
suggestive but not diagnostic for any one variant.

Precisely picking a protein gene's synonymous/missense/stop_gain row for a
given SNV needs codon-level translation against the vertebrate mitochondrial
genetic code (which differs from the standard code, and several mtDNA genes
have overlapping reading frames). This module doesn't do that translation --
it has no ground truth available to validate it against -- so gene-level
lookups return every consequence row for the gene rather than guessing which
one applies to a specific variant.
"""

from __future__ import annotations

import csv
from functools import lru_cache
from pathlib import Path
from typing import Optional, TypedDict

FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"
GENE_CONSTRAINT_TSV = FIXTURES_DIR / "gnomad_v4.1_mtdna_gene_constraint.tsv"
REGIONAL_CONSTRAINT_TSV = FIXTURES_DIR / "gnomad_v4.1_mtdna_regional_constraint.tsv"

SOURCE_LABEL = "gnomAD v4.1 mtDNA constraint metrics"


class GeneConstraintRow(TypedDict):
    symbol: str
    start: int
    end: int
    consequence: str
    observed: float
    expected: float
    obs_exp: float
    lower_ci: float
    upper_ci: float


class RegionalConstraintRow(TypedDict):
    symbol: str
    start: int
    end: int
    protein_residue_start: Optional[int]
    protein_residue_end: Optional[int]
    observed: float
    expected: float
    obs_exp: float
    lower_ci: float
    upper_ci: float


def _num(s: str) -> float:
    return float(s)


def _maybe_int(s: str) -> Optional[int]:
    return None if s == "NA" else int(s)


@lru_cache(maxsize=1)
def _gene_rows() -> list[GeneConstraintRow]:
    rows: list[GeneConstraintRow] = []
    with open(GENE_CONSTRAINT_TSV, encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if not r.get("Symbol"):
                continue
            rows.append(
                {
                    "symbol": r["Symbol"],
                    "start": int(r["Start_position"]),
                    "end": int(r["End_position"]),
                    "consequence": r["Consequence"],
                    "observed": _num(r["Observed"]),
                    "expected": _num(r["Expected"]),
                    "obs_exp": _num(r["obs:exp"]),
                    "lower_ci": _num(r["Lower_CI"]),
                    "upper_ci": _num(r["Upper_CI"]),
                }
            )
    return rows


@lru_cache(maxsize=1)
def _regional_rows() -> list[RegionalConstraintRow]:
    rows: list[RegionalConstraintRow] = []
    with open(REGIONAL_CONSTRAINT_TSV, encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if not r.get("Symbol"):
                continue
            rows.append(
                {
                    "symbol": r["Symbol"],
                    "start": int(r["Start_position"]),
                    "end": int(r["End_position"]),
                    "protein_residue_start": _maybe_int(r["Protein_residue_start"]),
                    "protein_residue_end": _maybe_int(r["Protein_residue_end"]),
                    "observed": _num(r["Observed"]),
                    "expected": _num(r["Expected"]),
                    "obs_exp": _num(r["obs:exp"]),
                    "lower_ci": _num(r["Lower_CI"]),
                    "upper_ci": _num(r["Upper_CI"]),
                }
            )
    return rows


def gene_constraint_for(gene: str) -> list[GeneConstraintRow]:
    """All consequence-class constraint rows for this gene symbol (e.g. all
    of synonymous/missense/stop_gain for a protein gene, or the single
    RNA_variant row for a tRNA/rRNA gene). Empty if the gene isn't in the
    table (e.g. it's outside mtDNA's 37 annotated genes)."""
    return [r for r in _gene_rows() if r["symbol"] == gene]


def regional_constraint_for(pos: int, gene: str) -> Optional[RegionalConstraintRow]:
    """The regional-constraint hotspot containing this genomic position, if
    any -- matched by gene symbol + a plain chrM coordinate range, so no
    codon-frame assumption is needed. None if the position doesn't fall in a
    flagged low-constraint sub-region."""
    for r in _regional_rows():
        if r["symbol"] == gene and r["start"] <= pos <= r["end"]:
            return r
    return None


def constraint_for_variant(pos: int, gene: str) -> Optional[dict]:
    """Constraint context for a called variant's locus, or None if the gene
    isn't covered by the constraint table (e.g. D-loop / intergenic)."""
    gene_rows = gene_constraint_for(gene)
    if not gene_rows:
        return None
    return {
        "source": SOURCE_LABEL,
        "gene_rows": gene_rows,
        "regional": regional_constraint_for(pos, gene),
    }
