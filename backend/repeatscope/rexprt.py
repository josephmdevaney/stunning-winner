"""RExPRT genome-wide tandem-repeat pathogenicity predictions for the genome-wide
repeat catalog (not the 82-locus STRchive disease panel, which has its own
curated, length-threshold-based clinical classification and is unaffected by
this module).

RExPRT (Danzi et al. 2024, PMC10832122, github.com/ZuchnerLab/RExPRT) is a
published ensemble SVM+XGBoost model scoring a tandem repeat's pathogenicity
potential from its reference-length genomic context (conservation, gene
constraint, regulatory overlap, expression). These are the lab's own real
precomputed scores for ~836k reference TRs, lifted from hg19 to hg38 and
joined onto this project's own catalog coordinates by
scripts/real_data/build_rexprt_catalog_scores.py -- see
project/uploads/README.md for full provenance, coverage stats, and a verified
real false-negative example (C9orf72) that any use of these scores must
respect: this is a research-grade candidate-prioritization signal, not a
diagnostic call, and a missing entry means "not evaluated", never "benign".
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

_BOOL_COLS = {
    "promoter", "utr_5", "utr_3", "tad", "estr", "open_reg_region",
    "rad21", "smc3", "ctcf", "ensemble_binary",
}
_FLOAT_COLS = {"gerp", "loeuf", "pli", "gc_content", "svm", "xgb", "ensemble_score", "ensemble_max"}
_INT_COLS = {"alu_distance", "hg19_start", "hg19_end"}


@dataclass(frozen=True)
class RexprtScore:
    motif: str
    gene: str
    gene_type: str
    region: str  # "exon" | "intron" | "intergenic"
    location: str  # "First" | "Middle" | "Last" (position among the gene's own exons/introns)
    gerp: float
    loeuf: float
    pli: float
    promoter: bool
    utr_5: bool
    utr_3: bool
    tad: bool
    estr: bool
    open_reg_region: bool
    rad21: bool
    smc3: bool
    ctcf: bool
    alu: str
    alu_distance: int
    gc_content: float
    tissue: str
    svm: float
    xgb: float
    ensemble_score: float
    ensemble_max: float
    ensemble_binary: bool
    hg19_chrom: str
    hg19_start: int
    hg19_end: int


def load_rexprt_scores(path: str | Path) -> dict[tuple[str, int, int], RexprtScore]:
    """Load the derived TSV, keyed by this project's own catalog region
    coordinates (chrom, start, stop) -- an exact-match lookup, since the join
    script already resolved each RExPRT source row to its best-overlapping
    catalog region at build time."""
    scores: dict[tuple[str, int, int], RexprtScore] = {}
    with open(path, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            fields = line.rstrip("\n").split("\t")
            row = dict(zip(header, fields))
            chrom, start, stop = row["chrom"], int(row["start"]), int(row["stop"])
            kwargs = {}
            for key, value in row.items():
                if key in ("chrom", "start", "stop"):
                    continue
                if key in _BOOL_COLS:
                    kwargs[key] = value == "1"
                elif key in _FLOAT_COLS:
                    kwargs[key] = float(value)
                elif key in _INT_COLS:
                    kwargs[key] = int(value)
                else:
                    kwargs[key] = value
            scores[(chrom, start, stop)] = RexprtScore(**kwargs)
    return scores
