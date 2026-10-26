"""Tests for RExPRT genome-wide pathogenicity scores -- real, published data
(Danzi et al. 2024, PMC10832122) joined onto this project's own genome-wide
repeat catalog. See project/uploads/README.md for full provenance."""

from pathlib import Path

from repeatscope.annotate import annotate_catalog_region
from repeatscope.catalog import CatalogRegion
from repeatscope.rexprt import load_rexprt_scores

REPO_ROOT = Path(__file__).resolve().parents[2]
REXPRT_TSV = REPO_ROOT / "project" / "uploads" / "rexprt_scores.hg38.tsv"


def test_loads_real_scores_file():
    scores = load_rexprt_scores(REXPRT_TSV)
    assert len(scores) == 727_004


def test_real_fmr1_locus_is_predicted_pathogenic():
    """The FMR1 CGG repeat's real catalog region (chrX:147,912,049-147,912,198,
    TRExplorer v2.0's variation-clusters-and-isolated-TRs catalog -- this
    project's genome-wide TRGT run was genotyped against this exact catalog,
    giving 100% real TRID-matched genotype coverage) is the same repeat RExPRT
    scored in its own training set -- correctly classified pathogenic
    (ensemble_binary=1), in the gene's first exon/5'UTR, a real
    OMIM/nervous-system-expressed, high-GERP, high-constraint locus."""
    scores = load_rexprt_scores(REXPRT_TSV)
    score = scores[("chrX", 147_912_049, 147_912_198)]
    assert score.gene == "FMR1"
    assert score.motif == "CGG"
    assert score.ensemble_binary is True
    assert score.region == "exon"
    assert score.utr_5 is True
    assert score.promoter is True
    assert score.gerp > 1.0


def test_real_c9orf72_locus_is_a_documented_false_negative():
    """Known, verified model limitation (see project/uploads/README.md): the
    real C9orf72 GGGGCC expansion region (STRchive FTDALS1_C9orf72 spans
    chr9:27,573,484-27,573,546; this catalog's variation-cluster consolidation
    combines that same repeat with adjacent sub-repeats into
    chr9:27,573,523-27,573,546, carrying the identical underlying RExPRT
    source row as every other coordinate framing of this locus tested so far)
    scores *below* RExPRT's 0.5 threshold at reference length --
    ensemble_binary=0, "benign" -- despite being one of the best-established
    pathogenic repeat expansions known. This must never be surfaced as
    "C9orf72 is benign"; it documents why RExPRT scores are a research
    candidate-prioritization signal, not a diagnostic override of the
    STRchive disease panel's own curated classification."""
    scores = load_rexprt_scores(REXPRT_TSV)
    score = scores[("chr9", 27_573_523, 27_573_546)]
    assert score.gene == "C9orf72"
    assert score.ensemble_binary is False
    assert score.ensemble_max < 0.5


def test_unmatched_region_is_absent_not_a_fabricated_default():
    scores = load_rexprt_scores(REXPRT_TSV)
    assert ("chr1", 1, 2) not in scores


def test_annotate_catalog_region_attaches_real_pathogenicity():
    region = CatalogRegion(
        chrom="chrX", start=147_912_049, stop=147_912_198,
        id="FXS_FMR1", motifs="CGG", struc="<TR>",
    )
    scores = load_rexprt_scores(REXPRT_TSV)
    out = annotate_catalog_region(region, call=None, rexprt_scores=scores)
    assert out.pathogenicity is not None
    assert out.pathogenicity.gene == "FMR1"
    assert out.pathogenicity.ensemble_binary is True


def test_annotate_catalog_region_without_match_leaves_pathogenicity_none():
    region = CatalogRegion(
        chrom="chr1", start=1, stop=2, id="x", motifs="A", struc="(A)n",
    )
    scores = load_rexprt_scores(REXPRT_TSV)
    out = annotate_catalog_region(region, call=None, rexprt_scores=scores)
    assert out.pathogenicity is None


def test_annotate_catalog_region_with_no_rexprt_scores_given_leaves_pathogenicity_none():
    """No rexprt_scores argument at all (e.g. the fixture file isn't present)
    must not error and must not fabricate a score."""
    region = CatalogRegion(
        chrom="chrX", start=147_912_049, stop=147_912_198,
        id="FXS_FMR1", motifs="CGG", struc="<TR>",
    )
    out = annotate_catalog_region(region, call=None)
    assert out.pathogenicity is None
