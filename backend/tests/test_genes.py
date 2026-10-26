from pathlib import Path

from repeatscope.genes import GeneIndex, GeneRegion, load_genes

REPO_ROOT = Path(__file__).resolve().parents[2]
GENES_BED = REPO_ROOT / "project" / "uploads" / "genes.hg38.bed"


def _idx(*genes: GeneRegion) -> GeneIndex:
    return GeneIndex(list(genes))


def test_overlap_inside_single_gene():
    idx = _idx(GeneRegion("chr1", 1000, 2000, "GENEA", "+"))
    result = idx.annotate("chr1", 1200, 1300)
    assert result.genes == ("GENEA",)
    assert result.nearest_gene is None


def test_overlap_partial_at_either_edge_counts():
    """A region only partially overlapping a gene body (not fully contained)
    still counts as overlapping -- real repeats can straddle a gene boundary."""
    idx = _idx(GeneRegion("chr1", 1000, 2000, "GENEA", "+"))
    assert idx.annotate("chr1", 900, 1100).genes == ("GENEA",)
    assert idx.annotate("chr1", 1900, 2100).genes == ("GENEA",)


def test_touching_but_not_overlapping_is_not_an_overlap():
    """BED half-open convention: a region starting exactly at a gene's stop (or
    ending exactly at its start) is adjacent, not overlapping."""
    idx = _idx(GeneRegion("chr1", 1000, 2000, "GENEA", "+"))
    result = idx.annotate("chr1", 2000, 2100)
    assert result.genes == ()
    assert result.nearest_gene == "GENEA"
    assert result.nearest_distance_bp == 0


def test_overlap_multiple_nested_genes():
    """Real biology, not a data error: antisense/readthrough transcripts really
    do overlap. A region inside both should report both, not just whichever the
    index happens to scan first."""
    idx = _idx(
        GeneRegion("chr1", 1000, 5000, "BIGGENE", "+"),
        GeneRegion("chr1", 2000, 2500, "NESTEDGENE", "-"),
    )
    result = idx.annotate("chr1", 2100, 2200)
    assert set(result.genes) == {"BIGGENE", "NESTEDGENE"}


def test_nearest_gene_plus_strand_direction():
    """Plus-strand gene: upstream is the lower-coordinate side (before its TSS),
    downstream is the higher-coordinate side (past its TES)."""
    idx = _idx(GeneRegion("chr1", 10_000, 20_000, "GENEA", "+"))

    before = idx.annotate("chr1", 9_000, 9_500)
    assert before.nearest_gene == "GENEA"
    assert before.nearest_distance_bp == 500
    assert before.nearest_direction == "upstream"

    after = idx.annotate("chr1", 20_500, 21_000)
    assert after.nearest_gene == "GENEA"
    assert after.nearest_distance_bp == 500
    assert after.nearest_direction == "downstream"


def test_nearest_gene_minus_strand_direction_is_flipped():
    """Minus-strand genes transcribe stop->start, so their TSS/upstream side is
    the *higher* genomic coordinate -- the exact opposite of a plus-strand gene.
    A region past a minus-strand gene's higher-coordinate end is upstream of it,
    not downstream; getting this backwards would mislabel real biology."""
    idx = _idx(GeneRegion("chr1", 10_000, 20_000, "GENEB", "-"))

    before = idx.annotate("chr1", 9_000, 9_500)
    assert before.nearest_direction == "downstream"

    after = idx.annotate("chr1", 20_500, 21_000)
    assert after.nearest_direction == "upstream"


def test_nearest_gene_picks_true_closest_not_largest_start():
    """Once overlap is ruled out, the nearest-before gene is the one with the
    largest *end*, not necessarily the one with the largest start -- an earlier
    -starting but longer gene can still end closer to the query than a later
    -starting, shorter one (real when genes are nested/overlapping elsewhere on
    the same chromosome)."""
    idx = _idx(
        GeneRegion("chr1", 1_000, 9_800, "LONGER_EARLIER", "+"),
        GeneRegion("chr1", 5_000, 5_100, "SHORTER_LATER", "+"),
    )
    result = idx.annotate("chr1", 10_000, 10_100)
    assert result.nearest_gene == "LONGER_EARLIER"
    assert result.nearest_distance_bp == 200


def test_no_genes_on_chromosome_returns_empty():
    idx = _idx(GeneRegion("chr1", 1000, 2000, "GENEA", "+"))
    result = idx.annotate("chrZ_doesnt_exist", 1, 100)
    assert result.genes == ()
    assert result.nearest_gene is None


def test_no_genes_at_all_on_edge_of_chromosome_but_none_before_or_after():
    """A single gene on a chromosome, queried from a position where only one
    side (before or after) has a candidate -- the other side must not crash."""
    idx = _idx(GeneRegion("chr1", 1_000_000, 1_000_100, "ONLYGENE", "+"))
    before_everything = idx.annotate("chr1", 0, 100)
    assert before_everything.nearest_gene == "ONLYGENE"
    assert before_everything.nearest_direction == "upstream"


# --- against the real gene reference file -----------------------------------


def test_real_htt_repeat_overlaps_htt_gene_body():
    """HD_HTT's real CAG repeat locus (chr4:3,074,876-3,074,933) falls well
    within RefSeq's own HTT gene body (chr4:3,074,680-3,243,960)."""
    idx = GeneIndex(load_genes(GENES_BED))
    result = idx.annotate("chr4", 3_074_876, 3_074_933)
    assert result.genes == ("HTT",)


def test_real_bclaf3_repeat_is_intergenic_not_inside_refseqs_bclaf3():
    """Regression/documentation test: FRAXG_BCLAF3's real repeat locus
    (chrX:19,990,922-19,990,973) sits ~20.7kb upstream of RefSeq's own curated
    BCLAF3 transcript body (chrX:19,912,861-19,970,264), not inside it -- the
    paper's "5' UTR" description doesn't match this specific RefSeq transcript
    model's annotated boundaries. The true nearest gene is a different,
    smaller, closer locus (LOC729609, ~1.1kb away) -- reported honestly rather
    than special-cased to say BCLAF3 just because that's the more familiar name.
    LOC729609 is minus-strand, so its TSS/upstream side is its higher-coordinate
    end (19,989,779) -- our repeat sits even further right than that, i.e.
    upstream of it, not downstream."""
    idx = GeneIndex(load_genes(GENES_BED))
    result = idx.annotate("chrX", 19_990_922, 19_990_973)
    assert result.genes == ()
    assert result.nearest_gene == "LOC729609"
    assert result.nearest_distance_bp == 1143
    assert result.nearest_direction == "upstream"


def test_real_gene_symbol_at_distant_paralogous_copies_not_merged():
    """MIR6859-1 is a real small RNA gene duplicated at multiple, far-apart loci
    on chr1 alone (17,368-17,436 and 187,890-187,958) -- a naive min/max-by
    -symbol collapse would bridge them into one bogus ~170kb mega-span. Confirm
    the reference file keeps them as separate rows: a region between the two
    copies should NOT be reported as overlapping (or spuriously close to) a
    170kb-wide "MIR6859-1", since no such single gene body actually exists."""
    idx = GeneIndex(load_genes(GENES_BED))
    result = idx.annotate("chr1", 100_000, 100_100)
    assert "MIR6859-1" not in result.genes
    if result.nearest_gene == "MIR6859-1":
        assert result.nearest_distance_bp < 100_000
