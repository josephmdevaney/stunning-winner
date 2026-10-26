from repeatscope.catalog import CatalogIndex, CatalogRegion, chrom_sort_key


def _region(chrom, start, stop, motifs="CAG", region_id=None):
    return CatalogRegion(
        chrom=chrom, start=start, stop=stop, id=region_id or f"{chrom}_{start}_{stop}",
        motifs=motifs, struc=f"({motifs})n",
    )


def test_chrom_sort_key_orders_autosomes_then_sex_then_mt():
    chroms = ["chr2", "chrM", "chr10", "chrX", "chr1", "chrY"]
    assert sorted(chroms, key=chrom_sort_key) == ["chr1", "chr2", "chr10", "chrX", "chrY", "chrM"]


def test_index_totals_and_unique_motifs():
    idx = CatalogIndex([
        _region("chr1", 100, 200, motifs="CAG"),
        _region("chr1", 300, 400, motifs="CAG"),
        _region("chr2", 50, 150, motifs="GCC"),
    ])
    assert idx.total == 3
    assert idx.unique_motifs == 2
    assert idx.chromosomes == ["chr1", "chr2"]


def test_index_all_is_sorted_by_chrom_then_start_regardless_of_input_order():
    idx = CatalogIndex([
        _region("chr2", 500, 600),
        _region("chr1", 300, 400),
        _region("chr1", 100, 200),
    ])
    assert [(r.chrom, r.start) for r in idx.all] == [("chr1", 100), ("chr1", 300), ("chr2", 500)]


def test_filter_by_chrom():
    idx = CatalogIndex([_region("chr1", 100, 200), _region("chr2", 100, 200)])
    assert [r.chrom for r in idx.filter(chrom="chr2")] == ["chr2"]
    assert idx.filter(chrom="chrZ_missing") == []


def test_filter_by_motif_substring_case_insensitive():
    idx = CatalogIndex([
        _region("chr1", 100, 200, motifs="CAG"),
        _region("chr1", 300, 400, motifs="GCC"),
    ])
    assert len(idx.filter(motif="cag")) == 1
    assert idx.filter(motif="cag")[0].motifs == "CAG"


def test_filter_by_min_len():
    idx = CatalogIndex([
        _region("chr1", 100, 110),  # length 10
        _region("chr1", 200, 250),  # length 50
    ])
    result = idx.filter(min_len=20)
    assert len(result) == 1
    assert result[0].length == 50


def test_filter_combines_all_criteria():
    idx = CatalogIndex([
        _region("chr1", 100, 110, motifs="CAG"),
        _region("chr1", 200, 260, motifs="CAG"),
        _region("chr2", 200, 260, motifs="CAG"),
    ])
    result = idx.filter(chrom="chr1", motif="CAG", min_len=20)
    assert [(r.chrom, r.start) for r in result] == [("chr1", 200)]


def test_from_bed_real_catalog_file_matches_expected_scale():
    """Sanity check against the real genome-wide catalog file, not a synthetic
    fixture -- confirms the loader/index actually works end-to-end. Currently
    scoped to RExPRT-predicted-pathogenic regions only (see
    project/uploads/README.md), so this is thousands, not millions, of rows --
    that's the point, not a regression."""
    from pathlib import Path

    bed_path = Path(__file__).resolve().parents[2] / "project" / "uploads" / "repeat_catalog.hg38.bed"
    idx = CatalogIndex.from_bed(bed_path)
    assert idx.total > 10_000
    assert set(idx.chromosomes) >= {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
