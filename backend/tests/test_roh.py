"""Tests for the ROH (Runs of Homozygosity) module: real bcftools roh RG-line
parsing, region classification, gene overlap, F-value/UPD derivation, and the
API routes -- all against the real HG002 fixture."""

from fastapi.testclient import TestClient

from repeatscope import roh
from repeatscope.api.main import app

client = TestClient(app)

HG002_ROH = roh.ROH_SAMPLES["HG002"]


# --- parsing ------------------------------------------------------------------

def test_parses_real_rg_lines():
    regions = roh.parse_roh_file(HG002_ROH)
    assert len(regions) == 8966
    assert all(r.sample == "HG002" for r in regions)
    chroms = {r.chrom for r in regions}
    assert chroms == {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}


def test_handles_stray_empty_field():
    # fixtures/roh/HG002.GRCh38.roh.txt:chr6 region has a genuine paste artifact
    # (an extra empty tab). The real length_bp (43767) must still be recovered,
    # and it must equal end - start + 1 like every other real RG row.
    regions = roh.parse_roh_file(HG002_ROH)
    hit = next(r for r in regions if r.chrom == "chr6" and r.start == 137568592)
    assert hit.length_bp == 43767
    assert hit.length_bp == hit.end - hit.start + 1


def test_ignores_comment_lines():
    regions = roh.parse_roh_file(HG002_ROH)
    assert all(r.chrom.startswith("chr") for r in regions)


# --- classification -------------------------------------------------------

def test_position_classification():
    assert roh.classify_position("chr1", 1, 500_000) == "telomeric"
    assert roh.classify_position("chr1", 248_900_000, 248_956_422) == "telomeric"
    assert roh.classify_position("chr1", 120_000_000, 125_000_000) == "pericentromeric"
    assert roh.classify_position("chr1", 50_000_000, 51_000_000) == "interstitial"
    assert roh.classify_position("chrM", 1, 100) is None


def test_length_band():
    assert roh.length_band(500_000) == "<1 Mb"
    assert roh.length_band(1_500_000) == "1–3 Mb"
    assert roh.length_band(3_500_000) == "3–5 Mb"
    assert roh.length_band(6_000_000) == "≥5 Mb"


def test_gene_overlap_real_hfe_region():
    hits = roh.gene_overlaps("chr6", 26_000_000, 26_200_000)
    assert [g["sym"] for g in hits] == ["HFE"]
    assert roh.gene_overlaps("chr6", 1, 100) == []


def test_upd_flag_only_applies_to_autosomes():
    assert roh.upd_flag("chrX", "telomeric", 10_000_000) is None
    assert roh.upd_flag("chr1", "telomeric", 100) is None  # too short
    assert roh.upd_flag("chr7", "telomeric", 5_000_000) is not None  # imprinted, at threshold
    assert roh.upd_flag("chr1", "interstitial", 15_000_000) is not None  # non-imprinted threshold
    assert roh.upd_flag("chr1", "interstitial", 14_000_000) is None  # below non-imprinted threshold


# --- overview / F-value / real findings ------------------------------------

def test_overview_stats_are_internally_consistent():
    ov = roh.build_overview("HG002", HG002_ROH)
    stats = ov["stats"]
    assert stats["n_regions_total"] == 8966
    assert (
        stats["n_regions_autosomal"] + stats["n_regions_chrX"]
        + stats["n_regions_chrY"] + stats["n_regions_chrM"]
        == stats["n_regions_total"]
    )
    assert stats["autosomal_genome_length_bp"] == roh_reference_autosomal_length()
    assert 0 < stats["f_value"] < 0.1  # HG002 is not consanguineous; low F expected
    assert stats["f_value_pct"] == round(stats["f_value"] * 100, 3)


def roh_reference_autosomal_length():
    from repeatscope import roh_reference as ref
    return ref.AUTOSOMAL_GENOME_LENGTH


def test_overview_finds_the_real_chr9_upd_flag():
    ov = roh.build_overview("HG002", HG002_ROH)
    flags = ov["upd_flags"]
    assert len(flags) == 1
    f = flags[0]
    assert f["chrom"] == "chr9"
    assert f["start"] == 43832600 and f["end"] == 60559564
    assert "interstitial" in f["upd_flag"]


def test_ar_candidate_genes_excludes_small_background_roh():
    # TJP2 (chr9) is overlapped by a real but tiny (~49 kb) ROH -- far below the
    # 1 Mb AR-candidate floor (AR_CANDIDATE_MIN_LENGTH_BP), so it must not be
    # surfaced as a candidate.
    ov = roh.build_overview("HG002", HG002_ROH)
    symbols = {g["symbol"] for g in ov["ar_candidate_genes"]}
    assert "TJP2" not in symbols


def test_ar_candidate_min_length_is_separate_from_f_value_floor():
    # AR candidate genes use a deliberately lower, separate floor (1 Mb) than
    # the F-value's own IBD-informative floor (3 Mb) -- a real AR
    # disease-causing variant can sit in a smaller genuinely-autozygous ROH
    # than the F-value estimate's stricter requirement.
    ov = roh.build_overview("HG002", HG002_ROH)
    assert ov["stats"]["ar_candidate_min_length_mb"] == 1
    assert ov["stats"]["f_value_min_length_mb"] == 3


def test_relationship_f_values_present_and_ordered_by_closeness():
    ov = roh.build_overview("HG002", HG002_ROH)
    vals = [r["f_value"] for r in ov["relationship_f_values"]]
    assert vals == sorted(vals, reverse=True)
    assert vals[0] == 0.25  # parent-child / full sibs


def test_by_chromosome_covers_all_autosomes_plus_sex_and_mt():
    ov = roh.build_overview("HG002", HG002_ROH)
    chroms = {c["chrom"] for c in ov["by_chromosome"]}
    assert chroms == {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
    chr9 = next(c for c in ov["by_chromosome"] if c["chrom"] == "chr9")
    assert chr9["imprinted"] is False
    chr7 = next(c for c in ov["by_chromosome"] if c["chrom"] == "chr7")
    assert chr7["imprinted"] is True


# --- query_regions filtering/pagination ------------------------------------

def test_query_regions_filters_by_chrom_and_length():
    result = roh.query_regions("HG002", chrom="chr9", min_len_mb=3)
    assert result["total_matching"] == 1
    assert result["rows"][0]["chrom"] == "chr9"


def test_query_regions_paginates():
    page0 = roh.query_regions("HG002", page=0, page_size=10)
    page1 = roh.query_regions("HG002", page=1, page_size=10)
    assert len(page0["rows"]) == 10
    assert page0["rows"][0]["chrom"] != page1["rows"][0]["chrom"] or (
        page0["rows"][0]["start"] != page1["rows"][0]["start"]
    )
    assert page0["total_all"] == page1["total_all"] == 8966


# --- API routes -------------------------------------------------------------

def test_roh_samples_endpoint():
    body = client.get("/api/roh/samples").json()
    assert body == {"samples": ["HG002"], "default": "HG002"}


def test_roh_overview_endpoint_defaults_to_hg002():
    body = client.get("/api/roh/overview").json()
    assert body["sample"] == "HG002"
    assert body["stats"]["n_regions_total"] == 8966


def test_roh_overview_unknown_sample_404s():
    assert client.get("/api/roh/overview?sample=bogus").status_code == 404


def test_roh_regions_endpoint_filters():
    body = client.get("/api/roh/regions?chrom=chr9&min_len_mb=3").json()
    assert body["total_matching"] == 1
    assert body["rows"][0]["chrom"] == "chr9"


def test_roh_regions_endpoint_unknown_sample_404s():
    assert client.get("/api/roh/regions?sample=bogus").status_code == 404


def test_repeatscope_routes_still_present():
    assert client.get("/api/modes").status_code == 200
