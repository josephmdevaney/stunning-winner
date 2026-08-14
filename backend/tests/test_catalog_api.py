from fastapi.testclient import TestClient

from repeatscope.api.main import app


def test_catalog_rows_include_real_rexprt_pathogenicity():
    # Filter by rexprt_pathogenic_only rather than a motif substring: this
    # catalog's variation-cluster consolidation labels FMR1's merged region
    # with just one constituent sub-repeat's motif ("TGG"), not "CGG", even
    # though the real repeat and RExPRT's own independent score are CGG-based
    # -- a motif-substring filter against the catalog's own MOTIFS field would
    # miss it. chrX has 1006 pathogenic-predicted hits; page_size covers all
    # of them rather than relying on FMR1's own sort position.
    with TestClient(app) as client:
        resp = client.get(
            "/api/catalog",
            params={"mode": "genome_wide", "chrom": "chrX", "rexprt_pathogenic_only": "true", "page_size": 2000},
        )
    assert resp.status_code == 200
    data = resp.json()
    assert data["total_matching"] == 1006
    fmr1 = next(r for r in data["rows"] if r["start"] == 147_912_049 and r["stop"] == 147_912_198)
    assert fmr1["genotyped"] is True
    assert fmr1["genes"] == ["FMR1"]
    assert fmr1["pathogenicity"]["gene"] == "FMR1"
    assert fmr1["pathogenicity"]["predicted_pathogenic"] is True


def test_catalog_summary_counts_are_internally_consistent():
    with TestClient(app) as client:
        resp = client.get("/api/catalog", params={"mode": "genome_wide", "page_size": 1})
    data = resp.json()
    assert 0 < data["n_rexprt_predicted_pathogenic"] <= data["n_rexprt_evaluated"] <= data["total_all"]


def test_rexprt_pathogenic_only_filter_excludes_unevaluated_and_benign_rows():
    with TestClient(app) as client:
        resp = client.get(
            "/api/catalog",
            params={"mode": "genome_wide", "rexprt_pathogenic_only": "true", "page_size": 200},
        )
    data = resp.json()
    assert data["total_matching"] > 0
    assert data["total_matching"] < data["total_all"]
    assert all(r["pathogenicity"] and r["pathogenicity"]["predicted_pathogenic"] for r in data["rows"])


def test_unevaluated_catalog_region_has_null_pathogenicity_not_fabricated():
    """chr1:20,797-20,863 has a real 18bp motif -- well outside RExPRT's 3-8bp
    trained range, so RExPRT's reference set never scored it. Its API row must
    show pathogenicity=None, not a made-up value."""
    with TestClient(app) as client:
        resp = client.get(
            "/api/catalog",
            params={"mode": "genome_wide", "chrom": "chr1", "page_size": 500},
        )
    data = resp.json()
    row = next(r for r in data["rows"] if r["start"] == 20797 and r["stop"] == 20863)
    assert row["pathogenicity"] is None
