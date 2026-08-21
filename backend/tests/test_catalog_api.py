from fastapi.testclient import TestClient

from repeatscope.api.main import app


def test_catalog_rows_include_real_rexprt_pathogenicity():
    # This catalog's variation-cluster consolidation labels FMR1's merged
    # region with just one constituent sub-repeat's motif ("TGG"), not "CGG",
    # even though the real repeat and RExPRT's own independent score are
    # CGG-based -- so filter by chrom only, not a motif substring, which
    # would miss it. chrX has 1006 rows total; page_size covers all of them
    # rather than relying on FMR1's own sort position.
    with TestClient(app) as client:
        resp = client.get(
            "/api/catalog",
            params={"mode": "genome_wide", "chrom": "chrX", "page_size": 2000},
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
    # n_rexprt_evaluated/n_rexprt_predicted_pathogenic report against the
    # *full* RExPRT reference table (727,004 real matches against the full
    # 5,339,786-region v2.0 catalog -- see project/uploads/README.md),
    # independent of how many regions the currently-deployed, rescoped
    # catalog (total_all) includes -- so the only universally-true relation
    # between them is pathogenic <= evaluated. total_all/total_matching are
    # a separate, catalog-scoped pair.
    with TestClient(app) as client:
        resp = client.get("/api/catalog", params={"mode": "genome_wide", "page_size": 1})
    data = resp.json()
    assert 0 < data["n_rexprt_predicted_pathogenic"] <= data["n_rexprt_evaluated"]
    assert 0 < data["total_matching"] <= data["total_all"]


def test_catalog_is_rescoped_to_rexprt_pathogenic_regions_only():
    """The deployed genome-wide catalog is scoped to just the regions RExPRT
    predicts pathogenic (see project/uploads/README.md) -- so every single
    row must show real, attached, positive pathogenicity. The
    rexprt_pathogenic_only filter is therefore a real no-op against this
    catalog specifically: it excludes nothing."""
    with TestClient(app) as client:
        unfiltered = client.get("/api/catalog", params={"mode": "genome_wide", "page_size": 500}).json()
        filtered = client.get(
            "/api/catalog",
            params={"mode": "genome_wide", "rexprt_pathogenic_only": "true", "page_size": 500},
        ).json()
    assert filtered["total_matching"] == unfiltered["total_matching"] == unfiltered["total_all"]
    assert all(r["pathogenicity"] and r["pathogenicity"]["predicted_pathogenic"] for r in unfiltered["rows"])


def test_c9orf72_is_absent_a_real_consequence_of_the_rescoping_not_a_bug():
    """C9orf72's real GGGGCC region (chr9:27,573,523-27,573,546) scores
    ensemble_binary=0 (benign) in RExPRT's own reference table -- see
    test_rexprt.py's test_real_c9orf72_locus_is_a_documented_false_negative
    -- so it is correctly absent from this RExPRT-positive-only catalog. This
    documents that real, visible cost directly, rather than leaving it as a
    silent gap."""
    with TestClient(app) as client:
        resp = client.get("/api/catalog", params={"mode": "genome_wide", "chrom": "chr9", "page_size": 5000})
    data = resp.json()
    assert not any(r["start"] == 27_573_523 and r["stop"] == 27_573_546 for r in data["rows"])
