from fastapi.testclient import TestClient

from repeatscope.api import main
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


def test_meta_caller_is_each_samples_own_trgt_version():
    # The real GV panel VCFs were genotyped with TRGT 3.0.0; the HG002 run used 5.1.0.
    with TestClient(app) as client:
        panel = client.get("/api/meta", params={"mode": "puretarget_v2", "sample": "GV-01"}).json()
        genome = client.get("/api/meta", params={"mode": "genome_wide"}).json()
    assert panel["caller"] == "TRGT v3.0.0"
    assert genome["caller"] == "TRGT v5.1.0"


# --- Pangenome Repeat Catalog (second catalog) --------------------------------------------

def test_catalogs_lists_the_repeat_catalog_first_then_the_pangenome_one():
    with TestClient(app) as client:
        data = client.get("/api/catalogs").json()
    assert data["default"] == "repeat_catalog"
    assert [(c["id"], c["label"], c["total"]) for c in data["catalogs"]] == [
        ("repeat_catalog", "Repeat Catalog", 24_956),
        ("pangenome", "Pangenome Repeat Catalog", 5_709),
    ]


def test_catalog_defaults_to_the_repeat_catalog_and_names_it():
    with TestClient(app) as client:
        data = client.get("/api/catalog", params={"mode": "genome_wide", "page_size": 1}).json()
    assert data["catalog"]["id"] == "repeat_catalog"
    assert data["total_all"] == 24_956
    assert data["genotypes_available"] is True


def test_unknown_catalog_is_a_404():
    with TestClient(app) as client:
        assert client.get("/api/catalog", params={"catalog": "nope"}).status_code == 404


def test_pangenome_catalog_serves_its_own_loci_with_real_hg002_genotypes():
    with TestClient(app) as client:
        data = client.get(
            "/api/catalog", params={"mode": "genome_wide", "catalog": "pangenome", "page_size": 6000}
        ).json()
    assert data["catalog"]["label"] == "Pangenome Repeat Catalog"
    assert data["total_all"] == data["total_matching"] == len(data["rows"]) == 5_709
    assert data["genotypes_available"] is True
    # TRGT called every locus but one (chr20:43,189,969-43,189,980, a no-call in HG002).
    assert sum(r["genotyped"] for r in data["rows"]) == 5_708

    # C9orf72's GGCCCC repeat is in this catalog although RExPRT's benign score for it kept it out
    # of the Repeat Catalog (see test_c9orf72_is_absent_...).
    c9 = next(r for r in data["rows"] if r["chrom"] == "chr9" and r["start"] == 27_573_484)
    assert c9["stop"] == 27_573_546 and c9["motifs"] == "GGCCCC"
    assert c9["genes"] == ["C9orf72"]
    assert c9["genotyped"] is True and c9["depth"] == 44


def test_pangenome_catalog_filters_like_the_repeat_catalog():
    with TestClient(app) as client:
        chr9 = client.get(
            "/api/catalog", params={"mode": "genome_wide", "catalog": "pangenome", "chrom": "chr9", "page_size": 6000}
        ).json()
        ggcccc = client.get(
            "/api/catalog", params={"mode": "genome_wide", "catalog": "pangenome", "motif": "ggcccc", "page_size": 6000}
        ).json()
    assert 0 < chr9["total_matching"] < chr9["total_all"] and {r["chrom"] for r in chr9["rows"]} == {"chr9"}
    assert any(r["start"] == 27_573_484 for r in ggcccc["rows"])
    assert all("GGCCCC" in r["motifs"].upper() for r in ggcccc["rows"])


def test_pangenome_genotypes_come_only_from_a_call_on_that_very_locus(monkeypatch):
    """Without a pangenome VCF the mode's own genome-wide VCF is used, which here holds only the
    Repeat Catalog's loci. A pangenome locus it doesn't hold must come back ungenotyped, not
    take the call of a different locus that happens to overlap it."""
    with TestClient(app) as client:
        mode = main._modes["genome_wide"]
        monkeypatch.setitem(mode, "extra_catalog_sample_paths", {})
        data = client.get(
            "/api/catalog", params={"mode": "genome_wide", "catalog": "pangenome", "page_size": 6000}
        ).json()
    assert data["genotypes_available"] is True
    assert sum(r["genotyped"] for r in data["rows"]) == 1_098


def test_pangenome_catalog_without_any_genome_wide_vcf_reports_no_genotypes(monkeypatch):
    with TestClient(app) as client:
        mode = main._modes["genome_wide"]
        monkeypatch.setitem(mode, "extra_catalog_sample_paths", {})
        monkeypatch.setitem(mode, "catalog_sample_paths", {})
        data = client.get(
            "/api/catalog", params={"mode": "genome_wide", "catalog": "pangenome", "page_size": 6000}
        ).json()
    assert data["genotypes_available"] is False
    assert data["total_all"] == 5_709
    assert not any(r["genotyped"] for r in data["rows"])


def test_catalog_rexprt_counts_cover_just_that_catalogs_loci():
    with TestClient(app) as client:
        pangenome = client.get(
            "/api/catalog", params={"mode": "genome_wide", "catalog": "pangenome", "page_size": 6000}
        ).json()
        default = client.get("/api/catalog", params={"mode": "genome_wide", "page_size": 1}).json()
    scored = [r["pathogenicity"] for r in pangenome["rows"] if r["pathogenicity"]]
    assert pangenome["n_catalog_rexprt_evaluated"] == len(scored)
    assert pangenome["n_catalog_rexprt_predicted_pathogenic"] == sum(p["predicted_pathogenic"] for p in scored)
    assert 0 < pangenome["n_catalog_rexprt_predicted_pathogenic"] < pangenome["n_catalog_rexprt_evaluated"] < 5_709
    # The Repeat Catalog holds only RExPRT-positive loci, so all of its evaluated loci are positive.
    assert default["n_catalog_rexprt_evaluated"] == default["n_catalog_rexprt_predicted_pathogenic"] == default["total_all"]
    # The whole-table counts are the same whichever catalog is shown.
    assert default["n_rexprt_evaluated"] == pangenome["n_rexprt_evaluated"] == 727_004
