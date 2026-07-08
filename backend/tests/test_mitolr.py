"""Tests for the MitoLR data generation and API routes."""

from fastapi.testclient import TestClient

from mitolr import generate, qc
from repeatscope.api.main import app

client = TestClient(app)


def test_variants_are_deterministic():
    a = generate.build_variants()
    b = generate.build_variants()
    assert [v["id"] for v in a] == [v["id"] for v in b]
    assert [v["heteroplasmy"] for v in a] == [v["heteroplasmy"] for v in b]


def test_panel_and_novel_counts():
    variants = generate.build_variants()
    # 21 curated pathogenic + 58 synthetic novel calls.
    assert len(variants) == 79
    assert sum(1 for v in variants if v["origin"] == "pathogenic") == 21
    assert sum(1 for v in variants if v["origin"] == "novel") == 58


def test_variants_sorted_by_position():
    variants = generate.build_variants()
    positions = [v["pos"] for v in variants]
    assert positions == sorted(positions)


def test_narrative_overrides_applied():
    by_id = {v["id"]: v for v in generate.build_variants()}
    assert by_id["m3243AG"]["heteroplasmy"] == 38.0
    assert by_id["m3243AG"]["status"] == "Pathogenic"
    assert by_id["m1555AG"]["heteroplasmy"] == 97.0
    assert by_id["kss_common_del"]["heteroplasmy"] == 11.0
    assert by_id["kss_common_del"]["status"] == "Low-level heteroplasmy"


def test_pileup_is_subsampled_and_preserves_ratio():
    by_id = {v["id"]: v for v in generate.build_variants()}
    v = by_id["m3243AG"]
    reads = generate.simulate_pileup(v)
    assert len(reads) == generate.PILEUP_DISPLAY_READS
    n_mut = sum(1 for r in reads if r["kind"] == "mut")
    # ~38% of 48 reads carry the mutant allele.
    assert n_mut == round(generate.PILEUP_DISPLAY_READS * v["heteroplasmy"] / 100)


def test_qc_shapes():
    data = qc.build_qc()
    assert len(data["read_length_histogram"]) == 18
    assert len(data["coverage_track"]) == 60
    assert data["coverage_track"][0]["pos"] == 0
    assert data["coverage_track"][-1]["pos"] == generate.genome_length()
    assert len(data["regions"]) == 8
    assert data["unavailable_metrics"]


def test_overview_endpoint():
    r = client.get("/api/mito/overview")
    assert r.status_code == 200
    body = r.json()
    assert body["stats"]["variants_called"] == 79
    assert body["stats"]["pathogenic_findings"] == 2
    assert body["haplogroup"]["name"] == "H2a5b1"
    assert "all" not in body["consequence_options"]


def test_variant_endpoint_and_404():
    r = client.get("/api/mito/variant/m3243AG")
    assert r.status_code == 200
    body = r.json()
    assert body["variant"]["id"] == "m3243AG"
    assert body["alt_base"] == "G"
    assert len(body["gene_map"]) == 15
    assert client.get("/api/mito/variant/does-not-exist").status_code == 404


def test_deletion_variant_detail():
    r = client.get("/api/mito/variant/kss_common_del")
    assert r.status_code == 200
    body = r.json()
    assert body["variant"]["isDeletion"] is True
    assert body["variant"]["stop"] == 13459


def test_repeatscope_routes_still_present():
    # Mounting MitoLR must not disturb the existing product.
    assert client.get("/api/modes").status_code == 200
