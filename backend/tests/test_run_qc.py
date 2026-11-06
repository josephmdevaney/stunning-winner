from fastapi.testclient import TestClient

from repeatscope.api.main import app


def test_run_qc_scoped_to_carrier_panel_not_genome_wide_catalog():
    """Carrier is a targeted 5-locus panel; the genome-wide 171k-region repeat
    catalog uses a different ID space and this VCF has no calls outside its own
    panel, so Run QC should report against the same panel Overview/Sample QC use
    -- not silently show "0 of 171,146 genotyped", which reads as broken."""
    with TestClient(app) as client:
        resp = client.get(
            "/api/qc/run", params={"mode": "carrier", "sample": "m84036_250809_035451_s4.bc2192"}
        )
    assert resp.status_code == 200
    data = resp.json()
    assert data["scoped_to_panel"] is True
    assert data["repeat_regions_total"] == 5
    assert data["repeat_regions_genotyped"] == 5
    assert "assay panel" in data["note"]


def test_run_qc_genome_wide_mode_still_uses_full_catalog():
    with TestClient(app) as client:
        resp = client.get("/api/qc/run", params={"mode": "genome_wide"})
    assert resp.status_code == 200
    data = resp.json()
    assert data["scoped_to_panel"] is False
    assert data["repeat_regions_total"] > 1000



# --- Sample QC depth floor: 20x for PureTarget data (PacBio's minimum target coverage), 10x genome-wide ---

def test_min_locus_depth_by_mode():
    from repeatscope.api.main import _min_locus_depth

    assert _min_locus_depth("puretarget_v2") == 20
    assert _min_locus_depth("carrier") == 20
    assert _min_locus_depth("genome_wide") == 10  # WGS coverage: PacBio's enrichment spec doesn't apply


def test_sample_qc_rows_carry_the_modes_floor():
    with TestClient(app) as client:
        data = client.get("/api/qc/sample", params={"mode": "puretarget_v2", "sample": "GV-01"}).json()
    assert data["rows"] and {r["min_required"] for r in data["rows"]} == {20}


def test_puretarget_loci_between_10x_and_20x_are_flagged(monkeypatch):
    import types

    from repeatscope.api import main

    fake = [types.SimpleNamespace(genotyped=True, depth=d, gene=f"G{d}", id=f"ID{d}") for d in (9, 15, 20, 25)]
    fake.append(types.SimpleNamespace(genotyped=False, depth=None, gene="NC", id="NC"))
    monkeypatch.setattr(main, "_get_panel", lambda m, s: fake)
    with TestClient(app) as client:
        data = client.get("/api/qc/sample", params={"mode": "puretarget_v2", "sample": "GV-01"}).json()
    assert {r["gene"]: r["flag"] for r in data["rows"]} == {
        "G9": "Low depth", "G15": "Low depth", "G20": "Pass", "G25": "Pass", "NC": "No call",
    }
    assert data["overall_status"] == "Review"


def test_genome_wide_loci_keep_the_10x_floor(monkeypatch):
    import types

    from repeatscope.api import main

    fake = [types.SimpleNamespace(genotyped=True, depth=d, gene=f"G{d}", id=f"ID{d}") for d in (9, 15)]
    monkeypatch.setattr(main, "_get_panel", lambda m, s: fake)
    with TestClient(app) as client:
        data = client.get("/api/qc/sample", params={"mode": "genome_wide"}).json()
    assert {r["gene"]: r["flag"] for r in data["rows"]} == {"G9": "Low depth", "G15": "Pass"}
    assert {r["min_required"] for r in data["rows"]} == {10}
