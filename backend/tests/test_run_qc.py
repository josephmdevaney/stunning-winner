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
