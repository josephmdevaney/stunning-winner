from fastapi.testclient import TestClient

from repeatscope.api.main import app


def test_gene_finding_gba_bc2108_named_variants_not_trusted():
    """Real bc2108 data: havanno flags 14 named variants on GBA, but a same-size pseudogene means
    that pattern-match can't be trusted (see classify_gene), so no variant list is exposed."""
    with TestClient(app) as client:
        resp = client.get(
            "/api/gene-finding",
            params={"gene": "GBA", "mode": "carrier", "sample": "m84036_250809_035451_s4.bc2108"},
        )
    assert resp.status_code == 200
    data = resp.json()
    assert data["status"] != "Pathogenic"
    assert data["variants"] == []


def test_gene_finding_cyp21_bc2108_named_variants_not_trusted():
    with TestClient(app) as client:
        resp = client.get(
            "/api/gene-finding",
            params={"gene": "CYP21A2", "mode": "carrier", "sample": "m84036_250809_035451_s4.bc2108"},
        )
    assert resp.status_code == 200
    data = resp.json()
    assert data["status"] != "Pathogenic"
    assert data["variants"] == []


def test_gene_finding_hba_bc2108_surfaces_genotype_and_deletion_without_variants():
    """Real bc2108 data: no named pathogenic variants for HBA, but Paraphase's own
    genotype ("-a/aa") and sv_called (a real 3.7kb deletion) are surfaced verbatim
    in raw_fields -- honestly, without upgrading to a pathogenic classification."""
    with TestClient(app) as client:
        resp = client.get(
            "/api/gene-finding",
            params={"gene": "HBA1/2", "mode": "carrier", "sample": "m84036_250809_035451_s4.bc2108"},
        )
    assert resp.status_code == 200
    data = resp.json()
    assert data["status"] == "Uncertain significance"
    assert data["variants"] == []
    assert data["raw_fields"]["genotype"] == "-a/aa"
    assert "3p7del" in data["raw_fields"]["sv_called"]


def test_gene_finding_unknown_gene_404s():
    with TestClient(app) as client:
        resp = client.get(
            "/api/gene-finding",
            params={"gene": "NOT_A_REAL_GENE", "mode": "carrier", "sample": "m84036_250809_035451_s4.bc2108"},
        )
    assert resp.status_code == 404


def test_gene_finding_caller_is_the_genes_own_tool_not_trgt():
    """Non-repeat genes never pass through TRGT, so the finding must not carry a TRGT label."""
    with TestClient(app) as client:
        resp = client.get(
            "/api/gene-finding",
            params={"gene": "GBA", "mode": "carrier", "sample": "m84036_250809_035451_s4.bc2108"},
        )
    assert resp.status_code == 200
    data = resp.json()
    assert data["caller"] == data["tool"]
    assert "TRGT" not in data["caller"]
