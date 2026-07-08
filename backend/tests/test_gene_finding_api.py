from fastapi.testclient import TestClient

from repeatscope.api.main import app


def test_gene_finding_gba_bc2108_lists_all_named_variants():
    """Real bc2108 data: GBA has 14 named pathogenic variants across 2 of 4
    haplotype-phase calls. The detail endpoint should expose the full structured
    list (not just embedded in the summary note text)."""
    with TestClient(app) as client:
        resp = client.get(
            "/api/gene-finding",
            params={"gene": "GBA", "mode": "carrier", "sample": "m84036_250809_035451_s4.bc2108"},
        )
    assert resp.status_code == 200
    data = resp.json()
    assert data["status"] == "Pathogenic"
    assert data["zygosity"] == "Heterozygous"
    assert len(data["variants"]) == 14
    assert "GBA1_c.115+1G>A" in data["variants"]
    assert data["raw_fields"]["haplotype_phase_calls"] == "2 of 4"


def test_gene_finding_cyp21_bc2108_lists_all_named_variants():
    with TestClient(app) as client:
        resp = client.get(
            "/api/gene-finding",
            params={"gene": "CYP21A2", "mode": "carrier", "sample": "m84036_250809_035451_s4.bc2108"},
        )
    assert resp.status_code == 200
    data = resp.json()
    assert data["status"] == "Pathogenic"
    assert len(data["variants"]) == 9
    assert "CYP21A2_p.Gln319Ter_p.Gln319Ter_c.955C>T" in data["variants"]


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
