"""Epigenetics routes, through the real app, on a synthetic registered sample (cpg reads are patched)."""

import json

import pytest
from fastapi.testclient import TestClient

from epigenetics import cpg, service
from epigenetics.cpg import CpG
from repeatscope.api.main import app


def _cpgs(n, pct, cov=30, start=1000):
    mod = round(cov * pct / 100)
    return [CpG(start + 10 * i, float(pct), cov, mod, cov - mod) for i in range(n)]


@pytest.fixture
def client(tmp_path, monkeypatch):
    for t in ("combined", "hap1", "hap2"):
        (tmp_path / f"SYN.{t}.bed.gz").write_bytes(b"")
        (tmp_path / f"SYN.{t}.bed.gz.tbi").write_bytes(b"")
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"samples": {"SYN": {"dir": str(tmp_path), "prefix": "SYN", "label": "Synthetic genome"}}}))
    monkeypatch.setenv("EPIGENETICS_CONFIG", str(cfg))

    # SNURF window -> both alleles methylated (a PWS-like pattern); everything else a normal imprint
    snurf_core = service.panel()["regions"]
    snurf_start = next(r for r in snurf_core if r["id"] == "snurf-tss-dmr")["core"]["start"]

    def fake_fetch(bed, chrom, start, end):
        track = bed.name.split(".")[1]
        if chrom == "chr15" and start == snurf_start:
            pct = {"combined": 96, "hap1": 96, "hap2": 97}[track]
        else:
            pct = {"combined": 48, "hap1": 95, "hap2": 3}[track]
        return _cpgs(25, pct, 30, start)

    monkeypatch.setattr(cpg, "fetch", fake_fetch)
    monkeypatch.setattr(cpg, "header", lambda bed: {"pb_cpg_tools_version": "3.0.0", "pileup_mode": "model", "modsites_mode": "denovo", "min_coverage": 4})
    service._cache.clear()
    service._cpg_cache.clear()
    with TestClient(app) as c:
        yield c


def test_samples(client):
    body = client.get("/api/epigenetics/samples").json()
    assert body["default"] == "SYN" and body["samples"][0]["tracks"] == ["combined", "hap1", "hap2"]


def test_imprinting_overview_reports_the_pattern_and_what_it_was_computed_from(client):
    o = client.get("/api/epigenetics/imprinting", params={"sample": "SYN"}).json()
    assert o["summary"]["regions"] == 48 and o["provenance"]["pileup_mode"] == "model"
    d = {x["id"]: x for x in o["disorders"]}
    assert d["pws"]["status"] == "pattern_present" and d["as"]["status"] == "not_observed"
    assert o["panel"]["pmcid"] == "PMC11295890" and "not a diagnosis" in o["caveats"]


def test_thresholds_can_be_overridden_and_are_echoed(client):
    o = client.get("/api/epigenetics/imprinting", params={"sample": "SYN", "gain_above": 99}).json()
    assert o["thresholds"]["gain_above"] == 99.0
    d = {x["id"]: x for x in o["disorders"]}
    assert d["pws"]["status"] == "not_observed"            # 96% is no longer "gain" at a 99% cutoff


def test_bad_parameters_are_rejected(client):
    assert client.get("/api/epigenetics/imprinting", params={"window": "nope"}).status_code == 422
    assert client.get("/api/epigenetics/imprinting", params={"loss_below": 90, "gain_above": 10}).status_code == 422
    assert client.get("/api/epigenetics/imprinting", params={"sample": "missing"}).status_code == 404
    assert client.get("/api/epigenetics/imprinting/region", params={"id": "nope"}).status_code == 404
    assert client.get("/api/epigenetics/imprinting/region", params={"id": "snurf-tss-dmr", "flank": 10**6}).status_code == 422


def test_region_endpoint(client):
    r = client.get("/api/epigenetics/imprinting/region", params={"id": "snurf-tss-dmr", "flank": 500}).json()
    assert r["name"] == "SNURF:TSS-DMR" and len(r["tracks"]["combined"]) == 25


def test_no_registered_samples_gives_an_actionable_message(tmp_path, monkeypatch):
    monkeypatch.setenv("EPIGENETICS_CONFIG", str(tmp_path / "none.json"))
    with TestClient(app) as c:
        r = c.get("/api/epigenetics/imprinting")
        assert r.status_code == 404 and "epigenetics.local.json.example" in r.json()["detail"]
        assert c.get("/api/epigenetics/samples").json() == {"samples": [], "default": ""}
