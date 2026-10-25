"""Paraphase JSON region table (paraphase_annotate/regions.py and its routes).
Synthetic JSON in both Paraphase schema generations; no real sample data."""

import json

import pytest
from fastapi.testclient import TestClient

from paraphase_annotate import calls, regions, samples

SAMPLE = "S1"


def _region(**kw):
    base = {
        "total_cn": 4, "genes_in_region": "GENE1,GENE2", "highest_total_cn": 4, "two_copy_haplotypes": [],
        "genome_depth": 20, "region_depth": {"median": 40.0, "percentile80": 44.0}, "sample_sex": "male",
        "phase_region": "38:chr1:100-200",
        "final_haplotypes": {"1121": "X_gene1hap1", "2212": "X_gene2hap1"},
        "haplotype_details": {
            "X_gene1hap1": {"variants": ["110_A_G", "120_C_T"], "boundary": [100, 200], "boundary_gene2": None, "is_truncated": []},
            "X_gene2hap1": {"variants": [], "boundary": [100, 150], "boundary_gene2": [100, 150], "is_truncated": ["3"], "haplogroup": "S1-2"},
        },
        "unique_supporting_reads": {"X_gene1hap1": ["r1", "r2"], "X_gene2hap1": ["r3"]},
    }
    base.update(kw)
    return base


V3 = {
    "plain": _region(gene_cn=2, fusions_called=None, alleles_final=None),
    "smn1": _region(smn1_cn=None, smn2_cn=None, smn_del78_cn=0, total_cn=None, fusions_called=None),
    "smn1_ok": _region(smn1_cn=2, smn2_cn=1, smn_del78_cn=1, total_cn=None),
    "hba": _region(genotype="aa/aa", sv_called={"3p7": "deletion"}, fusions_called=None),
    "rccx": _region(phasing_success=False, annotated_alleles=[["h1"], ["h2"]], alleles_final=[["a", "b"]]),
    "CYP2D6": _region(fusions_called={"CYP2D6_fusionhap1": {"type": "duplication", "sequence": "2211", "breakpoint": [1, 2]}}),
    "empty": {"total_cn": None, "region_depth": None, "final_haplotypes": {}, "genes_in_region": "NOPE"},
}


@pytest.fixture
def env(tmp_path, monkeypatch):
    path = tmp_path / f"{SAMPLE}.paraphase.json"
    path.write_text(json.dumps(V3))
    monkeypatch.setattr(samples, "SAMPLES", {SAMPLE: str(path)})
    regions._load.cache_clear()
    return path


def _row(table, name):
    return next(r for r in table["regions"] if r["region"] == name)


def test_table_summarizes_every_region(env):
    t = regions.table(SAMPLE)
    assert t["n_regions"] == len(V3) and t["source_schema"] == "v3" and t["sample_sex"] == "male"
    plain = _row(t, "plain")
    assert plain["genes"] == ["GENE1", "GENE2"] and plain["n_haplotypes"] == 2
    assert plain["depth_median"] == 40.0 and plain["depth_vs_genome"] == 2.0
    assert plain["headline"] == "2 gene copies of 4 total" and plain["flags"] == []


def test_headlines_restate_only_json_fields(env):
    t = regions.table(SAMPLE)
    assert _row(t, "smn1")["headline"] == "SMN1 ? · SMN2 ?"
    assert _row(t, "smn1_ok")["headline"] == "SMN1 2 · SMN2 1 · SMNΔ7-8 1"
    assert _row(t, "hba")["headline"] == "genotype aa/aa"
    assert _row(t, "rccx")["headline"] == "phasing failed"


def test_flags(env):
    t = regions.table(SAMPLE)
    kinds = lambda name: [f["kind"] for f in _row(t, name)["flags"]]
    assert kinds("smn1") == ["cn_unresolved"]
    assert kinds("smn1_ok") == []
    assert kinds("hba") == ["sv"]
    assert kinds("rccx") == ["phasing_failed"]
    assert kinds("CYP2D6") == ["fusion"]
    assert "duplication" in _row(t, "CYP2D6")["flags"][0]["text"]
    assert kinds("empty") == ["no_data"]
    assert kinds("plain") == []


def test_medical_flag_matches_the_annotation_regions_case_insensitively(env):
    t = regions.table(SAMPLE)
    assert _row(t, "smn1")["medical"] is True and _row(t, "hba")["medical"] is True
    assert _row(t, "plain")["medical"] is False


def test_v4_json_gives_the_same_table_as_v3(env, tmp_path, monkeypatch):
    v3_table = regions.table(SAMPLE)
    v4_path = tmp_path / "S2.paraphase.json"
    v4_path.write_text(json.dumps(calls.to_v4(V3)))
    monkeypatch.setattr(samples, "SAMPLES", {"S2": str(v4_path)})
    regions._load.cache_clear()
    v4_table = regions.table("S2")
    assert v4_table["source_schema"] == "v4"
    assert v4_table["regions"] == v3_table["regions"]


def test_detail_lists_haplotypes_with_groups_boundaries_and_variants(env):
    d = regions.detail(SAMPLE, "plain")
    g1, g2 = d["haplotypes"]
    assert (g1["name"], g1["group"], g1["kind"]) == ("X_gene1hap1", "gene1", "gene")
    assert g1["n_variants"] == 2 and g1["variants"] == ["110_A_G", "120_C_T"] and g1["unique_reads"] == 2
    assert (g2["kind"], g2["truncated"], g2["haplogroup"]) == ("paralog / pseudogene", ["3"], "S1-2")


def test_detail_reports_fusions_alleles_and_scalar_fields(env):
    d = regions.detail(SAMPLE, "CYP2D6")
    assert d["fusions"] == [{"name": "CYP2D6_fusionhap1", "type": "duplication", "breakpoint": [1, 2]}]
    r = regions.detail(SAMPLE, "rccx")
    assert r["alleles"] == [["a", "b"]] and r["region_fields"]["phasing_success"] is False
    assert "annotated_alleles" not in r["region_fields"]  # nested structures stay in the JSON
    assert regions.detail(SAMPLE, "nope") is None


def test_variants_are_capped_but_counted(env, tmp_path):
    big = _region()
    big["haplotype_details"]["X_gene1hap1"]["variants"] = [f"{i}_A_G" for i in range(400)]
    (tmp_path / f"{SAMPLE}.paraphase.json").write_text(json.dumps({"big": big}))
    regions._load.cache_clear()
    h = regions.detail(SAMPLE, "big")["haplotypes"][0]
    assert h["n_variants"] == 400 and len(h["variants"]) == 25


def test_a_changed_file_is_reread(env, tmp_path):
    assert regions.table(SAMPLE)["n_regions"] == len(V3)
    (tmp_path / f"{SAMPLE}.paraphase.json").write_text(json.dumps({"only": _region()}))
    assert regions.table(SAMPLE)["n_regions"] == 1


def test_routes(env):
    from repeatscope.api.main import app
    with TestClient(app) as client:
        t = client.get("/api/paraphase/json", params={"sample": SAMPLE})
        assert t.status_code == 200 and t.json()["n_regions"] == len(V3)
        d = client.get("/api/paraphase/json/region", params={"sample": SAMPLE, "region": "hba"})
        assert d.status_code == 200 and d.json()["summary"]["headline"] == "genotype aa/aa"
        assert client.get("/api/paraphase/json/region", params={"sample": SAMPLE, "region": "zzz"}).status_code == 404
        assert client.get("/api/paraphase/json", params={"sample": "unknown"}).status_code == 404
