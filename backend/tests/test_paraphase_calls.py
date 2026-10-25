"""Paraphase JSON schema normalization (paraphase_annotate/calls.py). Synthetic
data only -- real per-participant Paraphase output never goes into tests."""

import copy
import json

from paraphase_annotate import calls


def _v3_smn() -> dict:
    return {
        "smn1_cn": 2, "smn2_cn": 1, "smn_del78_cn": 0,
        "smn1_haplotypes": {"h1": "smn1hap1"},
        "highest_total_cn": 3, "final_haplotypes": {"h1": "1121"},
        "region_depth": {"median": 90.0, "percentile80": 95.0},
        "genome_depth": 30, "genes_in_region": "SMN1, SMN2",
        "gene_cn": 3, "fusions_called": None,
    }


def test_to_v4_moves_region_specific_fields_and_keeps_general_ones():
    out = calls.to_v4({"smn1": _v3_smn()})["smn1"]
    assert set(out) == {
        "highest_total_cn", "final_haplotypes", "region_depth", "genome_depth", "genes_in_region",
        "region_specific_info",
    }
    assert out["region_specific_info"] == {
        "smn1_cn": 2, "smn2_cn": 1, "smn_del78_cn": 0, "smn1_haplotypes": {"h1": "smn1hap1"},
        "gene_cn": 3, "fusions_called": None,
    }


def test_to_v4_gives_a_plain_region_an_empty_region_specific_info():
    plain = {"total_cn": 2, "final_haplotypes": {}, "region_depth": {"median": 40.0}}
    out = calls.to_v4({"AMY1A": plain})["AMY1A"]
    assert out["region_specific_info"] == {}
    assert out["total_cn"] == 2


def test_to_v4_treats_an_unlisted_field_as_region_specific():
    """A newer Paraphase may add a per-region field this module has never heard
    of; it must land in region_specific_info, where the v4 script looks for it."""
    region = {"total_cn": 2, "brand_new_field": [1, 2]}
    assert calls.to_v4({"x": region})["x"]["region_specific_info"] == {"brand_new_field": [1, 2]}


def test_to_v4_leaves_a_v4_region_alone_and_does_not_mutate_the_input():
    v4_region = {"total_cn": 2, "region_specific_info": {"gene_cn": 2}}
    v3 = {"smn1": _v3_smn()}
    snapshot = copy.deepcopy(v3)
    out = calls.to_v4({"a": v4_region, **v3})
    assert out["a"] is v4_region
    assert v3 == snapshot
    assert calls.to_v4(out) == out  # idempotent


def test_schema_detection():
    assert calls.schema({}) is None
    assert calls.schema({"smn1": _v3_smn()}) == "v3"
    assert calls.schema(calls.to_v4({"smn1": _v3_smn()})) == "v4"
    assert calls.schema({"a": {"total_cn": 2}, "b": {"region_specific_info": {}}}) == "v4"


def test_load_calls_rejects_a_non_object(tmp_path):
    p = tmp_path / "x.paraphase.json"
    p.write_text(json.dumps([1, 2]))
    try:
        calls.load_calls(p)
    except ValueError as e:
        assert "expected a JSON object" in str(e)
    else:
        raise AssertionError("expected ValueError")


def test_vcf_provenance_reads_version_and_command(tmp_path):
    d = tmp_path / "S1_paraphase_vcfs"
    d.mkdir()
    (d / "S1_smn1.vcf").write_text(
        "##fileformat=VCFv4.2\n"
        "##paraphase_version=4.1.0\n"
        "##paraphase_command=paraphase --bam S1.bam --out ./\n"
        "#CHROM\tPOS\n"
    )
    info = calls.vcf_provenance(d)
    assert info == {"paraphase_version": "4.1.0", "paraphase_command": "paraphase --bam S1.bam --out ./"}


def test_vcf_provenance_is_empty_when_the_directory_is_missing(tmp_path):
    assert calls.vcf_provenance(tmp_path / "nope") == {"paraphase_version": None, "paraphase_command": None}
