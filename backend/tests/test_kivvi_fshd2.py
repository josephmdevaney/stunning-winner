"""Tests for the FSHD2 research-evidence tiers (kivvi/fshd2.py) and how the sample registry feeds them.

Inputs are small synthetic Kivvi-shaped files and made-up variants -- no real per-participant data is read or committed.
"""

import json

import pytest
from fastapi.testclient import TestClient

from kivvi import fshd2, parse, samples
from repeatscope.api.main import app

client = TestClient(app)


def _allele(size, chromosome="chr4:Group1.1", distal="qAIntactPolyA", allele_type="assembled", methylation=0.5):
    n = int(str(size).lstrip(">="))
    return {
        "allele_name": "LeftFlank-" + "-".join(str(i) for i in range(1, n + 1)) + "-RightFlank",
        "chromosome": chromosome,
        "distal_haplotype": distal,
        "allele_type": allele_type,
        "allele_size": str(size),
        "methylation": methylation,
        "ending_in_qal": False,
    }


PERMISSIVE = lambda: _allele(14)  # noqa: E731
NOT_PERMISSIVE = [
    lambda: _allele(20, chromosome="chr10:Group2.1", distal="qADisruptedPolyA"),
    lambda: _allele(25, chromosome="chr4:Group1.2", distal="qB"),
    lambda: _allele(22, chromosome="chr10:Group2.2", distal="qB"),
]


def _fshd2_block(alleles, sample_methylation):
    data = {"allele_info": alleles, "additional": {"methylation_all_sites": sample_methylation, "qc_metrics": {}}}
    return parse.shape_d4z4(data)["fshd2"]


def _full_sample(sample_methylation, permissive=True):
    """Four resolved alleles, with or without a permissive one."""
    others = [make() for make in NOT_PERMISSIVE]
    return _fshd2_block(([PERMISSIVE()] if permissive else [_allele(30, chromosome="chr4:Group2.1")]) + others, sample_methylation)


LOF = {"gene": "SMCHD1", "hgvs_c": "c.181dup", "hgvs_p": "p.(C61Lfs*9)"}
VUS_MISSENSE = {"gene": "SMCHD1", "hgvs_c": "c.1844T>G", "hgvs_p": "p.(V615G)", "clinvar_classification": "Uncertain significance"}


# --- consequence ------------------------------------------------------------

@pytest.mark.parametrize(
    "hgvs_c, hgvs_p, expected",
    [
        ("c.3515-2A>G", None, "splice_site"),
        ("c.3515-2 A>G", "", "splice_site"),
        ("c.100+1G>A", None, "splice_site"),
        ("c.100+6T>C", None, "other"),                    # not a canonical splice position
        ("c.1844T>G", "p.(V615G)", "missense"),
        ("c.1844T>G", "p.Val615Gly", "missense"),
        ("c.181dup", "p.(C61Lfs*9)", "frameshift"),
        ("c.180_183delGTGT", "p.C61RfsX47", "frameshift"),
        ("c.2412del", "p.(Y804*)", "nonsense"),
        ("c.1030C>T", "p.R344X", "nonsense"),
        ("c.4267C>T", "p.(R1423*)", "nonsense"),
        ("c.4267C>T", "p.Arg1423Ter", "nonsense"),
        ("c.300del", None, "frameshift"),                # no protein change: a lone deletion shifts the frame
        ("c.300_302del", None, "other"),                 # three bases: in frame
        ("c.300_303del", None, "frameshift"),
        (None, None, "other"),
    ],
)
def test_consequence_from_hgvs(hgvs_c, hgvs_p, expected):
    assert fshd2.consequence(hgvs_c, hgvs_p) == expected


def test_a_splice_site_outranks_the_protein_change():
    assert fshd2.consequence("c.3515-2A>G", "p.?") == "splice_site"


# --- which variants count ---------------------------------------------------

def test_loss_of_function_counts_without_clinvar():
    assert fshd2.annotate_variant(LOF)["counts"] is True


def test_missense_of_uncertain_significance_is_shown_but_does_not_count():
    v = fshd2.annotate_variant(VUS_MISSENSE)
    assert v["consequence"] == "missense" and v["counts"] is False


def test_missense_counts_when_clinvar_says_pathogenic():
    v = fshd2.annotate_variant({**VUS_MISSENSE, "clinvar_classification": "Likely pathogenic"})
    assert v["counts"] is True and "ClinVar" in v["counts_reason"]


def test_a_variant_outside_the_fshd2_genes_does_not_count():
    v = fshd2.annotate_variant({**LOF, "gene": "TTN"})
    assert v["in_fshd2_gene"] is False and v["counts"] is False


@pytest.mark.parametrize("gene", ["SMCHD1", "smchd1", "DNMT3B", "LRIF1"])
def test_the_three_fshd2_genes_are_recognised_in_any_case(gene):
    assert fshd2.annotate_variant({**LOF, "gene": gene})["in_fshd2_gene"] is True


# --- tiers ------------------------------------------------------------------

def test_tier_1_needs_all_three_findings():
    e = fshd2.evidence(_full_sample(0.25), [LOF])
    assert (e["status"], e["tier"]) == ("all_three", 1)
    assert [e["criteria"][k]["met"] for k in ("permissive_allele", "hypomethylation", "variant")] == [True, True, True]
    assert e["pending_variant"] is False


def test_tier_2_is_the_permissive_allele_with_one_other_finding():
    low_methylation_only = fshd2.evidence(_full_sample(0.25), [])
    variant_only = fshd2.evidence(_full_sample(0.60), [LOF])
    assert (low_methylation_only["status"], low_methylation_only["tier"]) == ("two_of_three", 2)
    assert (variant_only["status"], variant_only["tier"]) == ("two_of_three", 2)


def test_tier_3_is_the_permissive_allele_alone():
    e = fshd2.evidence(_full_sample(0.60), [])
    assert (e["status"], e["tier"]) == ("permissive_only", 3)


def test_an_uncertain_variant_does_not_lift_the_tier():
    e = fshd2.evidence(_full_sample(0.25), [VUS_MISSENSE])
    assert (e["status"], e["tier"]) == ("two_of_three", 2)
    assert e["criteria"]["variant"]["state"] == "uncertain"
    assert any("does not count" in n for n in e["notes"])


def test_no_permissive_allele_means_no_tier_whatever_else_is_present():
    e = fshd2.evidence(_full_sample(0.10, permissive=False), [LOF])
    assert (e["status"], e["tier"]) == ("no_permissive", None)
    assert e["label"].startswith("Not supported")
    assert any("without a permissive allele" in n for n in e["notes"])


def test_the_4a166_haplotype_is_not_permissive():
    block = _fshd2_block([_allele(14, chromosome="chr4:Group2.1")] + [make() for make in NOT_PERMISSIVE], 0.2)
    assert fshd2.evidence(block, [LOF])["status"] == "no_permissive"


def test_unassessed_variant_is_not_the_same_as_none_found():
    block = _full_sample(0.25)
    not_assessed = fshd2.evidence(block, None)
    none_found = fshd2.evidence(block, [])
    assert not_assessed["criteria"]["variant"]["state"] == "not_assessed" and not_assessed["pending_variant"] is True
    assert none_found["criteria"]["variant"]["state"] == "none_reported" and none_found["pending_variant"] is False
    assert not_assessed["status"] == none_found["status"] == "two_of_three"


def test_no_pending_variant_when_there_is_no_permissive_allele():
    assert fshd2.evidence(_full_sample(0.25, permissive=False), None)["pending_variant"] is False


def test_missing_sample_methylation_is_noted_and_does_not_count():
    e = fshd2.evidence(_full_sample(None), [LOF])
    assert e["criteria"]["hypomethylation"]["met"] is False
    assert (e["status"], e["tier"]) == ("two_of_three", 2)
    assert any("could not be scored" in n for n in e["notes"])


def test_cutoff_boundary_is_strictly_below_30_percent():
    assert fshd2.evidence(_full_sample(0.299), [LOF])["status"] == "all_three"
    assert fshd2.evidence(_full_sample(0.30), [LOF])["status"] == "two_of_three"


def test_fewer_than_four_alleles_is_flagged_when_there_is_no_permissive_allele():
    block = _fshd2_block([_allele(20, chromosome="chr10:Group2.1", distal="qADisruptedPolyA"), _allele(25, chromosome="chr4:Group1.2", distal="qB")], 0.2)
    e = fshd2.evidence(block, [LOF])
    assert any("2 of 4 expected" in n for n in e["notes"])


def test_an_unplaced_qa_intact_allele_is_flagged_when_there_is_no_permissive_allele():
    block = _fshd2_block([_allele(">8", chromosome="unknown", allele_type="partial")] + [make() for make in NOT_PERMISSIVE], 0.2)
    e = fshd2.evidence(block, [])
    assert e["status"] == "no_permissive"
    assert any("could not be placed" in n for n in e["notes"])


def test_a_partial_permissive_allele_is_noted():
    block = _fshd2_block([_allele(">12", allele_type="merged")] + [make() for make in NOT_PERMISSIVE], 0.2)
    e = fshd2.evidence(block, [LOF])
    assert e["status"] == "all_three"
    assert any("lower bound" in n for n in e["notes"])


# --- the registry and API ---------------------------------------------------

@pytest.fixture
def root_with_variants(tmp_path, monkeypatch):
    for name, alleles in {"S1": [PERMISSIVE()] + [make() for make in NOT_PERMISSIVE], "S2": [PERMISSIVE()] + [make() for make in NOT_PERMISSIVE]}.items():
        d = tmp_path / name / "d4z4"
        d.mkdir(parents=True)
        data = {"allele_info": alleles, "additional": {"methylation_all_sites": 0.2, "qc_metrics": {}}}
        (d / f"{name}.kivvi.d4z4.json").write_text(json.dumps(data))
    # S1 has a variant, S2 was checked and has none; a third sample would be absent from the file altogether.
    (tmp_path / "fshd2_variants.json").write_text(json.dumps({"samples": {"S1": {"variants": [LOF]}, "S2": {"variants": []}}}))
    monkeypatch.setattr(samples, "RESULTS_ROOT", tmp_path)
    monkeypatch.delenv("KIVVI_FSHD2_VARIANTS", raising=False)
    return tmp_path


def test_the_variants_file_is_not_mistaken_for_a_sample(root_with_variants):
    assert samples.sample_keys() == ["S1", "S2"]


def test_variants_are_looked_up_per_sample(root_with_variants):
    assert [v["hgvs_c"] for v in samples.fshd2_variants("S1")] == ["c.181dup"]
    assert samples.fshd2_variants("S2") == []
    assert samples.fshd2_variants("S3") is None


def test_target_detail_carries_the_tier(root_with_variants):
    s1 = client.get("/api/kivvi/target/d4z4", params={"sample": "S1"}).json()["fshd2"]["evidence"]
    s2 = client.get("/api/kivvi/target/d4z4", params={"sample": "S2"}).json()["fshd2"]["evidence"]
    assert (s1["status"], s1["tier"]) == ("all_three", 1)
    assert (s2["status"], s2["tier"]) == ("two_of_three", 2)
    assert s1["variants"][0]["consequence"] == "frameshift" and s1["variants"][0]["counts"] is True


def test_overview_carries_the_tier_summary(root_with_variants):
    fshd2_overview = client.get("/api/kivvi/overview", params={"sample": "S1"}).json()["d4z4"]["fshd2"]
    assert fshd2_overview == {"status": "all_three", "tier": 1, "label": fshd2.TIER_LABELS["all_three"], "pending_variant": False}


def test_without_a_variants_file_every_sample_is_unassessed(root_with_variants):
    (root_with_variants / "fshd2_variants.json").unlink()
    e = client.get("/api/kivvi/target/d4z4", params={"sample": "S1"}).json()["fshd2"]["evidence"]
    assert e["criteria"]["variant"]["state"] == "not_assessed" and e["pending_variant"] is True


def test_the_variants_file_can_be_pointed_elsewhere(root_with_variants, tmp_path_factory, monkeypatch):
    elsewhere = tmp_path_factory.mktemp("elsewhere") / "calls.json"
    elsewhere.write_text(json.dumps({"samples": {"S1": {"variants": []}}}))
    monkeypatch.setenv("KIVVI_FSHD2_VARIANTS", str(elsewhere))
    assert samples.fshd2_variants("S1") == []
    assert samples.fshd2_variants("S2") is None
