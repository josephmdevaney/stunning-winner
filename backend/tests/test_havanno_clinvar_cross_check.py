"""A havanno-flagged "pathogenic" variant that ClinVar currently classifies as Benign is dropped.

PacBio's curated list (the source of havanno's names) goes stale. The case that surfaced it: RPGR's
p.Gly977_Glu983dup, in the ORF15 exon, is ClinVar-Benign (expert panel) yet curated PATHO, and flagged
6 of the 32 Atlanta carrier samples as RPGR "Pathogenic". The vendored list and ClinVar references
are public reference data; the sample regions here are synthetic.
"""
from pathlib import Path

import pytest

from repeatscope.api import main
from repeatscope.paraphase import (
    STATUS_PATHOGENIC_VARIANT,
    classify_gene,
    load_clinvar_reference,
    load_named_variant_positions,
)

UPLOADS = Path(__file__).resolve().parents[2] / "project" / "uploads"
STALE_RPGR = "RPGR_p.Gly977_Glu983dup_p.Gly977_Glu983dup_c.2919_2939dup"
STALE_RPGR_POSITION = ("chrX", 38286059, "C", "CTCCCCTTCTCCTTCCTCCCCT")


def _havanno(gene_key: str, *haplotype_variants: list[str]) -> dict:
    """Havanno annotations with one haplotype-phase call per variant list."""
    return {
        gene_key: {
            f"{gene_key}_hap{i}": {"num_pathogenic_variants": len(names), "variants": {"pathogenic": names}}
            for i, names in enumerate(haplotype_variants, start=1)
        }
    }


@pytest.fixture(scope="module")
def vendored_positions():
    return load_named_variant_positions(UPLOADS / "ptcp_variant_list.GRCh38.vcf")


@pytest.fixture(scope="module")
def vendored_rpgr_reference():
    return load_clinvar_reference(UPLOADS / "clinvar_rpgr.hg38.tsv")


def test_the_vendored_variant_list_maps_the_stale_rpgr_name_to_its_position(vendored_positions):
    assert vendored_positions[STALE_RPGR] == STALE_RPGR_POSITION


def test_the_vendored_rpgr_reference_has_clinvar_calling_that_variant_benign(vendored_rpgr_reference):
    entry = vendored_rpgr_reference[STALE_RPGR_POSITION]
    assert entry["clnsig"] == "Benign"
    assert entry["clnrevstat"] == "reviewed_by_expert_panel"


def test_without_the_cross_check_the_stale_variant_still_flags_pathogenic():
    finding = classify_gene("RPGR", {"rpgr": {"total_cn": 2}}, {}, _havanno("rpgr", [STALE_RPGR]))
    assert finding.status == STATUS_PATHOGENIC_VARIANT
    assert finding.variants == (STALE_RPGR,)


def test_the_stale_rpgr_variant_no_longer_flags_pathogenic_with_the_real_vendored_files(
    vendored_positions, vendored_rpgr_reference
):
    finding = classify_gene(
        "RPGR", {"rpgr": {"total_cn": 2}}, {}, _havanno("rpgr", [], [STALE_RPGR]),
        named_variant_positions=vendored_positions, clinvar_references={"rpgr": vendored_rpgr_reference},
    )
    assert finding.status == "Normal"  # falls through to the copy-number reading, 2 copies as expected
    assert not finding.variants


def test_a_variant_clinvar_does_not_call_benign_is_kept():
    positions = {"GENE_keep": ("chrX", 100, "A", "T"), "GENE_stale": ("chrX", 200, "A", "T")}
    reference = {
        ("chrX", 100, "A", "T"): {"clnsig": "Pathogenic"},
        ("chrX", 200, "A", "T"): {"clnsig": "Likely_benign"},
    }
    finding = classify_gene(
        "RPGR", {"rpgr": {"total_cn": 2}}, {}, _havanno("rpgr", ["GENE_keep", "GENE_stale"]),
        named_variant_positions=positions, clinvar_references={"rpgr": reference},
    )
    assert finding.status == STATUS_PATHOGENIC_VARIANT
    assert finding.variants == ("GENE_keep",)  # only the stale one is dropped


@pytest.mark.parametrize("clnsig", ["Uncertain_significance", "Conflicting_classifications_of_pathogenicity", "Pathogenic"])
def test_only_an_unambiguously_benign_classification_drops_a_variant(clnsig):
    positions = {"GENE_v": ("chrX", 100, "A", "T")}
    reference = {("chrX", 100, "A", "T"): {"clnsig": clnsig}}
    finding = classify_gene(
        "RPGR", {"rpgr": {"total_cn": 2}}, {}, _havanno("rpgr", ["GENE_v"]),
        named_variant_positions=positions, clinvar_references={"rpgr": reference},
    )
    assert finding.status == STATUS_PATHOGENIC_VARIANT


def test_a_name_missing_from_the_variant_list_or_clinvar_is_kept():
    reference = {("chrX", 100, "A", "T"): {"clnsig": "Benign"}}
    # not in the variant list at all
    kept = classify_gene(
        "RPGR", {"rpgr": {"total_cn": 2}}, {}, _havanno("rpgr", ["GENE_unmapped"]),
        named_variant_positions={"GENE_other": ("chrX", 100, "A", "T")}, clinvar_references={"rpgr": reference},
    )
    # in the list, but ClinVar has no entry for its position
    unknown = classify_gene(
        "RPGR", {"rpgr": {"total_cn": 2}}, {}, _havanno("rpgr", ["GENE_v"]),
        named_variant_positions={"GENE_v": ("chrX", 999, "A", "T")}, clinvar_references={"rpgr": reference},
    )
    assert kept.status == unknown.status == STATUS_PATHOGENIC_VARIANT


def test_the_hba_havanno_path_is_cross_checked_too():
    positions = {"HBA2_stale": ("chr16", 172000, "A", "G")}
    reference = {("chr16", 172000, "A", "G"): {"clnsig": "Benign"}}
    call = {"genotype": "aa/aa", "total_cn": 4}
    unchecked = classify_gene("HBA1/2", {"hba": call}, {}, _havanno("hba", ["HBA2_stale"]))
    checked = classify_gene(
        "HBA1/2", {"hba": call}, {}, _havanno("hba", ["HBA2_stale"]),
        named_variant_positions=positions, clinvar_references={"hba": reference},
    )
    assert unchecked.status == STATUS_PATHOGENIC_VARIANT
    assert checked.status == "Normal" and checked.zygosity == "Homozygous"  # the aa/aa genotype reading


def _rpgr_target():
    return {"gene": "RPGR", "chrom": "chrX", "start": 38283701, "stop": 38288763}


def test_the_api_loads_the_vendored_files_and_drops_the_stale_rpgr_variant():
    finding = main._gene_finding(_rpgr_target(), {"rpgr": {"total_cn": 2}}, {}, _havanno("rpgr", [STALE_RPGR]), [])
    assert finding.status == "Normal"


def test_the_api_keeps_flagging_when_the_reference_files_are_not_there(monkeypatch):
    monkeypatch.setattr(main, "NAMED_VARIANT_POSITIONS_VCF", str(UPLOADS / "does_not_exist.vcf"))
    monkeypatch.setattr(main, "CLINVAR_REFERENCE_TSV", {"rpgr": str(UPLOADS / "does_not_exist.tsv")})
    monkeypatch.setattr(main, "_named_variant_positions", None)
    monkeypatch.setattr(main, "_clinvar_references", {})
    finding = main._gene_finding(_rpgr_target(), {"rpgr": {"total_cn": 2}}, {}, _havanno("rpgr", [STALE_RPGR]), [])
    assert finding.status == STATUS_PATHOGENIC_VARIANT
