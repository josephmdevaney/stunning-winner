from pathlib import Path

from repeatscope.paraphase import classify_gene, load_f8_inversion_calls, load_paraphase_calls

REPO_ROOT = Path(__file__).resolve().parents[2]
PARAPHASE_DIR = REPO_ROOT / "backend" / "fixtures" / "carrier_panel" / "paraphase"
SAMPLE = "m84036_250714_174336_s1.bc2113"


def _load():
    calls = load_paraphase_calls(PARAPHASE_DIR / f"{SAMPLE}.paraphase.json")
    f8_calls = load_f8_inversion_calls(PARAPHASE_DIR / f"{SAMPLE}.f8inversion.vcf")
    return calls, f8_calls


def test_load_paraphase_calls_has_all_genes():
    calls, _ = _load()
    assert set(calls) == {"gba", "cyp21", "f8inv1", "f8inv22", "hba", "hbb", "rpgr", "smn1"}


def test_load_f8_inversion_calls_real_bc2113_no_inversion():
    _, f8_calls = _load()
    assert set(f8_calls) == {"f8inv1", "f8inv22"}
    # Real dedicated F8 inversion caller output: GT=0 (reference), SUPP=0 (no reads
    # support the inversion breakpoint) at both probes -- no inversion detected.
    assert f8_calls["f8inv1"]["gt"] == "0" and f8_calls["f8inv1"]["supp"] == 0
    assert f8_calls["f8inv22"]["gt"] == "0" and f8_calls["f8inv22"]["supp"] == 0


def test_classify_hba_normal_genotype():
    calls, _ = _load()
    finding = classify_gene("HBA1/2", calls, {})
    assert finding.status == "Normal"
    assert "aa/aa" in finding.note


def test_classify_smn1_normal_copy_number():
    calls, _ = _load()
    finding = classify_gene("SMN1/2", calls, {})
    assert finding.status == "Normal"
    assert "2 SMN1 copies" in finding.note


def test_classify_gba_normal_no_fusion():
    calls, _ = _load()
    finding = classify_gene("GBA", calls, {})
    assert finding.status == "Normal"
    assert "4 total" in finding.note


def test_classify_cyp21_atypical_copy_number_not_guessed_normal_or_pathogenic():
    """Real bc2113 data: total_cn=3 (not the expected 4), but gene_cn is null --
    Paraphase itself couldn't resolve which gene lost a copy. The correct, honest
    answer is "Uncertain significance", not a guessed Normal or Pathogenic."""
    calls, _ = _load()
    finding = classify_gene("CYP21A2", calls, {})
    assert finding.status == "Uncertain significance"
    assert "3 total" in finding.note


def test_classify_f8_no_inversion_detected():
    _, f8_calls = _load()
    finding = classify_gene("F8", {}, f8_calls)
    assert finding.status == "Normal"
    assert "No intron 1/22 inversion" in finding.note


def test_classify_f8_inversion_detected():
    finding = classify_gene("F8", {}, {"f8inv1": {"gt": "1", "supp": 42}})
    assert finding.status == "Full expansion"
    assert "f8inv1" in finding.note


def test_classify_tnxb_returns_none_no_caller_mapped():
    """TNXB has no Paraphase key mapped (PacBio's own pipeline has no caller for it
    yet) -- callers must fall back to the generic 'Not applicable' placeholder."""
    calls, _ = _load()
    assert classify_gene("TNXB", calls, {}) is None


def test_classify_gene_missing_from_calls_returns_none():
    assert classify_gene("GBA", {}, {}) is None
