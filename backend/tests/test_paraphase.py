from pathlib import Path

from repeatscope.paraphase import (
    classify_gene,
    load_f8_inversion_calls,
    load_havanno_annotations,
    load_paraphase_calls,
    load_sawfish_calls,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
PARAPHASE_DIR = REPO_ROOT / "backend" / "fixtures" / "carrier_panel" / "paraphase"
SAMPLE = "m84036_250714_174336_s1.bc2113"
SAMPLE_2 = "m84036_250809_035451_s4.bc2192"


def _load():
    calls = load_paraphase_calls(PARAPHASE_DIR / f"{SAMPLE}.paraphase.json")
    f8_calls = load_f8_inversion_calls(PARAPHASE_DIR / f"{SAMPLE}.f8inversion.vcf")
    return calls, f8_calls


def _load_bc2192():
    calls = load_paraphase_calls(PARAPHASE_DIR / f"{SAMPLE_2}.paraphase.json")
    f8_calls = load_f8_inversion_calls(PARAPHASE_DIR / f"{SAMPLE_2}.f8inversion.vcf")
    havanno = load_havanno_annotations(PARAPHASE_DIR / f"{SAMPLE_2}.havanno.json")
    sawfish_calls = load_sawfish_calls(PARAPHASE_DIR / f"{SAMPLE_2}.sv.vcf")
    return calls, f8_calls, havanno, sawfish_calls


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
    assert finding.zygosity == "Homozygous"


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
    assert finding.status == "Pathogenic"
    assert "f8inv1" in finding.note
    assert finding.zygosity == "Hemizygous"


def test_classify_tnxb_returns_none_no_caller_mapped():
    """TNXB has no Paraphase key mapped (PacBio's own pipeline has no caller for it
    yet) -- callers must fall back to the generic 'Not applicable' placeholder."""
    calls, _ = _load()
    assert classify_gene("TNXB", calls, {}) is None


def test_classify_gene_missing_from_calls_returns_none():
    assert classify_gene("GBA", {}, {}) is None


def test_classify_f8_diploid_reference_gt_not_misread_as_inversion():
    """Real bug found via bc2192 (a female sample): the F8 inversion caller reports
    a diploid GT ("0/0") on chrX for female samples but a single hemizygous "0" for
    male samples. The classifier used to only recognize the hemizygous form as
    reference, so a diploid "0/0" fell through and was misclassified as an
    inversion. Guard against regressing that fix."""
    finding = classify_gene("F8", {}, {"f8inv1": {"gt": "0/0", "supp": 0}, "f8inv22": {"gt": "0/0", "supp": 0}})
    assert finding.status == "Normal"


def test_load_havanno_annotations_bc2192_has_all_genes():
    _, _, havanno, _ = _load_bc2192()
    assert set(havanno) == {"cyp21", "f8inv1", "f8inv22", "gba", "hba", "hbb", "rpgr", "smn1"}


def test_load_sawfish_calls_bc2192_finds_real_hba_deletion():
    _, _, _, sawfish_calls = _load_bc2192()
    dels = [c for c in sawfish_calls if c["svtype"] == "DEL"]
    assert len(dels) == 1
    d = dels[0]
    assert d["chrom"] == "chr16" and d["pos"] == 163481 and d["end"] == 163802
    assert d["svlen"] == 321 and d["gt"] == "0/1"


def test_classify_cyp21_bc2192_real_pathogenic_named_variants():
    """Real bc2192 data: 3 of 8 haplotype-phase calls (cyp21_hap2/hap5/hap7) all
    carry the same 9 named CAH-associated pathogenic variants -- a real finding
    the havanno annotation itself asserts, not re-derived here."""
    calls, f8_calls, havanno, sawfish_calls = _load_bc2192()
    finding = classify_gene("CYP21A2", calls, f8_calls, havanno, sawfish_calls)
    assert finding.status == "Pathogenic"
    assert "1 distinct variant set" in finding.note
    assert "9 variant(s)" in finding.note
    assert "CYP21A2_p.Gln319Ter_p.Gln319Ter_c.955C>T" in finding.note
    # 3 of 8 total physical copies at this locus carry the variant set (Paraphase
    # itself couldn't phase this complex locus into 2 chromosomes for this
    # sample -- alleles_final/haplotype_links are empty -- so "Heterozygous" here
    # means "some but not all copies", not a classic single-gene het/hom call).
    assert finding.zygosity == "Heterozygous"


def test_classify_gba_bc2192_real_pathogenic_named_variants():
    """Real bc2192 data: 4 of 7 haplotype-phase calls carry overlapping sets of
    named Gaucher-associated pathogenic variants (13 or 14 variants each)."""
    calls, f8_calls, havanno, sawfish_calls = _load_bc2192()
    finding = classify_gene("GBA", calls, f8_calls, havanno, sawfish_calls)
    assert finding.status == "Pathogenic"
    assert "4 of 7 haplotype-phase calls" in finding.note
    assert "Gauchian_Exon9-11_1" in finding.note
    assert finding.zygosity == "Heterozygous"


def test_classify_hba_bc2192_uncertain_with_real_sawfish_deletion():
    """Real bc2192 data: Paraphase itself couldn't resolve an HBA1/HBA2 genotype
    (genotype is null), but sawfish independently found a real heterozygous 321bp
    deletion overlapping the region. The honest answer is Uncertain -- reporting
    the real structural evidence -- not a guessed Normal or Pathogenic."""
    calls, f8_calls, havanno, sawfish_calls = _load_bc2192()
    finding = classify_gene("HBA1/2", calls, f8_calls, havanno, sawfish_calls)
    assert finding.status == "Uncertain significance"
    assert "chr16:163481-163802" in finding.note
    assert "321 bp" in finding.note
    assert finding.zygosity == "Heterozygous"


def test_classify_f8_bc2192_real_no_inversion():
    calls, f8_calls, havanno, sawfish_calls = _load_bc2192()
    finding = classify_gene("F8", calls, f8_calls, havanno, sawfish_calls)
    assert finding.status == "Normal"


def test_classify_smn1_bc2192_real_normal():
    calls, f8_calls, havanno, sawfish_calls = _load_bc2192()
    finding = classify_gene("SMN1/2", calls, f8_calls, havanno, sawfish_calls)
    assert finding.status == "Normal"
    assert "4 SMN1 copies" in finding.note


def test_classify_rpgr_bc2192_real_normal():
    calls, f8_calls, havanno, sawfish_calls = _load_bc2192()
    finding = classify_gene("RPGR", calls, f8_calls, havanno, sawfish_calls)
    assert finding.status == "Normal"


def test_classify_hba_prefers_named_havanno_variant_over_genotype():
    """hba isn't exercised by either real fixture's havanno data (both report 0
    pathogenic variants there), but the same named-variant path used for
    cyp21/gba should still take priority over Paraphase's own genotype field if a
    future sample's havanno annotation does flag one."""
    calls = {"hba": {"genotype": "aa/aa"}}
    havanno = {"hba": {"hba_hba1hap1": {"num_pathogenic_variants": 1, "variants": {"pathogenic": ["HBA_fake_del"]}}}}
    finding = classify_gene("HBA1/2", calls, {}, havanno, [])
    assert finding.status == "Pathogenic"
    assert "HBA_fake_del" in finding.note


def test_havanno_zygosity_homozygous_when_all_copies_affected():
    calls = {"cyp21": {"total_cn": 2, "two_copy_haplotypes": []}}
    havanno = {
        "cyp21": {
            "cyp21_hap1": {"num_pathogenic_variants": 1, "variants": {"pathogenic": ["V1"]}},
            "cyp21_hap2": {"num_pathogenic_variants": 1, "variants": {"pathogenic": ["V1"]}},
        }
    }
    finding = classify_gene("CYP21A2", calls, {}, havanno, [])
    assert finding.status == "Pathogenic"
    assert finding.zygosity == "Homozygous"
