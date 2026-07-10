from pathlib import Path

from repeatscope.paraphase import (
    classify_gene,
    load_f8_inversion_calls,
    load_havanno_annotations,
    load_paraphase_calls,
    load_sawfish_calls,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
# Real per-sample layout: paraphase.json/sv.vcf.gz sit in a per-sample subdirectory
# ({sample}_paraphase/, {sample}_sawfish/); f8inversion.vcf/havanno.json sit flat.
PARAPHASE_DIR = REPO_ROOT / "backend" / "fixtures" / "carrier_panel_v2"
SAMPLE = "m84036_250809_035451_s4.bc2113"
SAMPLE_2 = "m84036_250809_035451_s4.bc2192"
SAMPLE_3 = "m84036_250809_035451_s4.bc2108"
SAMPLE_4 = "m84036_250809_035451_s4.bc2100"
SAMPLE_5 = "m84036_250809_035451_s4.bc2103"


def _load():
    calls = load_paraphase_calls(PARAPHASE_DIR / f"{SAMPLE}_paraphase" / f"{SAMPLE}.paraphase.json")
    f8_calls = load_f8_inversion_calls(PARAPHASE_DIR / f"{SAMPLE}.f8inversion.vcf")
    return calls, f8_calls


def _load_bc2192():
    calls = load_paraphase_calls(PARAPHASE_DIR / f"{SAMPLE_2}_paraphase" / f"{SAMPLE_2}.paraphase.json")
    f8_calls = load_f8_inversion_calls(PARAPHASE_DIR / f"{SAMPLE_2}.f8inversion.vcf")
    havanno = load_havanno_annotations(PARAPHASE_DIR / f"{SAMPLE_2}.havanno.json")
    sawfish_calls = load_sawfish_calls(PARAPHASE_DIR / f"{SAMPLE_2}_sawfish" / f"{SAMPLE_2}.sv.vcf.gz")
    return calls, f8_calls, havanno, sawfish_calls


def _load_bc2108():
    calls = load_paraphase_calls(PARAPHASE_DIR / f"{SAMPLE_3}_paraphase" / f"{SAMPLE_3}.paraphase.json")
    f8_calls = load_f8_inversion_calls(PARAPHASE_DIR / f"{SAMPLE_3}.f8inversion.vcf")
    havanno = load_havanno_annotations(PARAPHASE_DIR / f"{SAMPLE_3}.havanno.json")
    return calls, f8_calls, havanno


def _load_sample(sample: str):
    calls = load_paraphase_calls(PARAPHASE_DIR / f"{sample}_paraphase" / f"{sample}.paraphase.json")
    f8_calls = load_f8_inversion_calls(PARAPHASE_DIR / f"{sample}.f8inversion.vcf")
    havanno = load_havanno_annotations(PARAPHASE_DIR / f"{sample}.havanno.json")
    return calls, f8_calls, havanno


def test_load_paraphase_calls_has_all_genes():
    calls, _ = _load()
    assert set(calls) == {"gba", "cyp21", "f8inv1", "f8inv22", "hba", "hbb", "rpgr", "smn1"}


def test_load_f8_inversion_calls_real_bc2113_no_inversion():
    _, f8_calls = _load()
    assert set(f8_calls) == {"f8inv1", "f8inv22"}
    # Real dedicated F8 inversion caller output for this sample: GT="." (no call,
    # not reference "0") at both probes, SUPP=0 (no reads support an inversion
    # breakpoint either way).
    assert f8_calls["f8inv1"]["gt"] == "." and f8_calls["f8inv1"]["supp"] == 0
    assert f8_calls["f8inv22"]["gt"] == "." and f8_calls["f8inv22"]["supp"] == 0


def test_classify_hba_bc2113_real_unresolved_genotype():
    """Real bc2113 data: Paraphase did not resolve an HBA1/HBA2 genotype at all for
    this sample -- the honest answer is Uncertain, not a guessed Normal."""
    calls, _ = _load()
    finding = classify_gene("HBA1/2", calls, {})
    assert finding.status == "Uncertain significance"
    assert "did not report a resolved HBA1/HBA2 genotype" in finding.note


def test_classify_smn1_bc2113_real_unresolved_copy_number():
    calls, _ = _load()
    finding = classify_gene("SMN1/2", calls, {})
    assert finding.status == "Uncertain significance"
    assert "did not report an SMN1 copy number" in finding.note


def test_classify_gba_bc2113_real_unresolved_copy_number():
    calls, _ = _load()
    finding = classify_gene("GBA", calls, {})
    assert finding.status == "Uncertain significance"
    assert "did not report a copy number for this region" in finding.note


def test_classify_cyp21_bc2113_real_unresolved_copy_number_not_guessed():
    """Real bc2113 data: Paraphase reported no copy number at all for this region
    (unlike bc2192/bc2108/bc2100/bc2103, where it resolves cleanly). The correct,
    honest answer is "Uncertain significance", not a guessed Normal or Pathogenic."""
    calls, _ = _load()
    finding = classify_gene("CYP21A2", calls, {})
    assert finding.status == "Uncertain significance"
    assert "did not report a copy number for this region" in finding.note


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


def test_classify_cyp21_bc2108_real_heterozygous_carrier_normal_copy_number():
    """Real bc2108 data: unlike bc2192's atypical elevated total_cn=8, this sample
    has a completely normal total_cn=4 at this locus, with 2 of 4 haplotype-phase
    calls carrying the same 9 named CAH-associated pathogenic variants -- a clean
    textbook heterozygous carrier finding."""
    calls, f8_calls, havanno = _load_bc2108()
    finding = classify_gene("CYP21A2", calls, f8_calls, havanno)
    assert finding.status == "Pathogenic"
    assert finding.zygosity == "Heterozygous"
    assert "2 of 4 haplotype-phase calls" in finding.note


def test_classify_gba_bc2108_real_heterozygous_carrier_normal_copy_number():
    """Real bc2108 data: total_cn=4 (normal), 2 of 4 haplotype-phase calls carry
    the same 14 named Gaucher-associated pathogenic variants."""
    calls, f8_calls, havanno = _load_bc2108()
    finding = classify_gene("GBA", calls, f8_calls, havanno)
    assert finding.status == "Pathogenic"
    assert finding.zygosity == "Heterozygous"
    assert "2 of 4 haplotype-phase calls" in finding.note


def test_classify_hba_bc2108_real_alpha_thalassemia_genotype_not_overclaimed():
    """Real bc2108 data: Paraphase itself resolved a genotype of '-a/aa' (one
    single-alpha-globin-gene deletion, the classic silent-carrier alpha
    thalassemia allele) and its own sv_called field names a 3.7kb deletion
    (3p7del) at chr16:172240-176044. The honest answer is Uncertain -- the
    genotype and structural evidence are surfaced verbatim, but a specific
    clinical-significance call isn't fabricated for a non-aa/aa genotype."""
    calls, f8_calls, havanno = _load_bc2108()
    finding = classify_gene("HBA1/2", calls, f8_calls, havanno)
    assert finding.status == "Uncertain significance"
    assert "-a/aa" in finding.note


def test_classify_f8_bc2108_real_no_inversion():
    calls, f8_calls, havanno = _load_bc2108()
    finding = classify_gene("F8", calls, f8_calls, havanno)
    assert finding.status == "Normal"


def test_classify_smn1_bc2108_real_normal():
    calls, f8_calls, havanno = _load_bc2108()
    finding = classify_gene("SMN1/2", calls, f8_calls, havanno)
    assert finding.status == "Normal"


def test_classify_cyp21_bc2100_real_single_haplotype_carrier():
    """Real bc2100 data: unlike bc2108's 2-of-4 pattern, only 1 of 4 haplotype-phase
    calls carries the 9-variant CAH set here -- still resolves to Heterozygous
    since normal total_cn=4 means "not all copies affected" either way."""
    calls, f8_calls, havanno = _load_sample(SAMPLE_4)
    finding = classify_gene("CYP21A2", calls, f8_calls, havanno)
    assert finding.status == "Pathogenic"
    assert finding.zygosity == "Heterozygous"
    assert "1 of 4 haplotype-phase calls" in finding.note


def test_classify_gba_bc2100_real_heterozygous_carrier():
    calls, f8_calls, havanno = _load_sample(SAMPLE_4)
    finding = classify_gene("GBA", calls, f8_calls, havanno)
    assert finding.status == "Pathogenic"
    assert finding.zygosity == "Heterozygous"


def test_classify_hba_bc2100_real_normal_aa_aa():
    """Real bc2100 data: unlike bc2108's real '-a/aa' alpha-thalassemia finding,
    this sample's HBA1/HBA2 genotype resolved cleanly to aa/aa."""
    calls, f8_calls, havanno = _load_sample(SAMPLE_4)
    finding = classify_gene("HBA1/2", calls, f8_calls, havanno)
    assert finding.status == "Normal"
    assert finding.zygosity == "Homozygous"


def test_classify_cyp21_bc2103_real_single_haplotype_carrier():
    calls, f8_calls, havanno = _load_sample(SAMPLE_5)
    finding = classify_gene("CYP21A2", calls, f8_calls, havanno)
    assert finding.status == "Pathogenic"
    assert finding.zygosity == "Heterozygous"


def test_classify_gba_bc2103_real_heterozygous_carrier():
    calls, f8_calls, havanno = _load_sample(SAMPLE_5)
    finding = classify_gene("GBA", calls, f8_calls, havanno)
    assert finding.status == "Pathogenic"
    assert finding.zygosity == "Heterozygous"


def test_classify_hba_bc2103_real_normal_aa_aa():
    calls, f8_calls, havanno = _load_sample(SAMPLE_5)
    finding = classify_gene("HBA1/2", calls, f8_calls, havanno)
    assert finding.status == "Normal"
    assert finding.zygosity == "Homozygous"
