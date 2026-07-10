"""Tests for MitoLR: seeded demo generation, real MitoRSaw VCF parsing, the
sample registry, and the API routes."""

from fastapi.testclient import TestClient

from mitolr import annotate, constraint, generate, gnomad, qc, samples, vcf
from repeatscope.api.main import app

client = TestClient(app)

HG002_VCF = samples.SAMPLES["HG002"]["vcf"]


# --- seeded demo sample -----------------------------------------------------

def test_demo_variants_deterministic_and_counted():
    a = generate.build_variants()
    b = generate.build_variants()
    assert [v["id"] for v in a] == [v["id"] for v in b]
    assert len(a) == 79
    assert sum(1 for v in a if v["origin"] == "pathogenic") == 21
    assert sum(1 for v in a if v["origin"] == "novel") == 58


def test_demo_narrative_overrides():
    by_id = {v["id"]: v for v in generate.build_variants()}
    assert by_id["m3243AG"]["heteroplasmy"] == 38.0
    assert by_id["m3243AG"]["status"] == "Pathogenic"
    assert by_id["kss_common_del"]["heteroplasmy"] == 11.0


# --- real MitoRSaw VCF parsing ---------------------------------------------

def test_hg002_vcf_parses():
    vs = vcf.build_variants(HG002_VCF)
    assert len(vs) == 15
    assert [v["pos"] for v in vs] == sorted(v["pos"] for v in vs)
    # HG002 has no curated-pathogenic calls; all are novel/benign.
    assert {v["origin"] for v in vs} == {"novel"}


def test_hg002_heteroplasmy_from_vaf_and_depth_from_dp():
    by_pos = {v["pos"]: v for v in vcf.build_variants(HG002_VCF)}
    # m.263A>G is homoplasmic (VAF 1.0) at DP 1108.
    assert by_pos[263]["heteroplasmy"] == 100.0
    assert by_pos[263]["depth"] == 1108
    assert by_pos[263]["change"] == "A>G"
    # m.302 A>AC is a control-region insertion at VAF 0.9369.
    assert by_pos[302]["heteroplasmy"] == 93.7
    assert by_pos[302]["change"].startswith("ins")


def test_hg002_locus_annotation():
    by_pos = {v["pos"]: v for v in vcf.build_variants(HG002_VCF)}
    assert by_pos[750]["locus"] == "MT-RNR1"
    assert by_pos[4769]["locus"] == "MT-ND2"
    assert by_pos[15326]["locus"] == "MT-CYB"
    assert by_pos[263]["locus"] == "D-loop"
    assert by_pos[750]["consequence"] == "rRNA variant"


def test_pileup_ratio_matches_allele_depths():
    by_pos = {v["pos"]: v for v in vcf.build_variants(HG002_VCF)}
    v = by_pos[302]  # AD = 70 ref, 1040 alt
    reads = vcf.simulate_pileup(v)
    assert len(reads) == 48
    frac = sum(1 for r in reads if r["kind"] == "mut") / len(reads)
    assert abs(frac - 1040 / (70 + 1040)) < 0.03


def test_variant_notation_forms():
    assert annotate.variant_notation(263, "A", "G") == ("m.263A>G", "A>G", False)
    pos, change, is_del = annotate.variant_notation(302, "A", "AC")
    assert change == "insC" and not is_del
    pos, change, is_del = annotate.variant_notation(3106, "CN", "C")
    assert change == "del" and not is_del  # 1 bp deletion, not the large-arc kind


# --- sample registry / dispatch --------------------------------------------

def test_default_sample_is_hg002():
    assert samples.resolve(None) == "HG002"


def test_overview_dispatch_per_sample():
    hg = samples.overview("HG002")
    assert hg["stats"]["variants_called"] == 15
    assert hg["stats"]["pathogenic_findings"] == 0
    assert hg["haplogroup"]["name"] == "Unclassified"
    assert hg["caller"].startswith("mitorsaw")

    demo = samples.overview("MT-0417-P2")
    assert demo["stats"]["variants_called"] == 79
    assert demo["haplogroup"]["name"] == "H2a5b1"


def test_qc_uses_bam_when_available():
    # HG002 has a precomputed BAM-stats fixture -> rich, real QC.
    hg = samples.sample_qc("HG002")
    assert len(hg["read_length_histogram"]) == 18  # real, not empty
    assert len(hg["coverage_track"]) > 100  # continuous per-base track
    labels = [s["label"] for s in hg["stats"]]
    assert "Mean coverage" in labels and "chrM reads mapped" in labels
    assert hg["coverage_note"].startswith("Per-base read depth")

    demo = samples.sample_qc("MT-0417-P2")
    assert len(demo["read_length_histogram"]) == 18
    assert len(demo["coverage_track"]) == 60


def test_bam_qc_suppresses_trna_strand_balance():
    hg = samples.sample_qc("HG002")
    by_region = {r["region"]: r for r in hg["regions"]}
    # tiny tRNA locus -> strand balance shown as "—" (midpoint approximation).
    assert by_region["MT-TF"]["strandBalance"] == "—"
    # large rRNA/protein locus keeps its computed value.
    assert by_region["MT-RNR2"]["strandBalance"] != "—"


def test_vcf_qc_fallback_without_bam(monkeypatch):
    # If the BAM fixture is absent, HG002 falls back to VCF-derived QC.
    monkeypatch.setitem(samples.SAMPLES["HG002"], "bam_qc", "/nonexistent.json")
    samples._bam_qc.cache_clear()
    hg = samples.sample_qc("HG002")
    assert hg["read_length_histogram"] == []
    samples._bam_qc.cache_clear()


def test_snv_pileup_uses_bam_indel_uses_caller():
    snv = samples.variant_detail("HG002", "chrM_263_A_G")
    assert "aligned HiFi BAM" in snv["pileup_source"]
    assert all(r["kind"] == "mut" for r in snv["pileup_reads"])  # homoplasmic G

    indel = samples.variant_detail("HG002", "chrM_302_A_AC")
    assert "haplotype-aware" in indel["pileup_source"]


# --- API routes -------------------------------------------------------------

def test_samples_endpoint():
    body = client.get("/api/mito/samples").json()
    assert body["default"] == "HG002"
    keys = [s["key"] for s in body["samples"]]
    assert "HG002" in keys and "MT-0417-P2" in keys


def test_overview_endpoint_defaults_to_hg002():
    body = client.get("/api/mito/overview").json()
    assert body["stats"]["variants_called"] == 15


def test_overview_endpoint_demo_sample():
    body = client.get("/api/mito/overview?sample=MT-0417-P2").json()
    assert body["stats"]["variants_called"] == 79


def test_variant_endpoint_and_404():
    vs = client.get("/api/mito/overview").json()["variants"]
    vid = vs[0]["id"]
    r = client.get(f"/api/mito/variant/{vid}")
    assert r.status_code == 200
    assert r.json()["variant"]["id"] == vid
    assert client.get("/api/mito/variant/nope").status_code == 404
    assert client.get("/api/mito/overview?sample=bogus").status_code == 404


def test_repeatscope_routes_still_present():
    assert client.get("/api/modes").status_code == 200


# --- gnomAD v3.1 mtDNA population frequency ---------------------------------

def test_gnomad_lookup_known_common_variant():
    hit = gnomad.frequency_for(263, "A", "G")
    assert hit is not None
    assert hit["pass_qc"] is True
    assert hit["af_hom"] > 0.9  # near-fixed in the population


def test_gnomad_lookup_miss_returns_none():
    # rCRS position 3107 is an 'N' padding base; MitoRSaw's anchored deletion
    # call here (CN>C at pos 3106) has no real counterpart in gnomAD's callset.
    assert gnomad.frequency_for(3106, "CN", "C") is None


def test_gnomad_flags_artifact_prone_homopolymer_site():
    hit = gnomad.frequency_for(302, "A", "AC")
    assert hit is not None
    assert hit["pass_qc"] is False
    assert "artifact_prone_site" in hit["filters"]


# --- gnomAD v4.1 mtDNA constraint --------------------------------------------

def test_gene_constraint_covers_rna_and_protein_genes():
    rows = constraint.gene_constraint_for("MT-RNR1")
    assert len(rows) == 1 and rows[0]["consequence"] == "RNA_variant"

    nd1_rows = constraint.gene_constraint_for("MT-ND1")
    assert {r["consequence"] for r in nd1_rows} == {"synonymous", "missense", "stop_gain"}


def test_gene_constraint_empty_for_dloop():
    # The control region isn't a gene and has no row in the constraint table.
    assert constraint.gene_constraint_for("D-loop") == []


def test_regional_constraint_matches_known_hotspot():
    # MT-CYB 14804-14857 is a flagged low-constraint sub-region.
    hit = constraint.regional_constraint_for(14820, "MT-CYB")
    assert hit is not None
    assert hit["start"] == 14804 and hit["end"] == 14857

    # Same gene, a position outside any flagged hotspot.
    assert constraint.regional_constraint_for(15000, "MT-CYB") is None


# --- real HG002 calls carry gnomAD + constraint annotation ------------------

def test_hg002_variants_carry_real_gnomad_and_constraint():
    by_pos = {v["pos"]: v for v in vcf.build_variants(HG002_VCF)}

    # m.750A>G: homoplasmic, MT-RNR1, common in gnomAD -> Common variant + constraint present.
    v750 = by_pos[750]
    assert v750["status"] == "Common variant"
    assert v750["gnomad"]["af_hom"] > 0.9
    assert v750["constraint"]["gene_rows"][0]["symbol"] == "MT-RNR1"

    # m.6800A>G: homoplasmic but real AF_hom << 1% -> Rare, not Common.
    v6800 = by_pos[6800]
    assert v6800["status"] == "Rare variant"
    assert v6800["gnomad"]["af_hom"] < 0.01

    # The 3106 CN>C deletion has no gnomAD counterpart at all.
    v3106 = by_pos[3106]
    assert v3106["gnomad"] is None
    assert v3106["status"] == "Rare variant"
    assert "Not found" in v3106["note"]

    # D-loop positions (263, 16304) aren't covered by the gene constraint table.
    assert by_pos[263]["constraint"] is None


def test_demo_sample_never_queries_real_population_data():
    # The synthetic sample's positions are generated, not observed -- a
    # coincidental match against gnomAD would misrepresent a fake call as real.
    for v in generate.build_variants():
        assert "gnomad" not in v
        assert "constraint" not in v


def test_classify_backward_compatible_without_gnomad_args():
    # Old call signature (no gnomad_hit/gnomad_checked) keeps old semantics.
    result = annotate.classify(99999, "A", "G", 100.0)
    assert result["status"] == "Common variant"
    assert "not evaluated" in result["note"]


# --- ACMG/AMP evidence tags (McCormick et al. 2020) ------------------------

def _fake_hit(af_hom, af_het=0.0, pass_qc=True):
    return {
        "af_hom": af_hom, "af_het": af_het, "ac_hom": 1, "ac_het": 0, "an": 1,
        "filters": ["PASS"] if pass_qc else ["artifact_prone_site"],
        "pass_qc": pass_qc, "max_observed_heteroplasmy": 1.0,
    }


def test_acmg_ba1_stand_alone_benign_above_1pct():
    result = gnomad.acmg_frequency_evidence(_fake_hit(0.02))
    assert [e["code"] for e in result] == ["BA1"]


def test_acmg_bs1_strong_benign_range():
    result = gnomad.acmg_frequency_evidence(_fake_hit(0.007))
    assert [e["code"] for e in result] == ["BS1"]


def test_acmg_pm2_supporting_below_threshold_or_absent():
    assert [e["code"] for e in gnomad.acmg_frequency_evidence(None)] == ["PM2_supporting"]
    assert [e["code"] for e in gnomad.acmg_frequency_evidence(_fake_hit(0.00001))] == ["PM2_supporting"]


def test_acmg_gap_between_thresholds_yields_no_evidence():
    # 0.1% falls between PM2_supporting's 0.002% ceiling and BS1's 0.5% floor --
    # a real gap in McCormick et al.'s thresholds, not a bug.
    assert gnomad.acmg_frequency_evidence(_fake_hit(0.001)) == []
    # 0.995% falls just above BS1's 0.99% ceiling and below BA1's 1% floor.
    assert gnomad.acmg_frequency_evidence(_fake_hit(0.00995)) == []


def test_acmg_pm5_supporting_different_allele_same_trna_position():
    # m.3243A>G (MELAS, MT-TL1) is in the curated panel. A different allele at the
    # same position raises suspicion (PM5_supporting); the exact same allele is
    # already covered by match_pathogenic, a stronger line of evidence, not this one.
    hit = annotate.pm5_supporting_candidate(3243, "C")
    assert hit is not None and hit["id"] == "m3243AG"
    assert annotate.pm5_supporting_candidate(3243, "G") is None  # exact match, excluded


def test_acmg_pm5_supporting_scoped_to_trna_rrna_not_missense():
    # m.8993T>G (NARP/Leigh, MT-ATP6) is a protein-coding Missense entry -- PM5_supporting
    # only applies to tRNA/rRNA panel entries (no codon-level consequence annotation
    # exists to support the protein-coding form of this criterion).
    assert annotate.pm5_supporting_candidate(8993, "C") is None


def test_hg002_acmg_evidence_end_to_end():
    by_pos = {v["pos"]: v for v in vcf.build_variants(HG002_VCF)}

    # m.263A>G: common (99.07% hom) -> BA1.
    assert [e["code"] for e in by_pos[263]["acmg_evidence"]] == ["BA1"]
    # m.8557G>A: real frequency ~0.86% -> BS1.
    assert [e["code"] for e in by_pos[8557]["acmg_evidence"]] == ["BS1"]
    # The 3106 CN>C deletion is absent from gnomAD entirely -> PM2_supporting.
    assert [e["code"] for e in by_pos[3106]["acmg_evidence"]] == ["PM2_supporting"]
    # m.6800A>G: real frequency ~0.03%, in the gap between thresholds -> no evidence.
    assert by_pos[6800]["acmg_evidence"] == []


def test_demo_sample_never_gets_acmg_evidence_key():
    for v in generate.build_variants():
        assert "acmg_evidence" not in v
