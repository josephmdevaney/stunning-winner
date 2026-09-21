"""Tests for MitoLR: seeded demo generation, real MitoRSaw VCF parsing, the
sample registry, and the API routes."""

from fastapi.testclient import TestClient

from mitolr import annotate, constraint, generate, gnomad, haplogroup, qc, reference, samples, vcf
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
    assert by_id["m3243AG"]["heteroplasmy"] == 58.0
    assert by_id["m3243AG"]["status"] == "Pathogenic"
    assert "MELAS" in by_id["m3243AG"]["note"]
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


# --- overlapping mtDNA reading frames ---------------------------------------

def test_features_for_pos_returns_both_genes_in_overlap():
    # MT-ATP8 (8366-8572) and MT-ATP6 (8527-9207) genuinely overlap in the real
    # rCRS annotation -- a variant at 8550 belongs to both, not just whichever
    # gene happens to be listed first in MT_FEATURES.
    names = [f["name"] for f in reference.features_for_pos(8550)]
    assert names == ["MT-ATP8", "MT-ATP6"]
    # MT-ND4L (10470-10766) / MT-ND4 (10760-12137) likewise overlap at 10760-10766.
    names2 = [f["name"] for f in reference.features_for_pos(10763)]
    assert names2 == ["MT-ND4L", "MT-ND4"]
    # A non-overlapping position still returns exactly one feature.
    assert [f["name"] for f in reference.features_for_pos(750)] == ["MT-RNR1"]


def test_locus_and_consequence_joins_overlapping_genes():
    locus, consequence, genes = annotate.locus_and_consequence(8550)
    assert locus == "MT-ATP8/MT-ATP6"
    assert consequence == "Protein-coding variant"
    assert genes == ["MT-ATP8", "MT-ATP6"]


def test_constraint_for_variant_merges_both_overlapping_genes():
    # Neither gene's constraint context should be dropped just because they
    # share this genomic position.
    result = constraint.constraint_for_variant(8550, ["MT-ATP8", "MT-ATP6"])
    assert result is not None
    symbols = {r["symbol"] for r in result["gene_rows"]}
    assert symbols == {"MT-ATP8", "MT-ATP6"}


# --- real (non-synthetic) large-deletion matching ---------------------------

def test_match_pathogenic_deletion_real_vcf_style_call():
    # A real MitoRSaw-style deletion call (REF longer than ALT, anchored) at the
    # exact KSS/CPEO common-deletion breakpoints (m.8483_13459del4977) must be
    # recognized via the real-VCF matching path -- not just the synthetic demo's
    # hardcoded override, which used to be the only place this variant matched.
    ref = "A" * 4978  # anchored: pos=8482 + 4978 bases -> del span 8483-13459
    hit = annotate.match_pathogenic(8482, ref, "A")
    assert hit is not None and hit["id"] == "kss_common_del"


def test_match_pathogenic_deletion_different_breakpoint_not_matched():
    # A deletion of similar size but different breakpoints is a related but
    # distinct event -- not silently matched to the curated common deletion.
    ref = "A" * 4978
    assert annotate.match_pathogenic(9000, ref, "A") is None


def test_deletion_span_wraps_across_circular_origin():
    # A deletion starting near the D-loop whose REF extends past position
    # 16,569 wraps back into 1..GENOME_LENGTH, rather than reporting a
    # nonsensical position beyond the genome's real length.
    del_start, del_stop = annotate.deletion_span(16560, "A" * 20, "A")
    assert del_start == 16561
    assert del_stop == 10


def test_variant_notation_deletion_wraps_label_and_span():
    position, change, is_del = annotate.variant_notation(16560, "A" * 20, "A")
    assert position == "m.16561_10del"
    assert change == "del"
    assert is_del is False  # 19 bp span, below the >100 bp large-deletion cutoff


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

def test_default_sample_resolves_to_configured_default():
    # DEFAULT_SAMPLE is environment-configurable (points at whichever real
    # sample is loaded); this only checks resolve(None) honors it, not which
    # sample that happens to be.
    assert samples.resolve(None) == samples.DEFAULT_SAMPLE


def test_overview_dispatch_per_sample():
    hg = samples.overview("HG002")
    assert hg["stats"]["variants_called"] == 15
    assert hg["stats"]["pathogenic_findings"] == 0
    assert hg["caller"].startswith("mitorsaw")
    # Real Haplogrep 3 classification if installed (scripts/setup_haplogrep.py);
    # otherwise honestly "Unclassified" rather than crashing when it's not
    # available -- either is a valid state for this environment.
    if haplogroup.available():
        assert hg["haplogroup"]["name"] != "Unclassified"
        assert hg["haplogroup"]["confidence"] is not None
    else:
        assert hg["haplogroup"]["name"] == "Unclassified"

    demo = samples.overview("MT-0417-P2")
    assert demo["stats"]["variants_called"] == 79


def test_mtdna_coverage_uses_real_whole_genome_source():
    # HG002 has a real per-base BAM QC fixture -> genome-wide mean from that.
    hg = samples.overview("HG002")
    assert hg["stats"]["mtdna_mean_coverage"] == 1194
    assert "BAM" in hg["stats"]["mtdna_coverage_source"]

    # MT-0417-P2 is the synthetic demo -> seeded figure, honestly labeled.
    demo = samples.overview("MT-0417-P2")
    assert demo["stats"]["mtdna_mean_coverage"] == 4820
    assert "synthetic" in demo["stats"]["mtdna_coverage_source"]
    assert demo["haplogroup"]["name"] == "H2a5b1"


def test_haplogroup_classify_hg002_when_installed():
    if not haplogroup.available():
        import pytest
        pytest.skip("Haplogrep 3 not installed -- run scripts/setup_haplogrep.py")
    result = haplogroup.classify("HG002", samples.SAMPLES["HG002"]["vcf"])
    assert result["name"] != "Unclassified"
    assert result["confidence"] is not None and 0 <= result["confidence"] <= 100
    assert "PhyloTree" in result["region"]


def test_haplogroup_classify_unavailable_reports_honestly(tmp_path, monkeypatch):
    # Force the "not installed" path regardless of whether the real classifier
    # happens to be set up in this environment.
    monkeypatch.setattr(haplogroup, "JAR_PATH", tmp_path / "does-not-exist.jar")
    haplogroup.classify.cache_clear()
    result = haplogroup.classify("HG002", samples.SAMPLES["HG002"]["vcf"])
    assert result == {
        "name": "Unclassified",
        "confidence": None,
        "region": "Haplogrep 3 not installed -- run `python backend/scripts/setup_haplogrep.py`",
    }
    haplogroup.classify.cache_clear()


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
    assert body["default"] == samples.DEFAULT_SAMPLE
    keys = [s["key"] for s in body["samples"]]
    assert "HG002" in keys and "MT-0417-P2" in keys


def test_overview_endpoint_for_hg002():
    body = client.get("/api/mito/overview?sample=HG002").json()
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


def test_acmg_ignores_heteroplasmic_frequency_never_sums_into_homoplasmic():
    # McCormick et al. validated BA1/BS1/PM2_supporting against AF_hom alone and
    # give no equivalent AF_het threshold. A variant never seen homoplasmic but
    # common heteroplasmically (AF_hom=0, AF_het=20%) must NOT qualify for BA1
    # just because AF_hom + AF_het crosses 1% -- that was the bug.
    result = gnomad.acmg_frequency_evidence(_fake_hit(0.0, af_het=0.20))
    assert [e["code"] for e in result] == ["—"]
    assert "no equivalent" in result[0]["note"] or "no ACMG frequency evidence" in result[0]["note"]


def test_acmg_real_302_ACC_insertion_yields_no_frequency_code_not_ba1():
    # Real gnomAD v3.1 row: chrM 302 A>ACC has AF_hom=0.0, AF_het=0.3493 -- the
    # exact case that used to be wrongly reported as BA1 (stand-alone benign)
    # under the old AF_hom+AF_het summation.
    hit = gnomad.frequency_for(302, "A", "ACC")
    assert hit is not None and hit["af_hom"] == 0.0 and hit["af_het"] > 0.3
    result = gnomad.acmg_frequency_evidence(hit)
    assert [e["code"] for e in result] == ["—"]


def test_acmg_pm2_supporting_still_applies_when_het_also_absent():
    # A variant that's rare at both homoplasmic AND heteroplasmic frequency
    # still legitimately gets PM2_supporting -- only the "never homoplasmic but
    # common heteroplasmic" case should be suppressed.
    result = gnomad.acmg_frequency_evidence(_fake_hit(0.0, af_het=0.0))
    assert [e["code"] for e in result] == ["PM2_supporting"]


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


# --- real per-variant heteroplasmy thresholds -------------------------------

def test_m8993_pathogenic_only_above_90_percent():
    status, note = annotate._classify_m8993(89.9, "NARP / Leigh syndrome")
    assert status == "Low-level heteroplasmy"
    assert "90%" in note and "Thorburn" in note

    status, note = annotate._classify_m8993(90.0, "NARP / Leigh syndrome")
    assert status == "Pathogenic"
    assert "90%" in note and "Thorburn" in note

    status, note = annotate._classify_m8993(2.0, "NARP / Leigh syndrome")
    assert status == "Not detected"


def test_m3243_disease_specific_bands():
    # Below 5%: not detected, but MIDD caveat still disclosed (can occur very low).
    status, note = annotate._classify_m3243(2.0, "MELAS / MIDD")
    assert status == "Not detected"
    assert "MIDD" in note and "Laloi-Michelin" in note

    # 5-49%: low-level, MIDD still possible even though below MELAS range.
    status, note = annotate._classify_m3243(38.0, "MELAS / MIDD")
    assert status == "Low-level heteroplasmy"
    assert "MIDD" in note

    # 50-89%: MELAS range.
    status, note = annotate._classify_m3243(58.0, "MELAS / MIDD")
    assert status == "Pathogenic"
    assert "MELAS" in note

    # >=90%: Leigh range.
    status, note = annotate._classify_m3243(92.0, "MELAS / MIDD")
    assert status == "Pathogenic"
    assert "Leigh" in note


def test_lhon_penetrance_caveat_only_when_pathogenic():
    # Homoplasmic (typical for LHON) -> Pathogenic, with the penetrance caveat.
    status, note = annotate._classify_lhon(100.0, "Leber hereditary optic neuropathy (LHON)")
    assert status == "Pathogenic"
    assert "penetrance" in note and "Yu-Wai-Man" in note

    # Below the generic low-level band -> no penetrance addendum tacked onto
    # a status where it wouldn't apply.
    status, note = annotate._classify_lhon(2.0, "Leber hereditary optic neuropathy (LHON)")
    assert status == "Not detected"
    assert "penetrance" not in note


def test_classify_pathogenic_heteroplasmy_dispatches_known_overrides_only():
    # A curated override variant uses its real per-variant threshold...
    status, note = annotate.classify_pathogenic_heteroplasmy("m8993TG", 95.0, "NARP / Leigh syndrome")
    assert status == "Pathogenic" and "Thorburn" in note

    # ...while any other panel variant without a sourced threshold falls back
    # to the generic band convention, unchanged.
    status, note = annotate.classify_pathogenic_heteroplasmy("m8344AG", 38.0, "MERRF")
    assert status == "Pathogenic" and "generic" in note


# --- tissue-of-origin caveat (Parikh et al. 2015) ---------------------------

def test_tissue_caveat_silent_without_a_known_tissue():
    assert annotate.tissue_caveat("Not detected", "m3243AG", None) is None


def test_tissue_caveat_silent_for_tissues_the_guidance_does_not_address():
    assert annotate.tissue_caveat("Not detected", "m3243AG", "muscle") is None


def test_tissue_caveat_silent_when_finding_is_already_positive():
    assert annotate.tissue_caveat("Pathogenic", "m3243AG", "blood") is None


def test_tissue_caveat_fires_for_blood_and_buccal_on_negative_or_low_finding():
    for tissue in ("blood", "buccal"):
        for status in ("Not detected", "Low-level heteroplasmy"):
            note = annotate.tissue_caveat(status, "m8993TG", tissue)
            assert note is not None
            assert "Parikh" in note and tissue in note


def test_tissue_caveat_m3243_specific_wording_mentions_urine():
    note = annotate.tissue_caveat("Not detected", "m3243AG", "blood")
    assert note is not None
    assert "urine" in note and "de Laat" in note

    # A different variant on the same tissue gets the general caveat, not the
    # m.3243A>G-specific urine wording.
    note = annotate.tissue_caveat("Not detected", "m8993TG", "blood")
    assert "urine" not in note


def test_classify_appends_tissue_caveat_for_negative_pathogenic_panel_call():
    result = annotate.classify(3243, "A", "G", 2.0, tissue="blood")
    assert result["status"] in ("Not detected", "Low-level heteroplasmy")
    assert "Parikh" in result["note"] and "urine" in result["note"]


def test_classify_omits_tissue_caveat_when_tissue_unknown():
    result = annotate.classify(3243, "A", "G", 2.0)
    assert "Parikh" not in result["note"]


def test_classify_tissue_has_no_effect_on_novel_variants():
    with_tissue = annotate.classify(99999, "A", "G", 100.0, tissue="blood")
    without_tissue = annotate.classify(99999, "A", "G", 100.0)
    assert with_tissue["note"] == without_tissue["note"]


# --- sample-level tissue registry --------------------------------------------

def test_hg002_meta_reports_blood_tissue():
    assert samples.meta("HG002")["tissue"] == "blood"


def test_demo_sample_meta_has_no_fabricated_tissue():
    assert samples.meta("MT-0417-P2")["tissue"] is None


def test_vcf_build_variants_accepts_tissue_without_affecting_novel_calls():
    plain = vcf.build_variants(HG002_VCF)
    with_tissue = vcf.build_variants(HG002_VCF, tissue="blood")
    # HG002's real calls are all novel (no curated-pathogenic hits), so the
    # tissue caveat -- which only applies to pathogenic-panel findings --
    # should leave every note unchanged.
    assert [v["note"] for v in plain] == [v["note"] for v in with_tissue]
