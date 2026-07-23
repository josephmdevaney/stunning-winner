from pathlib import Path

from repeatscope.strchive import classify, load_strchive_loci

REPO_ROOT = Path(__file__).resolve().parents[2]
STRCHIVE_JSON = REPO_ROOT / "project" / "uploads" / "STRchive-loci.json"


def test_load_strchive_loci_count_and_fields():
    loci = load_strchive_loci(STRCHIVE_JSON)
    assert len(loci) == 81
    by_id = {l.id: l for l in loci}

    htt = by_id["HD_HTT"]
    assert htt.gene == "HTT"
    assert htt.chrom == "chr4"
    assert htt.motif == "CAG"
    assert htt.benign_min == 6 and htt.benign_max == 26
    assert htt.pathogenic_min == 36

    # multi-value inheritance/evidence lists are joined with "/"
    sca2 = by_id["SCA2_ATXN2"]
    assert sca2.inheritance == "AD/AR"

    rfc1 = by_id["CANVAS_RFC1"]
    assert set(rfc1.pathogenic_motifs) == {"AAGGG", "ACAGG", "AAAGG", "CAGGG"}
    assert set(rfc1.benign_motifs) == {"AAAAG", "AAAGGG"}
    assert len(rfc1.unknown_motifs) == 14
    assert rfc1.is_motif_dependent is True
    assert htt.is_motif_dependent is False


def test_classify_pathogenic_intermediate_normal():
    loci = load_strchive_loci(STRCHIVE_JSON)
    htt = next(l for l in loci if l.id == "HD_HTT")

    patho = classify(htt, 18, 44)
    assert patho.status == "Full expansion"

    inter = classify(htt, 18, 30)
    assert inter.status == "Intermediate"

    normal = classify(htt, 15, 18)
    assert normal.status == "Normal"


def test_classify_no_genotype():
    loci = load_strchive_loci(STRCHIVE_JSON)
    htt = next(l for l in loci if l.id == "HD_HTT")
    result = classify(htt, None, None)
    assert result.status == "Normal"
    assert "No genotype call" in result.note


def test_classify_contraction_locus_direction():
    """Regression test for a real bug found against real HG002 TRGT output:
    CHNG3_MIR7-2 is a contraction disease (benign=4, pathogenic=3 -- *fewer* repeats
    is what's pathogenic). The naive expansion-style check (mx >= pathogenic_min)
    would call HG002's real benign 4/4 genotype "pathogenic" just because 4 >= 3."""
    loci = load_strchive_loci(STRCHIVE_JSON)
    mir7 = next(l for l in loci if l.id == "CHNG3_MIR7-2")
    assert mir7.pathogenic_min == 3 and mir7.benign_max == 4

    normal = classify(mir7, 4, 4)
    assert normal.status == "Normal"

    contracted = classify(mir7, 4, 3)
    assert contracted.status == "Full expansion"
    assert "Below pathogenic threshold" in contracted.note


def test_classify_vwa1_deviation_not_contraction_direction():
    """Regression test for a real bug found against real HG002 TRGT output:
    HMNR7_VWA1's pathogenic range (1-3) contains its own single benign value (2), per
    GeneReviews: "any deviation from 2 motifs is thought to be pathogenic." The
    contraction-direction check (call pathogenic when the smaller allele <=
    pathogenic_max=3) misfires here -- it flagged HG002's real, completely normal 2/2
    genotype as "Full expansion" just because 2 <= 3."""
    loci = load_strchive_loci(STRCHIVE_JSON)
    vwa1 = next(l for l in loci if l.id == "HMNR7_VWA1")
    assert vwa1.benign_min == vwa1.benign_max == 2
    assert vwa1.pathogenic_min == 1 and vwa1.pathogenic_max == 3

    # HG002's real call.
    normal = classify(vwa1, 2, 2)
    assert normal.status == "Normal"

    contraction = classify(vwa1, 2, 1)
    assert contraction.status == "Full expansion"
    assert "deviates from the benign value" in contraction.note

    expansion = classify(vwa1, 3, 2)
    assert expansion.status == "Full expansion"


def test_classify_ci_straddle_flagged_normal_to_pathogenic():
    """ACMG's 2026 NGS repeat-expansion points-to-consider statement: a call whose point
    estimate lands in one band but whose confidence interval reaches into another should
    carry a follow-up-confirmation caveat. HTT is CAG (3bp); a called-allele point
    estimate of 26 (Normal, <=26) with a real TRGT ALLR range wide enough to reach 38
    repeats (>=36, pathogenic) should flag the overlap rather than reporting a clean
    "Normal" with no caveat."""
    loci = load_strchive_loci(STRCHIVE_JSON)
    htt = next(l for l in loci if l.id == "HD_HTT")

    result = classify(htt, 18, 26, a1_range_bp=None, a2_range_bp=(78, 114))
    assert result.status == "Normal"
    assert "confidence interval" in result.note
    assert "Full expansion" in result.note
    assert "78-114 bp" in result.note
    assert "26-38 repeats" in result.note


def test_classify_ci_tight_no_straddle():
    """A narrow confidence interval that stays within the called band gets no caveat."""
    loci = load_strchive_loci(STRCHIVE_JSON)
    htt = next(l for l in loci if l.id == "HD_HTT")

    result = classify(htt, 18, 30, a1_range_bp=None, a2_range_bp=(88, 92))
    assert result.status == "Intermediate"
    assert "confidence interval" not in result.note


def test_classify_ci_absent_is_backward_compatible():
    """Omitting the range params (the pre-existing call signature) behaves exactly as
    before -- no caveat, since there's nothing to compare."""
    loci = load_strchive_loci(STRCHIVE_JSON)
    htt = next(l for l in loci if l.id == "HD_HTT")

    result = classify(htt, 18, 26)
    assert result.status == "Normal"
    assert "confidence interval" not in result.note


def test_classify_rfc1_motif_aware():
    """RFC1/CANVAS: pathogenicity depends on *which* motif expanded, not just length.
    A naive length-only check would misclassify a long benign-motif (AAAAG) expansion
    as pathogenic just because it clears RFC1's raw repeat-count threshold."""
    loci = load_strchive_loci(STRCHIVE_JSON)
    rfc1 = next(l for l in loci if l.id == "CANVAS_RFC1")

    # Real GV-01 finding: 90 repeats of AAAAG (benign motif) -- not disease-causing,
    # even though 90 falls inside RFC1's raw "intermediate" numeric range (11-200).
    benign_long = classify(rfc1, 11, 90, "AAAAG", "AAAAG")
    assert benign_long.status == "Normal"

    # Same length, but the pathogenic motif AAGGG instead -> genuinely intermediate.
    patho_motif_intermediate = classify(rfc1, 11, 90, "AAAAG", "AAGGG")
    assert patho_motif_intermediate.status == "Intermediate"

    # Pathogenic motif at full pathogenic length.
    patho_motif_full = classify(rfc1, 11, 500, "AAAAG", "AAGGG")
    assert patho_motif_full.status == "Full expansion"

    # A motif STRchive doesn't characterize as either benign or pathogenic -> uncertain,
    # not silently assumed benign.
    uncertain = classify(rfc1, 11, 90, "AAAAG", "ACAAG")
    assert uncertain.status == "Uncertain significance"


def test_classify_rfc1_ci_straddle_intermediate_to_pathogenic():
    """The CI-straddle caveat also applies to motif-dependent loci: RFC1's AAGGG is
    5bp; a point estimate of 11 (Intermediate, >=11) with a real ALLR range wide enough
    to reach 401 repeats (>=400, pathogenic) should flag the overlap."""
    loci = load_strchive_loci(STRCHIVE_JSON)
    rfc1 = next(l for l in loci if l.id == "CANVAS_RFC1")

    result = classify(rfc1, 11, 90, "AAGGG", "AAAAG", a1_range_bp=(55, 2005))
    assert result.status == "Intermediate"
    assert "confidence interval" in result.note
    assert "Full expansion" in result.note


def test_bclaf3_locus_fields():
    """FRAXG_BCLAF3 (LaFlamme, Clarkson, Ibanez et al. 2026): real GRCh38 coordinates
    from the paper's own Additional File 1 (Data S1). Thresholds are the empirical
    benign/intermediate boundary from 50,747 Genomics England control genomes (Data
    S3 -- control-allele counts drop from 24 at 23 repeats to 2 at 25) and the largest
    confirmed expansion measured in the study (1,238 repeats, Data S9); the 24-172
    range stays methylation-gated (see test_classify_methylation_dependent_locus)."""
    loci = load_strchive_loci(STRCHIVE_JSON)
    bclaf3 = next(l for l in loci if l.id == "FRAXG_BCLAF3")

    assert bclaf3.gene == "BCLAF3"
    assert bclaf3.chrom == "chrX"
    assert (bclaf3.start, bclaf3.stop) == (19990922, 19990973)
    assert bclaf3.motif == "CCG"
    assert bclaf3.benign_min == 6
    assert bclaf3.benign_max == 23
    assert bclaf3.intermediate_min == 24
    assert bclaf3.intermediate_max == 172
    assert bclaf3.pathogenic_min == 173
    assert bclaf3.pathogenic_max == 1238
    assert bclaf3.methylation_dependent is True


def test_classify_methylation_dependent_locus():
    """A methylation-dependent locus (FRAXG_BCLAF3) reports Uncertain significance --
    not Intermediate -- for a call in its intermediate band (24-172 repeats), since the
    paper's own long-read data leaves the exact hypermethylation threshold undefined in
    that range (bounded between 117 and 172 repeats). Above the pathogenic threshold
    (>=173, the largest confirmed measurement in the study), length alone is enough to
    call Full expansion -- the methylation_dependent flag only softens the intermediate
    band, not an explicit pathogenic threshold. A locus with no `methylation_dependent`
    flag (HTT) is unaffected and still reports plain Intermediate."""
    loci = load_strchive_loci(STRCHIVE_JSON)
    bclaf3 = next(l for l in loci if l.id == "FRAXG_BCLAF3")
    htt = next(l for l in loci if l.id == "HD_HTT")

    normal = classify(bclaf3, 15, 17)
    assert normal.status == "Normal"

    in_uncertain_band = classify(bclaf3, 15, 117)
    assert in_uncertain_band.status == "Uncertain significance"
    assert "hypermethylation" in in_uncertain_band.note
    assert "Intermediate" not in in_uncertain_band.note

    above_pathogenic_threshold = classify(bclaf3, 15, 482)
    assert above_pathogenic_threshold.status == "Full expansion"

    # methylation_dependent doesn't leak into loci that don't opt into it.
    assert htt.methylation_dependent is False
    inter = classify(htt, 18, 30)
    assert inter.status == "Intermediate"
