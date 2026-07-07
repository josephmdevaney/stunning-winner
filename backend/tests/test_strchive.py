from pathlib import Path

from repeatscope.strchive import classify, load_strchive_loci

REPO_ROOT = Path(__file__).resolve().parents[2]
STRCHIVE_JSON = REPO_ROOT / "project" / "uploads" / "STRchive-loci.json"


def test_load_strchive_loci_count_and_fields():
    loci = load_strchive_loci(STRCHIVE_JSON)
    assert len(loci) == 80
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
