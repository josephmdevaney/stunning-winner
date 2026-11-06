from pathlib import Path

from repeatscope.trgt import TrgtCall, TrgtIndex, iter_trgt_records, sample_names

REPO_ROOT = Path(__file__).resolve().parents[2]
FIXTURE = REPO_ROOT / "backend" / "fixtures" / "synthetic_trgt.vcf"


def test_sample_names():
    assert sample_names(FIXTURE) == ["SYNTH-001"]


def test_iter_trgt_records_basic_fields():
    records = list(iter_trgt_records(FIXTURE))
    assert len(records) == 11
    htt = next(r for r in records if r.trid == "HD_HTT")
    assert htt.chrom == "chr4"
    assert htt.motifs == ["CAG"]


def test_allele_index_not_phase_order():
    """AL/MC etc. are indexed by declared-allele order, not by GT haplotype order."""
    records = list(iter_trgt_records(FIXTURE))
    cacna1a = next(r for r in records if r.trid == "SCA6_CACNA1A")
    # GT is "1|0": hap1 -> allele index 1, hap2 -> allele index 0.
    # (The reference span here rounds to 13 CTG units, coincidentally equal to the
    # configured a2, so index 0 ends up the 13-unit allele and index 1 the 11-unit one --
    # a naive "just take the first two MC values in field order" implementation would
    # get this backwards.)
    assert cacna1a.hap1_idx == 1 and cacna1a.hap2_idx == 0
    assert cacna1a.motif_count[0] == [13] and cacna1a.motif_count[1] == [11]
    units1, units2 = cacna1a.repeat_units()
    assert units1 == 11  # hap1 -> allele index 1
    assert units2 == 13  # hap2 -> allele index 0


def test_homozygous_call_units_and_depth():
    records = list(iter_trgt_records(FIXTURE))
    dmpk = next(r for r in records if r.trid == "DM1_DMPK")
    units1, units2 = dmpk.repeat_units()
    assert units1 == units2 == 5
    assert dmpk.total_depth == 47  # combined support for the single declared allele


def test_pathogenic_expansion_units():
    records = list(iter_trgt_records(FIXTURE))
    htt = next(r for r in records if r.trid == "HD_HTT")
    units1, units2 = htt.repeat_units()
    assert {units1, units2} == {18, 44}


def test_repeat_units_compound_motif_selects_target_not_sum():
    """Regression test for a real bug found against real HG002 TRGT output: HTT's VCF
    record has STRUC=(CAG)nCAACAG(CCG)n / MOTIFS=CAG,CCG, so MC reports one sub-count
    per segment (e.g. "18_8,25_8"). The clinically-relevant count is the CAG segment
    alone, not the sum of CAG+CCG -- summing would misclassify a benign CAG count as
    intermediate/pathogenic just because a downstream CCG tract adds to the total.

    Values below are taken verbatim from a real trgt v5.1.0 genotype call on HG002
    (GT=1/2, no haplotype calls REF, so the FORMAT arrays are exactly 2 slots wide --
    not a 3-slot REF+ALT1+ALT2 array)."""
    call = TrgtCall(
        chrom="chr4", pos=3074876, end=3074966, trid="HD_HTT",
        motifs=["CAG", "CCG"], struc="(CAG)nCAACAG(CCG)n", ref="N", alt=["<a>", "<b>"],
        hap1_idx=1, hap2_idx=2, phased=False, hemizygous=False,
        allele_length=[84, 105], allele_length_range=[(84, 88), (105, 110)],
        allele_support=[15, 13], motif_count=[[18, 8], [25, 8]],
        allele_purity=[0.97619, 0.980952], allele_methylation=[None, None],
    )
    # STRchive records this locus's pathogenic motif as "CAG" (matches MOTIFS[0] exactly).
    units1, units2 = call.repeat_units("CAG")
    assert (units1, units2) == (18, 25)

    # Without a target motif (e.g. generic catalog use), falls back to summing segments.
    sum1, sum2 = call.repeat_units()
    assert (sum1, sum2) == (26, 33)


def test_repeat_units_cyclic_rotation_matching():
    """STR motif reading-frame is arbitrary: STRchive may call FXN's motif "AAG" where
    the VCF's own MOTIFS says "GAA" -- these are the same repeat, just phase-shifted,
    and must still resolve to the same segment rather than falling back to summing.

    Values below are verbatim from a real trgt v5.1.0 call on HG002, with GT=2|1 (hap1
    is the *higher* original allele index) to also exercise non-ascending phase order:
    the FORMAT arrays are still ordered [value-for-idx1, value-for-idx2] regardless of
    which haplotype the phased GT lists first."""
    call = TrgtCall(
        chrom="chr9", pos=69037261, end=69037304, trid="FRDA_FXN",
        motifs=["A", "GAA"], struc="(A)n(GAA)n", ref="N", alt=["<a>", "<b>"],
        hap1_idx=2, hap2_idx=1, phased=True, hemizygous=False,
        allele_length=[52, 55], allele_length_range=[(50, 54), (54, 56)],
        allele_support=[12, 15], motif_count=[[23, 9], [26, 9]],
        allele_purity=[0.961538, 0.963636], allele_methylation=[None, None],
    )
    # hap1_idx=2 is the higher index -> gets slot 1 (bp=55); hap2_idx=1 -> slot 0 (bp=52).
    assert call.allele_lengths_bp == (55, 52)
    units1, units2 = call.repeat_units("AAG")
    assert (units1, units2) == (9, 9)


def test_methylation_present_and_absent():
    records = list(iter_trgt_records(FIXTURE))
    c9 = next(r for r in records if r.trid == "FTDALS1_C9orf72")
    m1, m2 = c9.allele_methylation_pct
    # one haplotype (the massive expansion) has methylation reported, the short allele does not
    assert (m1 is None) != (m2 is None)
    assert m1 == 0.71 or m2 == 0.71


def test_hemizygous_gt_single_allele_not_diploid_padded():
    """trgt's hemizygous calling mode (male chrX/chrY) emits GT="1", not "1/1" or
    "1/.": a real single-allele call, not a diploid call with a missing haplotype.
    FORMAT arrays for a hemizygous record have exactly one slot, not two."""
    call = TrgtCall(
        chrom="chrX", pos=147912049, end=147912111, trid="FXS_FMR1",
        motifs=["CGG", "AGG"], struc="<TR>", ref="N", alt=["<a>"],
        hap1_idx=1, hap2_idx=None, phased=False, hemizygous=True,
        allele_length=[2203], allele_length_range=[(95, 2535)],
        allele_support=[160], motif_count=[[731, 1]],
        allele_purity=[0.997732], allele_methylation=[0.67],
    )
    assert call.genotyped
    assert call.allele_lengths_bp == (2203, None)
    assert call.allele_depths == (160, None)
    assert call.allele_methylation_pct == (0.67, None)
    units1, units2 = call.repeat_units("CGG")
    assert (units1, units2) == (731, None)


def test_index_lookup_by_trid_and_overlap():
    index = TrgtIndex.from_vcf(FIXTURE)
    assert index.by_trid["HD_HTT"].chrom == "chr4"

    call = index.find_overlap("chr4", 3074876, 3074933)
    assert call is not None and call.trid == "HD_HTT"

    assert index.find_overlap("chr4", 999_999_999, 999_999_999 + 10) is None


def _synthetic_call(chrom: str, pos: int, end: int, trid: str) -> TrgtCall:
    return TrgtCall(
        chrom=chrom, pos=pos, end=end, trid=trid, motifs=["A"], struc="<TR>",
        ref="N", alt=["<a>"], hap1_idx=None, hap2_idx=None, phased=False, hemizygous=True,
        allele_length=[None], allele_length_range=[None], allele_support=[None],
        motif_count=[None], allele_purity=[None], allele_methylation=[None],
    )


def test_find_overlap_beyond_a_short_fixed_window():
    """Regression test: find_overlap must not assume a short fixed backward scan
    is enough to find a real overlap. A single wide call (early_call, pos=1000,
    end=2000) is the only one that overlaps a later query -- but 20 narrow calls
    sorted in between it and the query's own position push it more than 5 slots
    back from where a fixed-size scan would start, which previously caused a real
    overlap to be silently missed (e.g. TRExplorer v2.0's variation-cluster
    regions, which have runs of 1000+ consecutive nearby entries)."""
    calls = [_synthetic_call("chr1", 1000, 2000, "early_call")]
    for i in range(20):
        # Narrow, non-overlapping-with-query calls sorted after early_call by
        # start position, occupying index slots between it and the query.
        calls.append(_synthetic_call("chr1", 1001 + i, 1002 + i, f"filler_{i}"))
    index = TrgtIndex(calls)

    call = index.find_overlap("chr1", 1500, 1500)
    assert call is not None and call.trid == "early_call"


def test_dominant_motifs_real_rfc1_call():
    """Regression test using GV-01's real RFC1 call (real PacBio PureTarget panel data):
    GT=0/1, MOTIFS=AAGGG,ACAGG,AGGGC,AAGGC,AGAGG,AAAAG,AAAGG,AAGAG,AAAGGG,
    MC=0_0_0_0_0_11_0_0_0,0_0_0_0_0_90_0_0_0 -- both alleles' only nonzero segment is
    index 5 ("AAAAG", the benign reference motif), even though the second allele's
    repeat tract is long (90 units)."""
    call = TrgtCall(
        chrom="chr4", pos=39348424, end=39348483, trid="CANVAS_RFC1",
        motifs=["AAGGG", "ACAGG", "AGGGC", "AAGGC", "AGAGG", "AAAAG", "AAAGG", "AAGAG", "AAAGGG"],
        struc="<TR>", ref="N", alt=["<a>"],
        hap1_idx=0, hap2_idx=1, phased=False, hemizygous=False,
        allele_length=[59, 454], allele_length_range=[(55, 67), (430, 508)],
        allele_support=[180, 169],
        motif_count=[[0, 0, 0, 0, 0, 11, 0, 0, 0], [0, 0, 0, 0, 0, 90, 0, 0, 0]],
        allele_purity=[0.983333, 0.997802], allele_methylation=[None, 0.82],
    )
    m1, m2 = call.dominant_motifs()
    assert m1 == "AAAAG" and m2 == "AAAAG"

    units1, units2 = call.repeat_units()
    assert (units1, units2) == (11, 90)


# --- TRGT version comes from the VCF's own header ---------------------------------------------

def test_trgt_version_is_read_from_the_header(tmp_path):
    import gzip

    from repeatscope.trgt import trgt_version

    plain = tmp_path / "with_version.vcf"
    plain.write_text("##fileformat=VCFv4.2\n##trgtVersion=5.1.0-ec66463\n#CHROM\tPOS\n")
    assert trgt_version(plain) == "5.1.0"  # the commit suffix is dropped

    zipped = tmp_path / "with_version.vcf.gz"
    with gzip.open(zipped, "wt") as f:
        f.write("##fileformat=VCFv4.2\n##trgtVersion=3.0.0-d49a4dd\n#CHROM\tPOS\n")
    assert trgt_version(zipped) == "3.0.0"


def test_trgt_version_is_none_when_the_header_has_none(tmp_path):
    from repeatscope.trgt import trgt_version

    vcf = tmp_path / "no_version.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\n##trgtVersion=9.9.9\n")  # after #CHROM: not a header line
    assert trgt_version(vcf) is None
    assert trgt_version(FIXTURE) is None  # the synthetic fixture claims no version
