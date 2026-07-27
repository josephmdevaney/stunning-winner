from pathlib import Path

from repeatscope.annotate import (
    annotate_strchive_panel,
    iter_annotated_catalog,
    resolve_panel_locus_ids,
    simulate_pileup_reads,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
FIXTURE_VCF = REPO_ROOT / "backend" / "fixtures" / "synthetic_trgt.vcf"
STRCHIVE_JSON = REPO_ROOT / "project" / "uploads" / "STRchive-loci.json"
CATALOG_BED = REPO_ROOT / "project" / "uploads" / "repeat_catalog.hg38.bed"
GV01_VCF = REPO_ROOT / "backend" / "fixtures" / "gv_puretarget_panel_v2" / "trgt_vcfs" / "GV-01.trgt.vcf"
PURETARGET_V2_BED = REPO_ROOT / "backend" / "fixtures" / "panels" / "puretarget_v2.repeat_definitions.hg38.bed"
CARRIER_BED = REPO_ROOT / "backend" / "fixtures" / "panels" / "carrier_screen.repeat_definitions.hg38.bed"
HG002_VCF = REPO_ROOT / "backend" / "fixtures" / "hg002_real" / "HG002.strchive_loci.trgt.vcf.gz"


def _panel():
    return annotate_strchive_panel(FIXTURE_VCF, STRCHIVE_JSON)


def test_annotate_strchive_panel_covers_all_81_loci():
    panel = _panel()
    assert len(panel) == 81


def test_genotyped_loci_classified_correctly():
    panel = {l.id: l for l in _panel()}

    htt = panel["HD_HTT"]
    assert htt.genotyped
    assert {htt.a1, htt.a2} == {18, 44}
    assert htt.status == "Full expansion"
    assert htt.zygosity == "Heterozygous"

    fmr1 = panel["FXS_FMR1"]
    assert fmr1.a1 == fmr1.a2 == 29
    assert fmr1.zygosity == "Homozygous"
    assert fmr1.status == "Normal"

    atxn1 = panel["SCA1_ATXN1"]
    assert atxn1.status == "Intermediate"


def test_hg002_real_vwa1_normal_not_pathogenic():
    """End-to-end regression test: HG002's real HMNR7_VWA1 call is 2/2 -- the exact
    benign value. See test_strchive.py's test_classify_vwa1_deviation_not_contraction_direction
    for why a naive contraction-direction check misclassified this as pathogenic."""
    panel = {l.id: l for l in annotate_strchive_panel(HG002_VCF, STRCHIVE_JSON, sample="HG002")}
    vwa1 = panel["HMNR7_VWA1"]
    assert vwa1.genotyped
    assert (vwa1.a1, vwa1.a2) == (2, 2)
    assert vwa1.status == "Normal"


def test_ungenotyped_loci_marked_no_call_not_fabricated():
    panel = {l.id: l for l in _panel()}
    # a locus not present in the synthetic fixture at all
    other = panel["OPDM5_ABCD3"]
    assert other.genotyped is False
    assert other.a1 is None and other.a2 is None
    assert "No TRGT genotype call" in other.note


def test_simulate_pileup_reads_uses_real_depth_and_range():
    panel = {l.id: l for l in _panel()}
    htt = panel["HD_HTT"]
    reads = simulate_pileup_reads(htt)
    # total reads == real reported SD (16 + 15), not an arbitrary constant
    assert len(reads) == htt.a1_depth + htt.a2_depth
    cluster1_vals = [r["value"] for r in reads if r["cluster"] == 1]
    cluster2_vals = [r["value"] for r in reads if r["cluster"] == 2]
    assert all(10 <= v <= 26 for v in cluster1_vals)  # near a1=18
    assert all(35 <= v <= 53 for v in cluster2_vals)  # near a2=44


def test_simulate_pileup_empty_when_not_genotyped():
    panel = {l.id: l for l in _panel()}
    assert simulate_pileup_reads(panel["OPDM5_ABCD3"]) == []


def test_iter_annotated_catalog_matches_generic_and_strchive_regions():
    """Catalog regions and the sample's own VCF calls come from different tools/sources
    with different boundary-calling conventions, so this deliberately doesn't require an
    exact coordinate match -- see TrgtIndex.find_overlap. The fixture's synthetic call at
    chr1:44835-44867 (AAAT) falls entirely inside this real catalog's wider
    chr1:44804-44906 region for the same repeat, and should still be found by overlap."""
    regions = list(iter_annotated_catalog(CATALOG_BED, FIXTURE_VCF))
    assert len(regions) == 937_195

    by_coords = {(r.chrom, r.start, r.stop): r for r in regions}
    hit = by_coords[("chr1", 44804, 44906)]
    assert hit.genotyped
    assert {hit.a1, hit.a2} == {8, 9}

    miss = by_coords[("chr1", 102078, 102182)]
    assert miss.genotyped is False
    assert miss.a1 is None


def test_rfc1_real_gv01_call_not_misclassified_pathogenic():
    """End-to-end regression test: GV-01's real RFC1 call is a 90-unit expansion of the
    benign AAAAG motif. A length-only classifier would call this "Intermediate" (90
    falls in RFC1's raw 11-200 range) or worse; the correct, motif-aware answer is
    Normal, since AAAAG isn't disease-associated at any length."""
    panel = {l.id: l for l in annotate_strchive_panel(GV01_VCF, STRCHIVE_JSON, sample="GV-01")}
    rfc1 = panel["CANVAS_RFC1"]
    assert rfc1.genotyped
    assert {rfc1.a1, rfc1.a2} == {11, 90}
    assert rfc1.a1_motif == "AAAAG" and rfc1.a2_motif == "AAAAG"
    assert rfc1.status == "Normal"


def test_resolve_panel_locus_ids_handles_naming_drift():
    """Regression test: PacBio's real PureTarget V2 panel definition calls this locus
    "OPDM_ABCD3", but STRchive's own id for the same locus is "OPDM5_ABCD3" -- same
    coordinates, different name. A naive ID-set intersection would silently drop it
    from the panel's scope; coordinate-overlap fallback should still resolve it."""
    locus_ids = resolve_panel_locus_ids(PURETARGET_V2_BED, STRCHIVE_JSON)
    assert len(locus_ids) == 38
    assert "OPDM5_ABCD3" in locus_ids


def test_annotate_strchive_panel_scoped_to_puretarget_v2():
    locus_ids = resolve_panel_locus_ids(PURETARGET_V2_BED, STRCHIVE_JSON)
    panel = annotate_strchive_panel(GV01_VCF, STRCHIVE_JSON, sample="GV-01", locus_ids=locus_ids)
    assert len(panel) == 38
    assert all(l.genotyped for l in panel)


def test_hg002_real_fmr1_hemizygous_with_methylation():
    """A real hemizygous male chrX sample is available again: HG002's FXS_FMR1
    record was replaced with the true hemizygous call (GT="1", a single allele
    index) from a separate, genuinely genome-wide TRGT run for this donor
    (fixtures/hg002_real_genomewide/) -- swapped in specifically because that
    run also carries real 5mCpG methylation (TRGT's AM field), which the
    original padded-reference run for this fixture never captured (AM="." for
    every locus there). See fixtures/hg002_real/README.md's "FXS_FMR1 real
    methylation" section. A prior real carrier-panel run (bc2113, since
    replaced) demonstrated the same hemizygous-parsing code path; that
    regression case had no real fixture backing it for a while -- this
    replaces it with current real data rather than a since-removed one."""
    panel = {l.id: l for l in annotate_strchive_panel(HG002_VCF, STRCHIVE_JSON, sample="HG002")}
    fmr1 = panel["FXS_FMR1"]
    assert fmr1.genotyped
    assert fmr1.a1 == 29
    assert fmr1.a2 is None
    assert fmr1.zygosity == "Hemizygous"
    assert fmr1.status == "Normal"
    assert fmr1.a1_methylation == 0.12
    assert fmr1.a2_methylation is None
