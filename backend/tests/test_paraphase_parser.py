"""Parsing of interpret_paraphase.py's summary CSV (paraphase_annotate/parser.py). Synthetic rows only."""

from paraphase_annotate import enrich, parser

HEADER = (
    "Region,Region_median_depth,Gene,Total_haplotypes,Num_haplotypes_with_pathogenic_variants,"
    "Num_haplotypes_without_pathogenic_variants,Useful_information_from_Paraphase_JSON\n"
)


def _summary(tmp_path, *rows):
    path = tmp_path / "S_summary.csv"
    path.write_text(HEADER + "".join(r + "\n" for r in rows))
    return parser.parse_summary(path)


def test_median_depth_is_read_whether_paraphase_wrote_it_whole_or_with_a_decimal(tmp_path):
    rows = _summary(
        tmp_path,
        "rccx,50,CYP21A2,4,2,2,info",          # a whole number, as some samples have it
        "smn1,39.0,SMN1,2,0,2,info",           # the same, with a decimal point
        "pms2,41.5,PMS2,4,0,4,info",           # a real fraction
        "hba,NA,HBA1,2,0,2,",                  # no depth
    )
    assert [r["region_median_depth"] for r in rows] == [50, 39, 41.5, None]
    assert isinstance(rows[1]["region_median_depth"], int)        # 39.0 is shown as 39, not 39.0


def test_haplotype_counts_stay_integers_and_flag_a_region_with_no_functional_haplotype(tmp_path):
    row = _summary(tmp_path, "smn1,60,SMN1,2,2,0,info")[0]
    assert (row["total_haplotypes"], row["n_haplotypes_with_pathogenic"], row["n_haplotypes_without_pathogenic"]) == (2, 2, 0)
    assert row["zero_functional_haplotypes"] is True and row["info"] == "info"


DETAILED_HEADER = (
    "Region,Gene,Transcript,Transcript_start,Transcript_end,Haplotype,Haplotype_boundary,Boundary_check,"
    "Pathogenic_variants_annotated\n"
)
THREE_VARIANTS = (
    "1:g.100A>G=ClinVar-id1;Pathogenic;:ENST1;missense_variant;MODERATE;,"
    "1:g.200_201del=ENST1;frameshift_variant;HIGH;,"
    "1:g.300C>T=ClinVar-id2;Pathogenic;:ENST1;stop_gained;HIGH;"
)


def test_every_pathogenic_variant_of_a_haplotype_is_kept_not_just_the_first(tmp_path, monkeypatch):
    # the script writes the three variants comma-joined in one unquoted cell; csv reads them as three columns
    monkeypatch.setattr(enrich, "hgvs_for", lambda variant, transcript: (None, None))
    monkeypatch.setattr(enrich, "clinvar_classification_for", lambda clinvar_id: (None, None, None))
    path = tmp_path / "S_detailed.csv"
    path.write_text(
        DETAILED_HEADER
        + f"GBA,GBA1,ENST1,1,2,GBA_gene2hap1,10-20,fully_contained,{THREE_VARIANTS}\n"
        + "GBA,GBA1,ENST1,1,2,GBA_gene1hap1,10-20,fully_contained,None annotated\n"
        + "GBA,GBA1,ENST1,1,2,GBA_gene1hap2,10-20,fully_contained,known_nonfunctional\n"
    )
    by_haplotype = {r["haplotype"]: r for r in parser.parse_detailed(path)}
    assert [v["variant"] for v in by_haplotype["GBA_gene2hap1"]["pathogenic_variants"]] == ["1:g.100A>G", "1:g.200_201del", "1:g.300C>T"]
    assert by_haplotype["GBA_gene1hap1"]["pathogenic_variants"] == []
    assert by_haplotype["GBA_gene1hap2"]["known_nonfunctional"] is True and by_haplotype["GBA_gene1hap2"]["pathogenic_variants"] == []


def test_a_summary_info_value_with_commas_is_not_cut_off(tmp_path):
    rows = _summary(tmp_path, "rccx,39.0,CYP21A2,4,2,2,Phasing success: False; alleles: a,b,c")
    assert rows[0]["info"] == "Phasing success: False; alleles: a,b,c"
