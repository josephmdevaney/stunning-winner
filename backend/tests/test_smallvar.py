"""Tests for the small-variant suite (backend/smallvar).

Inputs are small synthetic files in the real formats (a VCF, a ClinVar table, slivar gnotate
zips), plus the committed public HG002 chr22 slice. Nothing here depends on the real
annotation resources or on a real sample."""

import csv
import gzip
import json
import struct
import time
import zipfile
from array import array
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from repeatscope.api.main import app
from smallvar import index as index_mod
from smallvar import samples, vcf
from smallvar.clinvar import COLUMNS, ClinvarTable, classification_group, review_stars
from smallvar.genes import GeneIndex
from smallvar.gnotate import Gnotate, encode_allele
from smallvar.normalize import canon_chrom, chrom_order, trim_alleles, variant_key, variant_kind

FIXTURE_VCF = (
    Path(__file__).resolve().parents[1] / "fixtures" / "hg002_smallvar" / "HG002.GRCh38.giab_v4.2.1.chr22.vcf.gz"
)
ACMG_TSV = Path(__file__).resolve().parents[1] / "smallvar" / "reference" / "acmg_sf_v3_3.tsv"


# ---------------------------------------------------------------- normalisation
def test_canonical_chromosomes():
    assert canon_chrom("7") == "chr7" and canon_chrom("chr7") == "chr7"
    assert canon_chrom("MT") == "chrM" and canon_chrom("chrMT") == "chrM" and canon_chrom("M") == "chrM"
    assert [chrom_order(c) for c in ("chr1", "22", "X", "chrY", "chrM", "chrUn")] == [1, 22, 23, 24, 25, 99]


def test_trimming_gives_the_parsimonious_form():
    assert trim_alleles(100, "A", "G") == (100, "A", "G")
    assert trim_alleles(100, "A", "AT") == (101, "", "T")  # insertion: padding base removed
    assert trim_alleles(100, "TC", "T") == (101, "C", "")  # deletion
    assert trim_alleles(100, "GGAT", "TGCTG") == (100, "GGAT", "TGCTG")  # a delins is left alone
    assert trim_alleles(10, "CAG", "CAGCAG") == (10, "", "CAG")  # suffix first, then prefix
    assert variant_key("1", 100, "a", "t") == ("chr1", 100, "A", "T")


def test_variant_kinds():
    assert variant_kind("A", "G") == "snv"
    assert variant_kind("", "T") == "ins" and variant_kind("C", "") == "del"
    assert variant_kind("AC", "GT") == "mnv"
    assert variant_kind("GGAT", "TGCTG") == "complex"


# ---------------------------------------------------------------- VCF reader
VCF_HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"


def write_vcf(path, records, gz=False):
    text = VCF_HEADER + "".join("\t".join(map(str, r)) + "\n" for r in records)
    if gz:
        with gzip.open(path, "wt") as f:
            f.write(text)
    else:
        Path(path).write_text(text)


def test_reader_splits_multiallelics_and_reads_fields(tmp_path):
    p = tmp_path / "s.vcf"
    write_vcf(p, [
        ("chr1", 100, ".", "A", "G", "50", "PASS", ".", "GT:GQ:DP:VAF", "0/1:48:30:0.47"),
        ("1", 200, ".", "C", "T", ".", "PASS", ".", "GT", "1/1"),
        ("chr1", 300, ".", "A", "AT,AG", "40", "PASS", ".", "GT", "1/2"),
        ("chr1", 400, ".", "G", "A", "10", "RefCall", ".", "GT", "0/1"),
        ("chr1", 500, ".", "A", "<DEL>", "30", "PASS", ".", "GT", "0/1"),
        ("chr1", 600, ".", "A", "G,*", "30", "PASS", ".", "GT", "0/2"),
        ("chr1", 700, ".", "T", "C", "30", "PASS", ".", "GT", "0/0"),
        ("chr1", 800, ".", "T", "C", "30", "PASS", ".", "GT", "./."),
    ])
    stats = vcf.ReadStats()
    calls = list(vcf.iter_calls(str(p), stats))
    assert [(c.chrom, c.pos, c.ref, c.alt, c.copies) for c in calls] == [
        ("chr1", 100, "A", "G", 1), ("chr1", 200, "C", "T", 2), ("chr1", 300, "A", "AT", 1), ("chr1", 300, "A", "AG", 1),
        ("chr1", 400, "G", "A", 1),
    ]
    first = calls[0]
    assert (first.dp, first.gq, first.vaf, first.qual, first.kind) == (30, 48, 0.47, 50.0, "snv")
    assert calls[2].kind == "ins" and (calls[2].norm_pos, calls[2].norm_ref, calls[2].norm_alt) == (301, "", "T")
    assert calls[4].filter == "RefCall"
    assert stats.records == 8 and stats.symbolic_or_star == 2 and stats.no_alt_call == 2
    assert stats.samples == 1 and stats.sample_name == "S1"


def test_reader_phase_haplotypes_and_gzip(tmp_path):
    p = tmp_path / "s.vcf.gz"
    write_vcf(p, [
        ("chr2", 10, ".", "A", "C", "50", "PASS", ".", "GT:PS", "0|1:1000"),
        ("chr2", 20, ".", "G", "T", "50", "PASS", ".", "GT:PS", "1|0:1000"),
        ("chr2", 30, ".", "G", "T", "50", "PASS", ".", "GT:PS", "1|1:1000"),
        ("chr2", 40, ".", "G", "T", "50", "PASS", ".", "GT", "0/1"),
    ], gz=True)
    calls = list(vcf.iter_calls(str(p)))
    assert [(c.phased, c.ps, c.haplotypes) for c in calls] == [
        (True, "1000", "2"), (True, "1000", "1"), (True, "1000", "12"), (False, None, ""),
    ]


# ---------------------------------------------------------------- gnotate reader
def test_packed_allele_layout_matches_the_real_resource():
    # the SNV table learned from the real file: ref and alt in 2 bits each, lengths in two nibbles
    assert encode_allele("A", "C") == 0x110000004
    assert encode_allele("C", "T") == 0x11000001C
    assert encode_allele("G", "A") == 0x110000020
    assert encode_allele("T", "G") == 0x110000038
    # a real deletion from gnomAD v3.0 chr1 (pos 10108 0-based): ref AACCCT, alt A
    assert encode_allele("AACCCT", "A") == 0x610000570
    assert encode_allele("A" * 7, "C" * 6) is not None  # 13 bases still fit
    assert encode_allele("A" * 7, "C" * 7) is None  # 14 do not
    assert encode_allele("A", "N") is None


def make_gnotate(path, fields, by_chrom, long_rows=None):
    """Write a gnotate zip. by_chrom: {chrom: [(pos_vcf, ref, alt, [values in field order])]}."""
    long_rows = long_rows or {}
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("sli.var/chroms.txt", "\n".join(by_chrom) + "\n")
        z.writestr("sli.var/fields.txt", "\n".join(fields) + "\n")
        z.writestr("sli.var/args.txt", "version:test\n")
        z.writestr("sli.var/message.txt", "created on:2020-01-01\n")
        for chrom, entries in by_chrom.items():
            keyed = sorted(((((pos - 1) << 36) | encode_allele(ref, alt)), vals) for pos, ref, alt, vals in entries)
            z.writestr(f"sli.var/{chrom}/gnotate-variant.bin", array("Q", [k for k, _ in keyed]).tobytes())
            for i, f in enumerate(fields):
                z.writestr(f"sli.var/{chrom}/gnotate-{f}.bin", array("f", [v[i] for _, v in keyed]).tobytes())
            rows = long_rows.get(chrom, [])
            z.writestr(
                f"sli.var/{chrom}/long-alleles.txt",
                "position\treference\talternate\tfilter\t" + "|".join(f.split("_")[-1] for f in fields) + "\n"
                + "".join(f"{p}\t{r}\t{a}\tfalse\t" + "|".join(map(str, v)) + "\n" for p, r, a, v in rows),
            )


GNOMAD_FIELDS = ["gnomad_af", "gnomad_nhomalt", "gnomad_ac"]


def test_gnotate_lookup_snv_indel_long_and_absent(tmp_path):
    long_ref, long_alt = "ACGTACGTAC", "ACGT"
    make_gnotate(
        tmp_path / "g.zip", GNOMAD_FIELDS,
        {"1": [(100, "A", "G", [0.3, 5, 40000]), (600, "TC", "T", [0.02, 1, 2800]), (700, "A", "AT", [0.0005, 0, 70])],
         "X": [(50, "C", "T", [0.1, 2, 9000])]},
        long_rows={"1": [(999, long_ref, long_alt, [0.004, 0, 280.0])]},  # 0-based position, like the arrays
    )
    g = Gnotate(str(tmp_path / "g.zip"))
    assert g.fields == GNOMAD_FIELDS and g.has_chrom("chr1") and g.has_chrom("X") and not g.has_chrom("chr2")
    found, vals = g.lookup(
        "chr1", [100, 600, 700, 100, 1000, 1000], ["A", "TC", "A", "A", long_ref, long_ref], ["G", "T", "AT", "T", long_alt, "A"]
    )
    assert found == [True, True, True, False, True, False]  # A>T absent; the long allele is found; another one is not
    assert vals["gnomad_af"][0] == pytest.approx(0.3, rel=1e-5) and vals["gnomad_ac"][1] == pytest.approx(2800)
    assert vals["gnomad_af"][4] == pytest.approx(0.004)  # from long-alleles.txt
    # 0-based positions: the same alleles one base to the right are not the variant
    assert g.lookup("chr1", [101], ["A"], ["G"])[0] == [False]
    assert g.lookup("chr2", [100], ["A"], ["G"])[0] == [False]  # chromosome not in the file


# ---------------------------------------------------------------- genes
def test_gene_overlap_including_a_gene_that_spans_others():
    idx = GeneIndex([("chr1", 100, 1000, "BIG"), ("chr1", 200, 300, "INNER"), ("chr1", 5000, 5100, "FAR"), ("1", 50, 60, "OTHER")])
    assert idx.genes_at("chr1", 250, 251) == ("BIG", "INNER")  # prefix-max keeps the long gene in play
    assert idx.genes_at("chr1", 400, 401) == ("BIG",)
    assert idx.genes_at("chr1", 1000, 1001) == ()  # half-open: 1000 is outside
    assert idx.genes_at("chr1", 55, 56) == ("OTHER",) and idx.genes_at("chr3", 1, 2) == ()


# ---------------------------------------------------------------- ClinVar
def test_clinvar_classification_groups_and_stars():
    assert classification_group("Pathogenic") == "plp"
    assert classification_group("Likely pathogenic") == "plp"
    assert classification_group("Pathogenic/Likely pathogenic") == "plp"
    assert classification_group("Pathogenic; risk factor") == "plp"
    assert classification_group("Likely pathogenic, low penetrance") == "plp"
    assert classification_group("Conflicting classifications of pathogenicity") == "conflicting"
    for other in ("Uncertain significance", "Benign", "Likely benign", "risk factor", "not provided"):
        assert classification_group(other) == "other"
    assert review_stars("practice guideline") == 4
    assert review_stars("reviewed by expert panel") == 3
    assert review_stars("criteria provided, multiple submitters, no conflicts") == 2
    assert review_stars("criteria provided, single submitter") == 1
    assert review_stars("criteria provided, conflicting classifications") == 1
    assert review_stars("no assertion criteria provided") == 0


def write_clinvar(path, rows):
    with gzip.open(path, "wt") as f:
        f.write("# release: 2099-01\n" + "\t".join(COLUMNS) + "\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")


CLINVAR_ROWS = [
    # chrom pos ref alt gene classification group stars review variation_id name conditions last_evaluated n_submitters
    ("chr1", 200, "C", "T", "GENEA", "Pathogenic", "plp", 2, "criteria provided, multiple submitters, no conflicts", "1111",
     "NM_000001.1(GENEA):c.10C>T (p.Arg4Ter)", "Disease A", "Jan 1, 2025", 3),
    ("chr1", 300, "G", "A", "GENEA", "Conflicting classifications of pathogenicity", "conflicting", 1,
     "criteria provided, conflicting classifications", "2222", "NM_000001.1(GENEA):c.20G>A", "Disease A; Disease B", "Feb 1, 2025", 2),
    # the insertion A>AG written the trimmed way a sample's A>AG normalises to
    ("chr1", 701, "", "G", "GENEB", "Likely pathogenic", "plp", 3, "reviewed by expert panel", "3333", "NM_000002.1(GENEB):c.5dup", "Disease B", "Mar 1, 2025", 5),
]


def test_clinvar_table_loads_release_and_looks_up_exact_alleles(tmp_path):
    write_clinvar(tmp_path / "cv.tsv.gz", CLINVAR_ROWS)
    t = ClinvarTable.load(tmp_path / "cv.tsv.gz")
    assert t.release == "2099-01" and len(t) == 3
    rec = t.lookup("chr1", 200, "C", "T")
    assert (rec.gene, rec.group, rec.stars, rec.variation_id, rec.n_submitters) == ("GENEA", "plp", 2, "1111", 3)
    assert t.lookup("chr1", 200, "C", "G") is None and t.lookup("chr1", 201, "C", "T") is None


# ---------------------------------------------------------------- the index
SAMPLE_RECORDS = [
    ("chr1", 100, ".", "A", "G", "50", "PASS", ".", "GT:GQ:DP:VAF", "0/1:50:30:0.5"),
    ("chr1", 200, ".", "C", "T", "60", "PASS", ".", "GT:GQ:DP:VAF", "1/1:60:28:1.0"),
    ("chr1", 300, ".", "G", "A", "55", "PASS", ".", "GT:PS", "0|1:300"),
    ("chr1", 600, ".", "TC", "T", "45", "PASS", ".", "GT:GQ:DP", "0/1:40:25"),
    ("chr1", 700, ".", "A", "AG,AT", "50", "PASS", ".", "GT:GQ", "1/2:45"),
    ("chr1", 800, ".", "G", "C", "10", "RefCall", ".", "GT:GQ", "0/1:5"),
    ("chr1", 900, ".", "A", "<DEL>", "30", "PASS", ".", "GT", "0/1"),
    ("chr1", 1000, ".", "ACGTACGTAC", "ACGT", "50", "PASS", ".", "GT", "0/1"),
    ("chr2", 10, ".", "A", "C", "50", "PASS", ".", "GT:GQ", "0/1:30"),
    ("chrX", 3000000, ".", "A", "T", "50", "PASS", ".", "GT:GQ", "1/1:50"),
    ("chrX", 1000000, ".", "A", "T", "50", "PASS", ".", "GT:GQ", "1/1:50"),  # pseudo-autosomal
]
GENE_ROWS = [("chr1", 50, 500, "GENEA"), ("chr1", 650, 800, "GENEB"), ("chrX", 2900000, 3100000, "GENEX")]


@pytest.fixture
def env(tmp_path, monkeypatch):
    """Isolated resources and cache; returns a helper to register a sample and build its index."""
    root = tmp_path / "samples"
    (root / "S1").mkdir(parents=True)
    write_vcf(root / "S1" / "S1.vcf", SAMPLE_RECORDS)
    (root / "S1" / "sample.json").write_text(json.dumps({"label": "Synthetic S1", "sex": "male", "sex_source": "test"}))
    bed = tmp_path / "genes.bed"
    bed.write_text("".join(f"{c}\t{s}\t{e}\t{n}\t+\n" for c, s, e, n in GENE_ROWS))
    write_clinvar(tmp_path / "cv.tsv.gz", CLINVAR_ROWS)
    make_gnotate(
        tmp_path / "gnomad.zip", GNOMAD_FIELDS,
        {"1": [(100, "A", "G", [0.30, 5, 40000]), (300, "G", "A", [0.00005, 0, 4]), (600, "TC", "T", [0.02, 1, 2800]), (700, "A", "AT", [0.5, 9, 99999])]},
        long_rows={"1": [(999, "ACGTACGTAC", "ACGT", [0.004, 0, 280.0])]},
    )
    make_gnotate(tmp_path / "hprc.zip", ["hprc_af", "hprc_nhomalt", "hprc_ac"], {"2": [(10, "A", "C", [0.2, 1, 20])]})
    monkeypatch.setattr(samples, "FIXTURE_DIR", tmp_path / "no_fixture")
    monkeypatch.setattr(samples, "RESULTS_ROOT", root)
    monkeypatch.setattr(samples, "CACHE_DIR", tmp_path / "cache")
    monkeypatch.setattr(samples, "GENES_BED", str(bed))
    monkeypatch.setattr(samples, "CLINVAR_PATH", str(tmp_path / "cv.tsv.gz"))
    monkeypatch.setattr(samples, "GNOMAD_PATH", str(tmp_path / "gnomad.zip"))
    monkeypatch.setattr(samples, "HPRC_PATH", str(tmp_path / "hprc.zip"))
    for registry in (samples._resources, samples._overview_cache, samples._threads, samples._progress, samples._errors):
        registry.clear()
    return tmp_path


def build(key="S1"):
    s = samples.resolve(key)
    assert samples.ensure_index(s)["state"] == "building"
    samples._threads[s.key].join(timeout=60)
    assert samples.status(s)["state"] == "ready", samples.status(s)
    return s


def run(s, **kw):
    f = samples.Filters(**{k: v for k, v in kw.items() if k in samples.Filters._fields})
    rest = {k: v for k, v in kw.items() if k not in samples.Filters._fields}
    return samples.query_variants(s, f, **rest)


def test_index_annotates_every_call(env):
    s = build()
    meta = index_mod.read_meta(samples.db_path(s))
    # 11 records: the symbolic <DEL> is skipped (-1) and the 1/2 site gives two calls (+1)
    assert meta["n_variants"] == "11"
    assert meta["clinvar_release"] == "2099-01" and meta["n_clinvar_hits"] == "3"
    assert meta["n_gnomad_hits"] == "5" and meta["n_hprc_hits"] == "1"  # incl. the long allele and A>AT
    r = run(s, page_size=100)
    by = {(x["chrom"], x["pos"], x["alt"]): x for x in r["rows"]}
    assert by[("chr1", 100, "G")]["gnomad"]["af"] == pytest.approx(0.30, rel=1e-5)
    assert by[("chr1", 100, "G")]["genes"] == ["GENEA"] and by[("chr1", 100, "G")]["clinvar"] is None
    assert by[("chr1", 200, "T")]["clinvar"]["classification"] == "Pathogenic" and by[("chr1", 200, "T")]["gnomad"] is None
    assert by[("chr1", 700, "AG")]["clinvar"]["stars"] == 3  # matched through normalisation (A>AG == ''>G at 701)
    assert by[("chr1", 700, "AT")]["gnomad"]["af"] == pytest.approx(0.5) and by[("chr1", 700, "AT")]["clinvar"] is None
    assert by[("chr1", 1000, "ACGT")]["gnomad"]["af"] == pytest.approx(0.004)  # long allele
    assert by[("chr2", 10, "C")]["hprc"]["af"] == pytest.approx(0.2) and by[("chr2", 10, "C")]["gnomad"] is None
    assert by[("chr1", 600, "T")]["genes"] == []  # intergenic


def test_qc_counts_and_sex(env):
    s = build()
    qc = samples.overview(s)["qc"]
    assert qc["records"] == 11 and qc["calls_passing"] == 10 and qc["calls_filtered"] == 1
    # SNVs: chr1:100/200/300, chr2:10 and two chrX; insertions: the 1/2 site; deletions: TC>T and the 6-base trim
    assert (qc["snv"], qc["insertions"], qc["deletions"]) == (6, 2, 2)
    assert qc["het"] + qc["hom"] == 10 and qc["hom"] == 3
    assert qc["symbolic_or_star_alleles"] == 1
    assert qc["inferred_sex"] == "undetermined"  # far too few chrX calls to infer
    assert samples.overview(s)["sex"] == {"sex": "male", "source": "test"}  # sample.json wins, with its source
    assert qc["phased_het_fraction"] == pytest.approx(1 / 7, abs=1e-3)  # one of the seven het calls is phased
    assert qc["ti_tv"] is not None


def test_sex_inference_from_the_call_set():
    assert index_mod.infer_sex(20, 4000, 60000, 40000)[0] == "male"  # X is almost all hom-alt
    assert index_mod.infer_sex(2400, 1600, 60000, 40000)[0] == "female"
    assert index_mod.infer_sex(1200, 2800, 60000, 40000)[0] == "undetermined"
    sex, basis = index_mod.infer_sex(1, 5, 60000, 40000)
    assert sex == "undetermined" and "too few" in basis


def test_filters_funnel_and_genome_wide_scope(env):
    s = build()
    # genome-wide: calls on chr1, chr2 and chrX all come back, not a gene list's worth
    assert {x["chrom"] for x in run(s, page_size=100)["rows"]} == {"chr1", "chr2", "chrX"}
    assert run(s)["total_matching"] == 10 and run(s, pass_only=False)["total_matching"] == 11
    assert run(s, clinvar="plp")["total_matching"] == 2
    assert run(s, clinvar="plp_conflicting")["total_matching"] == 3
    assert run(s, clinvar="no_match")["total_matching"] == 7
    assert run(s, clinvar="plp", min_stars=3)["total_matching"] == 1
    assert run(s, zygosity="hom")["total_matching"] == 3
    assert run(s, kinds=("ins",))["total_matching"] == 2
    assert run(s, chrom="chr2")["total_matching"] == 1
    assert run(s, in_gene=True)["total_matching"] == 6  # GENEA holds 3 calls, GENEB 2 and GENEX 1; the rest are intergenic
    assert run(s, genes=("genea",))["total_matching"] == 3  # case-insensitive; chr1:50-500 holds three calls
    assert run(s, genes=("GENEA", "GENEX"))["total_matching"] == 4
    # frequency: absent counts as rare by default, and can be excluded
    rare = run(s, max_af=0.01)
    assert rare["total_matching"] == 7  # 10 calls minus the 0.30, 0.02 and 0.5 ones
    assert run(s, max_af=0.01, af_missing="exclude")["total_matching"] == 2  # only the 5e-5 and 0.004 calls have a frequency
    assert run(s, max_hprc_af=0.1)["total_matching"] == 9
    # GQ is applied only where the VCF has one; calls without GQ are kept
    assert run(s, min_gq=45)["total_matching"] == 8  # GQ 40 and 30 drop out; the two with no GQ stay
    funnel = run(s, clinvar="plp", genes=("GENEB",), explain=True)["funnel"]
    assert [(f["label"], f["count"]) for f in funnel] == [
        ("all calls in the VCF", 11), ("passes the VCF FILTER", 10), ("in GENEB", 2), ("ClinVar pathogenic / likely pathogenic", 1),
    ]
    # an empty result always explains itself
    empty = run(s, clinvar="plp", genes=("GENEX",))
    assert empty["total_matching"] == 0 and empty["funnel"][-1]["count"] == 0 and empty["funnel"][-2]["count"] == 1


def test_sorting_and_paging(env):
    s = build()
    pos = [(x["chrom"], x["pos"]) for x in run(s, page_size=100)["rows"]]
    assert pos == sorted(pos, key=lambda c: (samples.index_mod.chrom_order(c[0]), c[1]))
    page2 = run(s, page_size=4, page=2)
    assert page2["page"] == 2 and len(page2["rows"]) == 4 and page2["total_matching"] == 10
    by_af = run(s, sort="af", order="asc", max_af=1.0)["rows"]
    assert by_af[0]["gnomad"]["af"] == pytest.approx(0.00005)  # rarest first
    assert [x["gnomad"] is None for x in by_af][-1] in (True, False) and by_af[-1]["gnomad"] is None  # unknown frequency last
    top = run(s, sort="stars", order="desc", clinvar="has_match")["rows"]
    assert [x["clinvar"]["stars"] for x in top] == [3, 2, 1]


def test_zygosity_labels_use_sample_sex_outside_the_par(env):
    s = build()
    rows = {x["pos"]: x for x in run(s, chrom="chrX", page_size=10)["rows"]}
    assert rows[3000000]["genotype"]["zygosity"] == "hemizygous"  # male, chrX non-PAR
    assert rows[1000000]["genotype"]["zygosity"] == "hom"  # PAR1: diploid in both sexes
    assert run(s, chrom="chr1")["rows"][0]["genotype"]["zygosity"] in ("het", "hom")


def test_overview_funnel_resources_and_evidence(env):
    s = build()
    o = samples.overview(s)
    labels = [f["label"] for f in o["funnel"]]
    assert labels == [
        "All calls in the VCF", "Pass the VCF filter", "Rare: gnomAD AF ≤ 1% or absent", "In a gene body",
        "ClinVar pathogenic / likely pathogenic",
    ]
    counts = [f["count"] for f in o["funnel"]]
    assert counts == sorted(counts, reverse=True) and counts[0] == 11 and counts[1] == 10
    assert o["clinvar"] == {"plp": 2, "conflicting": 1}
    assert {g["gene"] for g in o["top_plp_genes"]} == {"GENEA", "GENEB"}
    installed = {r["key"]: r["installed"] for r in o["resources"]}
    assert installed == {"clinvar": True, "gnomad": True, "hprc": True, "genes": True, "consequence": False, "insilico": False, "hpo": False, "parents": False}
    ev = {e["criteria"]: e for e in o["evidence"]}
    assert ev["BA1 · BS1 · PM2"]["available"] and not ev["PVS1 · PM4 · BP3"]["available"] and not ev["PP4"]["available"]
    assert "not evaluated" in next(r["detail"] for r in o["resources"] if r["key"] == "consequence")


def test_missing_resources_read_as_not_evaluated(env, monkeypatch):
    monkeypatch.setattr(samples, "CLINVAR_PATH", str(env / "missing_cv.tsv.gz"))
    monkeypatch.setattr(samples, "GNOMAD_PATH", str(env / "missing_gnomad.zip"))
    s = build()
    r = run(s, clinvar="plp", min_stars=2, max_af=0.01)
    assert set(r["ignored_filters"]) == {"clinvar", "min_stars", "max_af"}  # the filters cannot act, and the response says so
    assert not r["gnomad_installed"] and not r["clinvar_installed"]
    assert r["total_matching"] == 10  # nothing was silently filtered out
    assert all(x["clinvar"] is None and x["gnomad"] is None for x in r["rows"])
    o = samples.overview(s)
    assert o["clinvar"] is None and o["gnomad_coverage"] is None
    assert [f["label"] for f in o["funnel"]] == ["All calls in the VCF", "Pass the VCF filter", "In a gene body"]


def test_index_rebuilds_when_the_vcf_or_a_resource_changes(env):
    s = build()
    before = index_mod.read_meta(samples.db_path(s))["built"]
    assert samples.status(s)["state"] == "ready"
    write_vcf(Path(s.vcf), SAMPLE_RECORDS[:3])  # the VCF changes
    assert samples.status(s)["state"] == "stale"
    s = build()
    assert index_mod.read_meta(samples.db_path(s))["n_variants"] == "3"
    write_clinvar(env / "cv.tsv.gz", CLINVAR_ROWS[:1])  # a resource changes
    assert samples.status(s)["state"] == "stale"


def test_a_failed_build_is_reported(env, monkeypatch):
    s = samples.resolve("S1")
    monkeypatch.setattr(index_mod, "build_index", lambda *a, **k: (_ for _ in ()).throw(RuntimeError("boom")))
    samples.ensure_index(s)
    samples._threads[s.key].join(timeout=30)
    st = samples.status(s)
    assert st["state"] == "error" and "boom" in st["error"]


# ---------------------------------------------------------------- API
def test_api_flow_status_variants_detail_and_errors(env):
    with TestClient(app) as client:
        assert client.get("/api/smallvar/samples").json() == {"samples": [{"key": "S1", "label": "Synthetic S1"}], "default": "S1"}
        assert client.get("/api/smallvar/variants", params={"sample": "nope"}).status_code == 404
        status = client.get("/api/smallvar/status", params={"sample": "S1"}).json()
        assert status["state"] in ("building", "ready")
        samples._threads["S1"].join(timeout=60)
        assert client.get("/api/smallvar/status").json() == {"state": "ready"}
        meta = client.get("/api/smallvar/meta").json()
        assert meta["label"] == "Synthetic S1" and meta["n_variants"] == 11 and meta["sex"]["sex"] == "male"
        r = client.get("/api/smallvar/variants", params={"clinvar": "plp", "max_af": 0.01, "genes": "GENEA,GENEB"}).json()
        assert r["total_matching"] == 2 and {x["alt"] for x in r["rows"]} == {"T", "AG"}
        one = client.get(f"/api/smallvar/variants/{r['rows'][0]['id']}").json()
        assert one["clinvar"]["name"].startswith("NM_")
        assert client.get("/api/smallvar/variants/999999").status_code == 404
        assert client.get("/api/smallvar/variants", params={"zygosity": "weird"}).status_code == 422
        assert client.get("/api/smallvar/variants", params={"clinvar": "weird"}).status_code == 422
        assert client.get("/api/smallvar/overview").json()["clinvar"]["plp"] == 2
        presets = client.get("/api/smallvar/gene-presets").json()["presets"]
        assert presets[0]["key"] == "acmg_sf_v3_3" and len(presets[0]["genes"]) == 84


def test_api_data_routes_answer_503_until_the_index_is_built(env, monkeypatch):
    release = __import__("threading").Event()
    real = index_mod.build_index

    def slow(*a, **k):
        release.wait(timeout=30)
        return real(*a, **k)

    monkeypatch.setattr(index_mod, "build_index", slow)
    with TestClient(app) as client:
        assert client.get("/api/smallvar/status").json()["state"] == "building"
        resp = client.get("/api/smallvar/variants")
        assert resp.status_code == 503 and resp.json()["detail"]["state"] == "building"
        release.set()
        samples._threads["S1"].join(timeout=60)
        assert client.get("/api/smallvar/variants").status_code == 200


# ---------------------------------------------------------------- reference data and the committed fixture
def test_acmg_sf_reference_table_matches_the_spec():
    with open(ACMG_TSV, encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert len(rows) == 84 and len({r["gene"] for r in rows}) == 84
    count = lambda key: {v: sum(1 for r in rows if r[key] == v) for v in {r[key] for r in rows}}
    assert count("category") == {"Cancer predisposition": 28, "Cardiovascular": 41, "Inborn errors of metabolism": 7, "Other": 8}
    assert count("inheritance") == {"AD": 71, "AR": 9, "XL": 3, "SD": 1}
    by = {r["gene"]: r for r in rows}
    assert by["HFE"]["rule_ids"] == "R6" and by["TTN"]["rule_ids"] == "R2, R7" and by["ABCD1"]["inheritance"] == "XL"
    assert {g for g, r in by.items() if r["added_in"] == "v3.3"} == {"ABCD1", "CYP27A1", "PLN"}


def test_committed_hg002_chr22_fixture_builds_and_reconciles(tmp_path, monkeypatch):
    """Public GIAB data: 50,284 records, of which 577 are multi-allelic (1/2 or 2/1) and give two calls."""
    assert FIXTURE_VCF.is_file()
    monkeypatch.setattr(samples, "CACHE_DIR", tmp_path / "cache")
    monkeypatch.setattr(samples, "CLINVAR_PATH", str(tmp_path / "none"))
    monkeypatch.setattr(samples, "GNOMAD_PATH", str(tmp_path / "none"))
    monkeypatch.setattr(samples, "HPRC_PATH", str(tmp_path / "none"))
    monkeypatch.setattr(samples, "RESULTS_ROOT", tmp_path / "none")
    for registry in (samples._resources, samples._overview_cache, samples._threads, samples._progress, samples._errors):
        registry.clear()
    assert [s.key for s in samples.discover()] == ["HG002_chr22"]
    s = build("HG002_chr22")
    qc = samples.overview(s)["qc"]
    assert qc["records"] == 50284 and qc["calls_passing"] == 50861
    assert (qc["het"], qc["hom"]) == (33520, 17341)
    assert qc["snv"] + qc["insertions"] + qc["deletions"] + qc["complex"] + qc["mnv"] == 50861
    assert samples.overview(s)["sex"]["sex"] == "male"
