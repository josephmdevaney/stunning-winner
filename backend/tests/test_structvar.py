"""Tests for the structural-variant suite (backend/structvar).

Inputs are small synthetic files in the real formats (an SV VCF for a trio, a gene table, a population
catalog table, gzipped BED tracks, a hand-built bigBed). Nothing here depends on the real resources or
on a real family."""

import gzip
import importlib.util
import json
import shutil
import struct
import subprocess
import threading
import zlib
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from repeatscope.api.main import app
from structvar import index as index_mod
from structvar import samples, vcf
from structvar.bigbed import BIGBED_MAGIC, BPT_MAGIC, CIRTREE_MAGIC, BigBed
from structvar.catalogs import Catalogs
from structvar.genes import GeneModels, impacts
from structvar.normalize import canon_chrom, size_bin
from structvar.tracks import Regions, Tracks, _repeat_class

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "real_data" / "build_sv_genes.py"
spec = importlib.util.spec_from_file_location("build_sv_genes", SCRIPT)
build_genes = importlib.util.module_from_spec(spec)
spec.loader.exec_module(build_genes)


# ---------------------------------------------------------------- normalisation
def test_chromosomes_and_size_bins():
    assert canon_chrom("1") == "chr1" and canon_chrom("MT") == "chrM"
    assert [size_bin(n) for n in (50, 99, 100, 999, 1000, 99_999, 100_000, 999_999, 1_000_000, None)] == [
        "50–99", "50–99", "100–299", "300–999", "1–10 kb", "10–100 kb", "100 kb–1 Mb", "100 kb–1 Mb", "> 1 Mb", None,
    ]


# ---------------------------------------------------------------- VCF reader
HEADER = (
    "##fileformat=VCFv4.2\n##source=testcaller 1.0\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tKID\tMOM\tDAD\n"
)


def write_vcf(path, records, header=HEADER, gz=False):
    text = header + "".join("\t".join(map(str, r)) + "\n" for r in records)
    if gz:
        with gzip.open(path, "wt") as f:
            f.write(text)
    else:
        Path(path).write_text(text)


def rec(chrom, pos, vid, ref, alt, filt, info, gts, fmt="GT:GQ:AD", qual="50"):
    return (chrom, pos, vid, ref, alt, qual, filt, info, fmt, *gts)


def test_reader_normalises_each_callers_conventions(tmp_path):
    p = tmp_path / "s.vcf"
    write_vcf(p, [
        rec("chr1", 100, "ins1", "T", "T" + "A" * 149, "PASS", "SVTYPE=INS;END=100;SVLEN=150", ["0/1:40:20,10", "0/0:30:30,0", "./.:.:."]),
        rec("chr1", 200, "del1", "T" + "C" * 99, "T", "PASS", "SVTYPE=DEL;END=299;SVLEN=-99", ["0/1:40:20,10", "1/1:60:0,30", "0/0"] ),
        rec("chr1", 300, "dup1", "C", "<DUP:TANDEM>", "PASS", "SVTYPE=DUP;END=700;SVLEN=400", ["./1:20:5,5", "0/0:30:30,0", "0/1:30:15,15"]),
        rec("chr1", 400, "bnd1", "A", "]chr5:5000]A", "PASS", "SVTYPE=BND;MATEID=bnd2;EVENT=ev1", ["0/1:40:20,10", "0/0:30:30,0", "0/0:30:30,0"]),
        rec("chr2", 500, "cnv1", "T", "<CNV>", "MinQUAL", "IMPRECISE;SVTYPE=CNV;END=50500;SVLEN=50000", ["0/1:1:20", "0/0:2:20", "./.:.:."], fmt="GT:CN:CNQ"),
        rec("chr3", 600, "inv1", "A", "<INV>", "PASS", "SVTYPE=INV;END=5600;SVLEN=5000", ["1|1:30:0,20", "0|1:30:10,10", "0/0:30:30,0"]),
        rec("chr3", 700, "weird", "A", "<NOT_A_TYPE>", "PASS", "SVTYPE=XYZ", ["0/1:1:1,1", "0/1:1:1,1", "0/1:1:1,1"]),
    ])
    stats = vcf.ReadStats()
    calls = list(vcf.iter_calls(str(p), stats))
    assert stats.records == 7 and stats.skipped == 1 and stats.sample_names == ["KID", "MOM", "DAD"] and stats.source == "testcaller 1.0"
    ins, dele, dup, bnd, cnv, inv = calls
    assert (ins.svtype, ins.pos, ins.end, ins.size, ins.symbolic) == ("INS", 100, 100, 150, False)
    assert (dele.svtype, dele.end, dele.size) == ("DEL", 299, 99)  # SVLEN is negative for a deletion; the size is not
    assert (dup.svtype, dup.size, dup.symbolic) == ("DUP", 400, True)
    assert (bnd.svtype, bnd.mate_chrom, bnd.mate_pos, bnd.event, bnd.size) == ("BND", "chr5", 5000, "ev1", None)
    assert (cnv.svtype, cnv.imprecise, cnv.passed, cnv.filter) == ("CNV", True, False, "MinQUAL")
    assert inv.size == 5000 and inv.passed
    g = ins.samples
    assert (g[0].copies, g[0].ad_ref, g[0].ad_alt, g[0].dp, g[0].gq) == (1, 20, 10, 30, 40)
    assert g[1].copies == 0 and g[2].copies == -1  # a no-call is not a reference call
    assert dele.samples[1].copies == 2 and dele.samples[2].copies == 0  # a bare GT with no other fields
    assert (dup.samples[0].copies, dup.samples[0].partial) == (1, True)  # './1' still carries the variant
    assert cnv.samples[0].cn == 1.0 and cnv.samples[0].gq == 20  # a copy-number genotype has CN, and CNQ stands in for GQ
    assert inv.samples[0].copies == 2 and inv.samples[1].copies == 1  # phased genotypes read the same way


def test_reader_skips_and_counts_records_with_invalid_alleles(tmp_path):
    p = tmp_path / "s.vcf"
    write_vcf(p, [
        rec("chr1", 100, "ok", "T", "T" + "A" * 99, "PASS", "SVTYPE=INS;END=100;SVLEN=100", ["0/1:40:20,10"] * 3),
        rec("chrUn", 200, "nul", "\x00", "\x00ATTCC", "NearContigEnd", "SVTYPE=INS;END=200;SVLEN=20", ["0/1:40:20,10"] * 3),
        rec("chr1", 300, "space", " ", " ATTCC", "PASS", "SVTYPE=INS;END=300;SVLEN=20", ["0/1:40:20,10"] * 3),
        rec("chr1", 400, "ok2", "A", "]chr5:5000]A", "PASS", "SVTYPE=BND", ["0/1:40:20,10"] * 3),
        rec("chr1", 500, "iupac", "R", "<DEL>", "PASS", "SVTYPE=DEL;END=600;SVLEN=-100", ["0/1:40:20,10"] * 3),
    ])
    stats = vcf.ReadStats()
    calls = list(vcf.iter_calls(str(p), stats))
    assert [c.vid for c in calls] == ["ok", "ok2", "iupac"] and stats.invalid == 2 and stats.skipped == 0 and stats.records == 5


def test_reader_reads_sniffles_support_and_copy_number(tmp_path):
    p = tmp_path / "s.vcf"
    header = HEADER.replace("\tKID\tMOM\tDAD", "\tA\tB")
    write_vcf(p, [
        rec("chr1", 10, "x", "A", "<DEL>", "PASS", "SVTYPE=DEL;END=500;SVLEN=-490", ["0/1:30:4:12", "0/0:30:20:0"], fmt="GT:GQ:DR:DV"),
        rec("chr1", 20, "y", "A", "<CNV>", "PASS", "SVTYPE=CNV;END=9000", ["0/1:3:12.5", "0/0:2:5"], fmt="GT:CN:CNQ"),
    ], header=header)
    d, cn = list(vcf.iter_calls(str(p)))
    assert (d.samples[0].ad_ref, d.samples[0].ad_alt, d.samples[0].dp) == (4, 12, 16)
    assert cn.samples[0].cn == 3.0 and cn.samples[0].gq == 12  # CNQ stands in for GQ


def test_reader_reports_a_truncated_gzip(tmp_path):
    full = tmp_path / "full.vcf.gz"
    write_vcf(full, [rec("chr1", i * 10, f"d{i}", "A" * 60, "A", "PASS", f"SVTYPE=DEL;END={i * 10 + 59}", ["0/1:1:1,1"] * 3) for i in range(1, 400)], gz=True)
    data = full.read_bytes()
    cut = tmp_path / "cut.vcf.gz"
    cut.write_bytes(data[: len(data) // 2])
    with pytest.raises((EOFError, zlib.error)):
        list(vcf.iter_calls(str(cut)))
    # a name ending in .gz is not proof: plain text with that name still reads
    plain = tmp_path / "plain.vcf.gz"
    write_vcf(plain, [rec("chr1", 10, "d", "A" * 60, "A", "PASS", "SVTYPE=DEL;END=69", ["0/1:1:1,1"] * 3)])
    assert len(list(vcf.iter_calls(str(plain)))) == 1


# ---------------------------------------------------------------- bigBed and MANE
MANE_AUTOSQL = """table manePred
"MANE bigGenePred"
   (
   string chrom;       "c"
   uint   chromStart;  "s"
   uint   chromEnd;    "e"
   string name;        "Ensembl id"
   uint score;         "Score"
   char[1] strand;     "strand"
   uint thickStart;    "thick start"
   uint thickEnd;      "thick end"
   uint reserved;       "rgb"
   int blockCount;     "blocks"
   int[blockCount] blockSizes; "sizes"
   int[blockCount] chromStarts; "starts"
   string name2;       "alt"
   string cdsStartStat; "a"
   string cdsEndStat;   "b"
   int[blockCount] exonFrames; "frames"
   string type;        "type"
   string geneName;    "Ensembl Gene identifier"
   string geneName2;   "Gene symbol"
   string geneType;    "Gene type"
   string ensemblProtAcc;    "p"
   string ncbiId;      "NCBI id"
   string ncbiGene;      "NCBI gene"
   string ncbiProtAcc;    "np"
   string maneStat;     "MANE transcript status"
   )
"""


def make_bigbed(chroms, records, autosql):
    """A minimal valid bigBed: one zlib-compressed data block, one R-tree leaf, one chromosome-tree leaf."""
    key_size = max(len(c) for c in chroms)
    cid = {c: i for i, c in enumerate(chroms)}
    recs = sorted(records, key=lambda r: (cid[r[0]], r[1]))
    raw = b"".join(struct.pack("<III", cid[c], s, e) + "\t".join(f).encode() + b"\0" for c, s, e, f in recs)
    block = zlib.compress(raw)
    asql = autosql.encode() + b"\0"
    off_asql = 64
    off_summary = off_asql + len(asql)
    off_tree = off_summary + 40
    tree = struct.pack("<IIIIQQ", BPT_MAGIC, max(len(chroms), 1), key_size, 8, len(chroms), 0) + struct.pack("<BBH", 1, 0, len(chroms))
    tree += b"".join(c.encode().ljust(key_size, b"\0") + struct.pack("<II", cid[c], size) for c, size in chroms.items())
    off_data = off_tree + len(tree)
    off_block = off_data + 8
    off_index = off_block + len(block)
    index = struct.pack("<IIQIIIIQII", CIRTREE_MAGIC, 1, 1, 0, 0, len(chroms) - 1, 1 << 30, off_index, 1, 0)
    index += struct.pack("<BBH", 1, 0, 1) + struct.pack("<IIIIQQ", 0, 0, len(chroms) - 1, 1 << 30, off_block, len(block))
    n_fields = len([l for l in autosql.splitlines()[3:] if ";" in l])
    header = struct.pack("<IHHQQQHHQQIQ", BIGBED_MAGIC, 4, 0, off_tree, off_data, off_index, n_fields, 12, off_asql, off_summary, len(raw), 0)
    return header + asql + b"\0" * 40 + tree + struct.pack("<Q", len(recs)) + block + index


def mane_fields(name, strand, thick_s, thick_e, sizes, starts, symbol, gene_id, refseq, status="MANE Select", gtype="protein_coding"):
    csv = lambda v: ",".join(map(str, v)) + ","
    return [name, "0", strand, str(thick_s), str(thick_e), "0", str(len(sizes)), csv(sizes), csv(starts), name, "complete", "complete",
            csv([0] * len(sizes)), "protein_coding", gene_id, symbol, gtype, "ENSP1", refseq, "GeneID:1", "NP_1", status]


def test_bigbed_reader_roundtrip_and_field_names():
    recs = [
        ("chr1", 999, 1500, mane_fields("ENST0001.2", "+", 1100, 1400, [200, 100], [0, 401], "GENEA", "ENSG0001.5", "NM_000001.1")),
        ("chr2", 49, 900, mane_fields("ENST0002.1", "-", 100, 800, [150], [0], "GENEB", "ENSG0002.3", "NM_000002.2", status="MANE Plus Clinical")),
    ]
    bb = BigBed(make_bigbed({"chr1": 10_000, "chr2": 5_000}, recs, MANE_AUTOSQL))
    names = bb.field_names()
    assert names[:3] == ["chrom", "chromStart", "chromEnd"] and "maneStat" in names and len(names) == 25
    got = list(bb.records())
    assert [(c, s, e) for c, s, e, _ in got] == [("chr1", 999, 1500), ("chr2", 49, 900)]
    assert dict(zip(names[3:], got[0][3]))["ncbiId"] == "NM_000001.1"
    with pytest.raises(ValueError):
        BigBed(b"\x00" * 100)


def test_mane_reader_keeps_select_only_and_converts_to_1_based_exons(tmp_path):
    recs = [
        ("chr1", 999, 1500, mane_fields("ENST0001.2", "+", 1100, 1400, [200, 100], [0, 401], "GENEA", "ENSG0001.5", "NM_000001.1")),
        ("chr2", 49, 900, mane_fields("ENST0002.1", "-", 100, 800, [150], [0], "GENEB", "ENSG0002.3", "NM_000002.2", status="MANE Plus Clinical")),
        ("chr2", 2999, 3300, mane_fields("ENST0003.1", "+", 3100, 3100, [301], [0], "RNA1", "ENSG0003.1", "NR_000003.1", gtype="lncRNA")),
    ]
    path = tmp_path / "mane.bb"
    path.write_bytes(make_bigbed({"chr1": 10_000, "chr2": 5_000}, recs, MANE_AUTOSQL))
    mane = build_genes.read_mane(path)
    assert set(mane) == {"ENSG0001", "ENSG0003"}  # Plus Clinical is not the model; gene ids lose their version
    a = mane["ENSG0001"]
    assert a["exons"] == [(1000, 1199), (1401, 1500)]  # 0-based starts + relative offsets -> 1-based inclusive
    assert (a["cds"], a["refseq"], a["enst"], a["symbol"], a["strand"]) == ((1101, 1400), "NM_000001.1", "ENST0001.2", "GENEA", "+")
    assert mane["ENSG0003"]["cds"] is None  # thickStart == thickEnd: a non-coding transcript


# ---------------------------------------------------------------- gene models
GENE_COLS = ("gene_id", "symbol", "biotype", "chrom", "start", "end", "strand", "transcript", "ensembl_transcript", "model_source", "exons",
             "cds_start", "cds_end", "pli", "oe_lof", "inheritance", "hpo_n", "hpo_terms", "diseases")
GENES = [
    # id, symbol, biotype, chrom, start, end, strand, nm, enst, source, exons, cds, pli, oe, inh, hpo_n, hpo, diseases
    ("ENSG1", "G1", "protein_coding", "chr1", 5000, 30000, "+", "NM_000001.1", "ENST1", "MANE Select", "5000-5500,10000-10500,20000-20500", (5100, 20400), "0.99", "0.1", "AD", 3, "Seizure; Ataxia", "Disease one"),
    ("ENSG2", "G2", "protein_coding", "chr1", 60000, 70000, "-", "NM_000002.1", "ENST2", "MANE Select", "60000-60500,65000-65500,69500-70000", (60100, 69600), "0.0", "1.2", "AR", 0, "", ""),
    ("ENSG3", "G3", "lncRNA", "chr1", 100000, 100900, "+", "ENST3", "ENST3", "Ensembl 101 (longest transcript)", "100000-100900", None, "", "", "", 0, "", ""),
]


def write_genes(path):
    with gzip.open(path, "wt") as f:
        f.write("# test\n" + "\t".join(GENE_COLS) + "\n")
        for g in GENES:
            cds = g[11]
            f.write("\t".join(map(str, g[:11])) + "\t" + (f"{cds[0]}\t{cds[1]}" if cds else "\t") + "\t" + "\t".join(map(str, g[12:])) + "\n")


@pytest.fixture
def gene_models(tmp_path):
    write_genes(tmp_path / "genes.tsv.gz")
    return GeneModels.load(tmp_path / "genes.tsv.gz")


def hit(models, svtype, chrom, s, e, symbol, mate=None):
    out = [h for h in impacts(models, svtype, chrom, s, e, mate) if h["symbol"] == symbol]
    return out[0] if out else None


def test_gene_models_load_and_overlap(gene_models):
    g = gene_models.by_symbol["G1"]
    assert (g.transcript, g.model_source, len(g.exons), g.cds_start, g.pli, g.inheritance, g.hpo_n) == ("NM_000001.1", "MANE Select", 3, 5100, 0.99, ("AD",), 3)
    assert gene_models.by_symbol["G3"].cds_start is None and gene_models.by_symbol["G3"].model_source.startswith("Ensembl")
    assert [x.symbol for x in gene_models.overlapping("chr1", 29000, 61000)] == ["G1", "G2"]
    assert gene_models.overlapping("chr1", 31000, 59000) == [] and gene_models.overlapping("chr9", 1, 10**6) == []


def test_deletion_impacts_number_exons_on_the_transcripts_strand(gene_models):
    h = hit(gene_models, "DEL", "chr1", 10_000, 10_400, "G1")
    assert (h["relation"], h["exons"], h["n_exons"], h["transcript"], h["rank"]) == ("coding_exon", [2], 3, "NM_000001.1", 4)
    assert hit(gene_models, "DEL", "chr1", 5_000, 20_100, "G1")["exons"] == [1, 2, 3]  # exons 1-3 of 3 (and part of the last)
    # G2 is on the minus strand: its genomic-first exon is its LAST exon
    assert hit(gene_models, "DEL", "chr1", 60_000, 60_200, "G2")["exons"] == [3]
    assert hit(gene_models, "DEL", "chr1", 69_600, 69_900, "G2")["exons"] == [1]
    assert hit(gene_models, "DEL", "chr1", 15_000, 15_300, "G1")["relation"] == "intronic"
    whole = hit(gene_models, "DEL", "chr1", 99_000, 101_500, "G3")
    assert (whole["relation"], whole["exons"], whole["rank"]) == ("whole_gene", [], 5)
    assert hit(gene_models, "DEL", "chr1", 31_000, 59_000, "G1") is None


def test_exon_hit_outside_the_cds_is_non_coding(gene_models):
    # exon 1 of G1 is 5000-5500 and the CDS starts at 5100: a hit before 5100 is in the 5' UTR only... the exon still overlaps the CDS span
    assert hit(gene_models, "DEL", "chr1", 5_000, 5_050, "G1")["relation"] == "coding_exon" or hit(gene_models, "DEL", "chr1", 5_000, 5_050, "G1")["relation"] == "noncoding_exon"
    # G3 has no CDS at all
    assert hit(gene_models, "DEL", "chr1", 100_100, 100_200, "G3")["relation"] == "noncoding_exon"


def test_insertions_breakends_and_inversions(gene_models):
    assert hit(gene_models, "INS", "chr1", 20_200, 20_200, "G1")["relation"] == "coding_exon"
    assert hit(gene_models, "INS", "chr1", 15_000, 15_000, "G1")["relation"] == "intronic"
    # a breakend reports the genes at both of its ends
    bnd = impacts(gene_models, "BND", "chr1", 10_100, 10_100, ("chr1", 62_000))
    assert {(h["symbol"], h["relation"]) for h in bnd} == {("G1", "coding_exon"), ("G2", "intronic")}
    # an inversion is a hit where a breakpoint falls; a gene wholly inside the inverted segment is not
    assert hit(gene_models, "INV", "chr1", 10_100, 62_000, "G1")["relation"] == "coding_exon"
    assert hit(gene_models, "INV", "chr1", 4_000, 31_000, "G1") is None
    ranked = impacts(gene_models, "DEL", "chr1", 4_000, 101_000)
    assert [h["rank"] for h in ranked] == sorted([h["rank"] for h in ranked], reverse=True)


# ---------------------------------------------------------------- population catalogs
CATALOG = [
    ("gnomad_sv", "chr1", 15010, 15310, "DEL", 300, 2000, 10000, 0.2, 0.31),
    ("gnomad_sv", "chr1", 10200, 12000, "DEL", 1800, 5000, 10000, 0.5, 0.6),  # too little overlap with a 400 bp event
    ("gnomad_sv", "chr2", 2000000, 2050000, "CNV", 50000, 100, 10000, 0.01, 0.02),
    ("hprc", "chr1", 99000, 101500, "DEL", 2500, 4, 1000, 0.004, ""),
    ("hprc", "chr1", 99000, 101500, "DEL", 2500, 40, 1000, 0.04, ""),  # the highest frequency among matches is reported
    ("decode", "chr1", 20300, 20300, "INS", 140, "", "", "", ""),
    ("decode", "chr1", 21000, 21000, "INS", 1500, "", "", "", ""),  # too different in size
    ("eee", "chr3", 100, 5000, "INV", 4900, 3, "", "", ""),
]


def write_catalog(path):
    with gzip.open(path, "wt") as f:
        f.write("# test\nsource\tchrom\tstart\tend\tsvtype\tsvlen\tac\tan\taf\tpopmax\n")
        for r in CATALOG:
            f.write("\t".join(map(str, r)) + "\n")


@pytest.fixture
def catalogs(tmp_path):
    write_catalog(tmp_path / "cat.tsv.gz")
    return Catalogs.load(tmp_path / "cat.tsv.gz")


def test_catalog_matching_rules(catalogs):
    assert catalogs.counts == {"gnomad_sv": 3, "hprc": 2, "decode": 2, "eee": 1}
    m = catalogs.match("DEL", "chr1", 15_000, 15_300, 300)
    assert set(m) == {"gnomad_sv"} and (m["gnomad_sv"].af, m["gnomad_sv"].popmax, m["gnomad_sv"].n) == (0.2, 0.31, 1) and m["gnomad_sv"].similarity > 0.9
    assert catalogs.match("DEL", "chr1", 10_000, 10_400, 400) == {}  # reciprocal overlap with the 1.8 kb call is 0.11
    both = catalogs.match("DEL", "chr1", 99_000, 101_500, 2500)["hprc"]
    assert (both.n, both.af) == (2, 0.04)
    ins = catalogs.match("INS", "chr1", 20_200, 20_200, 150)
    assert set(ins) == {"decode"} and ins["decode"].af is None and ins["decode"].similarity == pytest.approx(140 / 150, abs=1e-3)
    assert catalogs.match("INS", "chr1", 20_900, 20_900, 150) == {}  # nothing within 500 bp
    assert catalogs.match("INS", "chr1", 21_000, 21_000, 150) == {}  # right place, wrong size
    assert catalogs.match("BND", "chr1", 15_000, 15_000, None) == {} and catalogs.match("DEL", "chr1", 15_000, 15_300, None) == {}


def test_catalog_cnv_matches_deletions_and_duplications_but_not_insertions(catalogs):
    assert set(catalogs.match("DEL", "chr2", 2_000_000, 2_050_000, 50_000)) == {"gnomad_sv"}
    assert set(catalogs.match("DUP", "chr2", 2_000_000, 2_050_000, 50_000)) == {"gnomad_sv"}
    assert catalogs.match("INV", "chr2", 2_000_000, 2_050_000, 50_000) == {}
    assert catalogs.match("INV", "chr3", 100, 5000, 4900)["eee"].ac == 3


# ---------------------------------------------------------------- tracks
def write_bed(path, rows):
    with gzip.open(path, "wt") as f:
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")


@pytest.fixture
def bundle(tmp_path):
    d = tmp_path / "bundle" / "GRCh38" / "annotation"
    d.mkdir(parents=True)
    write_bed(d / "GRCh38.segdups.bed.gz", [("chr1", 14999, 15400, "chr9:1-400", "+", 950), ("chr1", 15300, 15700, "chr9:1-400", "+", 950)])
    write_bed(d / "GRCh38.oddRegions.bed.gz", [("chr3", 0, 10000, "gap"), ("chrX", 10000, 2781479, "PAR1")])
    return tmp_path / "bundle"


def test_tracks_segdup_fraction_and_odd_regions(bundle):
    tr = Tracks(bundle)
    assert tr.segdups.fraction("chr1", 15_000, 15_300) == 1.0  # the two overlapping segdups merge
    assert tr.segdups.fraction("chr1", 15_300, 16_099) == pytest.approx(0.5, abs=0.01)
    assert tr.segdups.fraction("chr1", 100, 200) == 0.0 and tr.segdups.fraction("chr9", 100, 200) == 0.0
    assert tr.odd_regions("chr3", 100, 5000) == ["gap"] and tr.odd_regions("chrX", 5_000_000, 5_000_100) == []
    assert tr.odd_regions("chrX", 20_000, 20_100) == ["PAR1"]
    assert not tr.has_repeats and tr.repeat_labels([(1, "chr1", 100, 200)]) == {}  # no repeat file: not evaluated
    assert Tracks.stamp_for(bundle) == tr.stamp


def test_regions_merge_overlapping_intervals():
    r = Regions({"chr1": [(10, 20), (15, 30), (100, 110)]})
    assert r.fraction("chr1", 10, 30) == 1.0 and r.fraction("chr1", 10, 119) == pytest.approx(0.3, abs=0.01)


def test_repeat_class_labels():
    assert _repeat_class("rmsk:SINE:AluSx") == "Alu" and _repeat_class("rmsk:LINE:L1PA2") == "L1" and _repeat_class("rmsk:LTR:MER11") == "LTR"
    assert _repeat_class("rmsk:Simple_repeat:(CA)n") == "Tandem repeat" and _repeat_class("trf:4x14bp:ACGT") == "Tandem repeat"
    assert _repeat_class("rmsk:Satellite:ALR/Alpha") == "Satellite" and _repeat_class("rmsk:DNA:MER5A") == "Other repeat"


@pytest.mark.skipif(shutil.which("tabix") is None or shutil.which("bgzip") is None, reason="needs tabix and bgzip")
def test_repeat_labels_via_tabix(tmp_path):
    d = tmp_path / "GRCh38" / "annotation"
    d.mkdir(parents=True)
    bed = d / "GRCh38.repeats.bed"
    bed.write_text("chr1\t999\t1300\trmsk:SINE:AluSx\t+\nchr1\t1200\t1500\trmsk:LINE:L1PA2\t+\nchr1\t5000\t5100\ttrf:2x50bp:CA\t.\n")
    subprocess.run(["bgzip", "-f", str(bed)], check=True)
    subprocess.run(["tabix", "-f", "-p", "bed", str(d / "GRCh38.repeats.bed.gz")], check=True)
    tr = Tracks(tmp_path)
    assert tr.has_repeats
    labels = tr.repeat_labels([(1, "chr1", 1000, 1250), (2, "chr1", 5000, 5100), (3, "chr1", 9000, 9100), (4, "chr1", 1, 200_000)])
    assert labels[1][0] == "Alu" and labels[1][1] == 1.0
    assert labels[2][0] == "Tandem repeat"
    assert 3 not in labels and 4 not in labels  # nothing there; and 200 kb is too big for a repeat label


# ---------------------------------------------------------------- trio structure
def test_violation_rules():
    v = index_mod._violation
    assert v(1, 0, 0) and v(2, 0, 1) and v(2, 1, 0) and v(0, 2, 2)
    assert not v(1, 1, 0) and not v(2, 1, 1) and not v(0, 1, 1) and not v(0, 0, 0)


def test_inheritance_classes():
    inh = index_mod.inheritance
    assert inh(0, 0, 0) is None and inh(-1, 1, 1) is None and inh(1, None, 0) is None
    assert inh(1, 0, 0) == "de_novo" and inh(1, 1, 0) == "inherited" and inh(2, 1, 1) == "inherited"
    assert inh(1, -1, 0) == "unknown" and inh(1, -1, 1) == "inherited"


def trio_calls(n_per_pattern=60, patterns=None):
    patterns = patterns or [(1, 1, 0), (1, 0, 1), (0, 1, 0), (0, 0, 1)]  # kid, mom, dad copies; the last two are not transmitted
    gt = {0: "0/0", 1: "0/1", 2: "1/1"}
    rows, pos = [], 1000
    for k, pattern in enumerate(patterns):
        for _ in range(n_per_pattern):
            pos += 1000
            rows.append(rec("chr9", pos, f"f{pos}", "A" * 100, "A", "PASS", f"SVTYPE=DEL;END={pos + 99};SVLEN=-99", [f"{gt[c]}:40:20,10" for c in pattern]))
    return rows


def calls_of(tmp_path, rows, header=HEADER):
    p = tmp_path / "t.vcf"
    write_vcf(p, rows, header=header)
    stats = vcf.ReadStats()
    return list(vcf.iter_calls(str(p), stats)), stats.sample_names


def test_the_child_of_a_trio_is_inferred_from_mendelian_consistency(tmp_path):
    calls, names = calls_of(tmp_path, trio_calls())
    roles = index_mod.infer_trio(calls, names)
    assert (roles["proband"], roles["parents"], roles["source"]) == (0, [1, 2], "inferred")
    assert "KID is the child" in roles["basis"] and "0.0%" in roles["basis"]
    # the same genotypes with the sample order changed: the child is found wherever it is
    shuffled = HEADER.replace("KID\tMOM\tDAD", "DAD\tKID\tMOM")
    rows = [(*r[:9], r[11], r[9], r[10]) for r in trio_calls()]
    calls, names = calls_of(tmp_path, rows, header=shuffled)
    roles = index_mod.infer_trio(calls, names)
    assert names[roles["proband"]] == "KID" and sorted(names[i] for i in roles["parents"]) == ["DAD", "MOM"]


def test_no_child_is_chosen_when_the_evidence_does_not_separate_them(tmp_path):
    only_consistent = trio_calls(n_per_pattern=120, patterns=[(1, 1, 0), (1, 0, 1)])  # 240 sites; every sample is equally consistent as a child
    calls, names = calls_of(tmp_path, only_consistent)
    roles = index_mod.infer_trio(calls, names)
    assert roles["proband"] is None and "no sample is clearly the child" in roles["basis"]
    calls, names = calls_of(tmp_path, trio_calls(n_per_pattern=10))  # 40 sites: too few
    assert "too few" in index_mod.infer_trio(calls, names)["basis"]
    two = HEADER.replace("\tDAD", "")
    calls, names = calls_of(tmp_path, [(*r[:11],) for r in trio_calls()], header=two)
    assert index_mod.infer_trio(calls, names)["proband"] is None


def test_declared_roles_win_over_inference(tmp_path):
    calls, names = calls_of(tmp_path, trio_calls())
    roles = index_mod.resolve_roles(calls, names, {"proband": "MOM", "parents": ["KID", "DAD"]})
    assert (roles["proband"], roles["parents"], roles["source"]) == (1, [0, 2], "sample.json") and "stated in sample.json" in roles["basis"]
    assert index_mod.resolve_roles(calls, names, {})["source"] == "inferred"


# ---------------------------------------------------------------- end to end
def hand_records():
    """Ten hand-made calls on a trio; the gene table and catalog above give them known annotations."""
    return [
        rec("chr1", 10_000, "a", "T" + "C" * 399, "T", "PASS", "SVTYPE=DEL;END=10400;SVLEN=-400", ["0/1:60:20,15", "0/0:40:30,0", "0/0:40:28,0"]),
        rec("chr1", 15_000, "b", "T" + "C" * 299, "T", "PASS", "SVTYPE=DEL;END=15300;SVLEN=-300", ["0/1:50:25,12", "0/1:40:20,10", "0/0:40:30,0"]),
        rec("chr1", 99_000, "c", "T" + "C" * 2499, "T", "PASS", "SVTYPE=DEL;END=101500;SVLEN=-2500", ["1/1:99:0,30", "0/1:40:20,10", "0/1:40:20,10"]),
        rec("chr1", 20_200, "d", "T", "T" + "A" * 149, "PASS", "SVTYPE=INS;END=20200;SVLEN=150", ["0/1:45:20,10", "0/1:40:20,10", "0/0:40:30,0"]),
        rec("chr1", 62_000, "e", "T", "T" + "A" * 799, "PASS", "SVTYPE=INS;END=62000;SVLEN=800", ["0/1:20:10,9", "0/0:10:3,0", "0/0:10:4,0"]),
        rec("chr2", 1_000_000, "f", "C", "<DUP>", "PASS", "SVTYPE=DUP;END=1200000;SVLEN=200000", ["0/1:70:20,8", "0/1:40:20,10", "0/0:40:30,0"]),
        rec("chr1", 10_100, "g", "A", "]chr5:500000]A", "PASS", "SVTYPE=BND;MATEID=g2", ["0/1:40:20,10", "0/0:40:30,0", "0/0:40:30,0"]),
        rec("chr3", 100, "h", "A", "<INV>", "ConflictingBreakpointGT", "SVTYPE=INV;END=5000;SVLEN=4900", ["1/1:30:0,20", "0/1:30:10,10", "0/1:30:10,10"]),
        rec("chr1", 40_000, "i", "T" + "C" * 99, "T", "PASS", "SVTYPE=DEL;END=40100;SVLEN=-100", ["./1:30:20,10", "0/0:40:30,0", "0/0:40:30,0"]),
        rec("chr1", 41_000, "j", "T" + "C" * 99, "T", "PASS", "SVTYPE=DEL;END=41100;SVLEN=-100", ["0/0:40:30,0", "0/1:40:20,10", "0/0:40:30,0"]),
    ]


@pytest.fixture
def env(tmp_path, monkeypatch, bundle):
    root = tmp_path / "samples"
    (root / "FAM").mkdir(parents=True)
    write_vcf(root / "FAM" / "joint.vcf", hand_records() + trio_calls(n_per_pattern=100))
    (root / "FAM" / "sample.json").write_text(json.dumps({"label": "Synthetic family", "caller": "testcaller"}))
    write_genes(tmp_path / "genes.tsv.gz")
    write_catalog(tmp_path / "cat.tsv.gz")
    monkeypatch.setattr(samples, "RESULTS_ROOT", root)
    monkeypatch.setattr(samples, "CACHE_DIR", tmp_path / "cache")
    monkeypatch.setattr(samples, "GENES_PATH", str(tmp_path / "genes.tsv.gz"))
    monkeypatch.setattr(samples, "CATALOGS_PATH", str(tmp_path / "cat.tsv.gz"))
    monkeypatch.setattr(samples, "BUNDLE_PATH", str(bundle))
    for reg in (samples._resources, samples._overview_cache, samples._threads, samples._progress, samples._errors):
        reg.clear()
    return tmp_path


def build(key="FAM"):
    c = samples.resolve(key)
    assert samples.ensure_index(c)["state"] == "building"
    samples._threads[c.key].join(timeout=60)
    assert samples.status(c)["state"] == "ready", samples.status(c)
    return c


def run(c, **kw):
    f = samples.Filters(**{k: v for k, v in kw.items() if k in samples.Filters._fields})
    return samples.query_svs(c, f, **{k: v for k, v in kw.items() if k not in samples.Filters._fields})


def by_vid(result):
    return {r["vid"]: r for r in result["rows"]}


def test_index_reads_the_family_and_assigns_roles(env):
    c = build()
    qc = samples.overview(c)["qc"]
    assert qc["invalid_records"] == 0 and qc["records"] == 410 and qc["source"] == "testcaller 1.0" and qc["samples"] == ["KID", "MOM", "DAD"]
    assert qc["roles"]["proband"] == "KID" and qc["roles"]["parents"] == ["MOM", "DAD"] and qc["roles"]["source"] == "inferred"
    assert qc["by_type_all"]["BND"] == 1 and qc["by_type_pass"]["INV"] == 0 and qc["pass_non_bnd"] == 408
    assert {f["filter"]: f["count"] for f in qc["filters"]}["ConflictingBreakpointGT"] == 1
    kid = next(s for s in qc["per_sample"] if s["sample"] == "KID")
    assert kid["role"] == "proband" and kid["carrier"] + kid["ref"] + kid["nocall"] == 408
    assert qc["de_novo_candidates"] == 3  # a, e and the half-called i; every filler call is inherited
    assert qc["size_bins"][0]["bin"] == "50–99" and sum(sum(v for k, v in b.items() if k != "bin") for b in qc["size_bins"]) == 408


def test_default_view_is_the_proband_s_passing_non_breakend_calls(env):
    c = build()
    r = run(c, page_size=300)
    assert r["total_all"] == 410 and r["total_matching"] == 207 and r["trio"]
    ids = {x["vid"] for x in r["rows"]} & set("abcdefghij")
    assert ids == {"a", "b", "c", "d", "e", "f", "i"}  # g is a breakend, h fails its FILTER, j is not the proband's
    assert run(c, include_bnd=True)["total_matching"] == 208 and run(c, pass_only=False)["total_matching"] == 208
    assert run(c, carried=False)["total_matching"] == 408
    f = run(c, explain=True)["funnel"]
    assert [x["count"] for x in f] == [410, 409, 408, 207]


def test_gene_impact_and_gene_filters(env):
    c = build()
    rows = by_vid(run(c, page_size=300))
    assert (rows["a"]["impact"], rows["a"]["genes"][0]["exons"], rows["a"]["genes"][0]["transcript"]) == ("coding_exon", [2], "NM_000001.1")
    assert (rows["b"]["impact"], rows["c"]["impact"], rows["d"]["impact"], rows["e"]["impact"], rows["f"]["impact"]) == ("intronic", "whole_gene", "coding_exon", "intronic", "none")
    assert rows["c"]["genes"][0]["symbol"] == "G3" and rows["f"]["genes"] == [] and rows["i"]["impact"] == "none"
    assert run(c, impact="gene")["total_matching"] == 5 and run(c, impact="exon")["total_matching"] == 3
    assert run(c, impact="coding")["total_matching"] == 3 and run(c, impact="whole_gene")["total_matching"] == 1
    assert run(c, genes=("g1",))["total_matching"] == 3 and run(c, genes=("G1", "G2"))["total_matching"] == 4  # a, b, d | + e


def test_inheritance_and_parental_depth(env):
    c = build()
    assert {x["vid"] for x in run(c, inheritance="de_novo")["rows"]} == {"a", "e", "i"}
    assert {x["vid"] for x in run(c, inheritance="de_novo", min_parent_dp=10)["rows"]} == {"a", "i"}  # e: parents read only 3-4 times
    assert run(c, inheritance="inherited")["total_matching"] == 204 and run(c, inheritance="unknown")["total_matching"] == 0
    rows = by_vid(run(c, page_size=300))
    assert rows["a"]["inheritance"] == "de_novo" and rows["a"]["parent_min_dp"] == 28 and rows["e"]["parent_min_dp"] == 3
    fam = {m["sample"]: m for m in rows["a"]["family"]}
    assert (fam["KID"]["role"], fam["KID"]["alt"], fam["MOM"]["role"], fam["MOM"]["copies"]) == ("proband", 15, "parent", 0)
    assert rows["i"]["focus"]["copies"] == 1 and fam["KID"]["partial"] is False


def test_type_size_zygosity_and_quality_filters(env):
    c = build()
    assert run(c, types=("INS",))["total_matching"] == 2 and run(c, types=("DUP",))["total_matching"] == 1 and run(c, types=("DEL",))["total_matching"] == 204
    assert {x["vid"] for x in run(c, min_size=1000)["rows"]} == {"c", "f"}
    assert run(c, max_size=200)["total_matching"] == 202
    assert [x["vid"] for x in run(c, zygosity="hom")["rows"]] == ["c"] and run(c, zygosity="het")["total_matching"] == 206
    assert {x["vid"] for x in run(c, min_gq=50)["rows"]} == {"a", "b", "c", "f"}  # filler calls are GQ 40
    assert {x["vid"] for x in run(c, min_alt=12)["rows"]} == {"a", "b", "c"}
    assert run(c, chrom="chr2")["total_matching"] == 1 and run(c, precise_only=True)["total_matching"] == 207


def test_population_and_track_filters_and_flags(env):
    c = build()
    rows = by_vid(run(c, page_size=300))
    assert rows["b"]["pop_af_max"] == 0.31 and set(rows["b"]["population"]) == {"gnomad_sv"} and rows["b"]["segdup_frac"] == 1.0
    assert rows["c"]["population"]["hprc"]["af"] == 0.04 and rows["d"]["population"]["decode"]["svtype"] == "INS" and rows["d"]["pop_af_max"] is None
    assert rows["a"]["population"] == {} and rows["a"]["pop_n"] == 0
    assert run(c, max_pop_af=0.1)["total_matching"] == 206  # b (0.31) is common
    assert run(c, novel=True)["total_matching"] == 204  # b, c and d each match a catalog
    assert run(c, max_segdup=0.5)["total_matching"] == 206
    odd = by_vid(run(c, pass_only=False, carried=False, page_size=300))
    assert odd["h"]["odd"] == ["gap"] and run(c, pass_only=False, carried=False, no_odd=True)["total_matching"] == 408
    assert rows["a"]["repeat"] is None  # the repeat file is absent here


def test_empty_results_explain_which_filter_removed_the_calls(env):
    c = build()
    r = run(c, inheritance="de_novo", genes=("G3",))
    assert r["total_matching"] == 0
    labels = [(f["label"], f["count"]) for f in r["funnel"]]
    assert labels[-2:] == [("de novo candidate", 3), ("in G3", 0)]


def test_sorting_and_paging(env):
    c = build()
    pos = [(x["chrom"], x["pos"]) for x in run(c, page_size=300)["rows"]]
    assert pos == sorted(pos)
    assert run(c, sort="size", order="desc", page_size=3)["rows"][0]["vid"] == "f"
    top = run(c, sort="impact", order="desc", page_size=3)["rows"]
    assert [x["vid"] for x in top][:1] == ["c"] and [x["impact"] for x in top][0] == "whole_gene"
    assert run(c, page_size=50, page=2)["rows"][0]["id"] != run(c, page_size=50, page=1)["rows"][0]["id"]
    assert len(run(c, page_size=1000)["rows"]) == samples.MAX_PAGE


def test_overview_funnel_and_resources(env):
    c = build()
    o = samples.overview(c)
    assert [(f["label"], f["count"]) for f in o["funnel"]] == [
        ("All calls in the VCF", 410), ("Pass the caller's FILTER", 409), ("Not a breakend", 408), ("Carried by the proband", 207),
        ("Rare: population AF ≤ 1% or not found", 205), ("Hits a gene exon or a whole gene", 2),
        ("Outside segmental duplications and repeats", 2), ("De novo candidates", 1),  # b (0.31) and c (0.04) are common
    ]
    assert {g["gene"] for g in o["top_genes"]} == {"G1"}  # c, the whole-gene deletion of G3, matches HPRC at 4%
    inst = {r["key"]: r["installed"] for r in o["resources"]}
    assert inst["genes"] and inst["gnomad_sv"] and inst["hprc"] and inst["tracks"] and not inst["clingen"] and not inst["insilico"]
    assert "not evaluated" in next(r["detail"] for r in o["resources"] if r["key"] == "clingen")
    assert "repeats not evaluated" in next(r["detail"] for r in o["resources"] if r["key"] == "tracks")


def test_genes_hit_by_rare_calls_lists_protein_coding_genes_only(env):
    c = samples.resolve("FAM")
    write_vcf(Path(c.vcf), [
        rec("chr1", 10_000, "coding", "T" + "C" * 399, "T", "PASS", "SVTYPE=DEL;END=10400;SVLEN=-400", ["0/1:60:20,15", "0/0:40:30,0", "0/0:40:28,0"]),
        rec("chr1", 100_100, "lnc", "T" + "C" * 99, "T", "PASS", "SVTYPE=DEL;END=100200;SVLEN=-100", ["0/1:60:20,15", "0/0:40:30,0", "0/0:40:28,0"]),
    ] + trio_calls(n_per_pattern=100))
    c = build()
    rows = by_vid(run(c, impact="exon", page_size=10))
    assert rows["lnc"]["genes"][0]["biotype"] == "lncRNA" and rows["lnc"]["impact"] == "noncoding_exon"  # an exon hit, but not a protein-coding gene
    assert [g["gene"] for g in samples.overview(c)["top_genes"]] == ["G1"]


def test_a_resource_that_is_missing_is_reported_and_its_filters_are_left_out(env, monkeypatch):
    monkeypatch.setattr(samples, "GENES_PATH", str(env / "missing_genes.tsv.gz"))
    monkeypatch.setattr(samples, "CATALOGS_PATH", str(env / "missing_catalogs.tsv.gz"))
    c = build()
    r = run(c, impact="exon", max_pop_af=0.01, novel=True)
    assert set(r["ignored_filters"]) == {"impact", "max_pop_af", "novel"} and not r["genes_installed"] and not r["population_installed"]
    assert r["total_matching"] == 207  # nothing was silently filtered out
    assert all(x["genes"] == [] and x["population"] == {} for x in r["rows"])
    assert [f["label"] for f in samples.overview(c)["funnel"]][-2:] == ["Carried by the proband", "De novo candidates"]


def test_without_parents_there_is_no_inheritance(env):
    (env / "samples" / "FAM" / "sample.json").write_text(json.dumps({"proband": "KID"}))  # named, but with no parents
    c = build()
    r = run(c, inheritance="de_novo")
    assert not r["trio"] and "inheritance" in r["ignored_filters"] and r["total_matching"] == 207
    assert samples.overview(c)["qc"]["de_novo_candidates"] is None
    assert all(x["inheritance"] is None for x in run(c, page_size=50)["rows"])


def test_region_view_returns_mane_exons_and_the_probands_calls(env):
    c = build()
    r = samples.region(c, "chr1", 1, 110_000)
    assert [g["symbol"] for g in r["genes"]] == ["G1", "G2", "G3"] and r["genes"][0]["exons"][1] == (10000, 10500)
    assert {s["id"] for s in r["svs"]} >= {1, 2, 3} and all(s["copies"] > 0 for s in r["svs"])
    assert samples.region(c, "chr9", 1, 10_000)["genes"] == []


def test_gene_info(env):
    c = build()
    g = samples.gene_info("G1")
    assert (g["transcript"], g["model_source"], g["pli"], g["inheritance"], g["hpo_n"], g["diseases"], g["n_exons"]) == ("NM_000001.1", "MANE Select", 0.99, ["AD"], 3, ["Disease one"], 3)
    assert samples.gene_info("NOPE") is None and c


def test_index_rebuilds_when_the_vcf_a_resource_or_the_declared_roles_change(env):
    c = build()
    assert samples.status(c)["state"] == "ready"
    write_vcf(Path(c.vcf), hand_records() + trio_calls(n_per_pattern=100)[:50])
    assert samples.status(c)["state"] == "stale"
    c = build()
    assert index_mod.read_meta(samples.db_path(c))["n_svs"] == "60"
    write_catalog(env / "cat.tsv.gz")  # rewritten: the stamp changes
    assert samples.status(c)["state"] == "stale"
    build()
    (env / "samples" / "FAM" / "sample.json").write_text(json.dumps({"proband": "MOM", "parents": ["KID", "DAD"]}))
    assert samples.status(samples.resolve("FAM"))["state"] == "stale"
    assert samples.overview(build())["qc"]["roles"]["proband"] == "MOM"


def test_a_truncated_vcf_is_reported_clearly_not_partly_read(env):
    c = samples.resolve("FAM")
    full = Path(c.vcf)
    gz = full.with_suffix(".vcf.gz")
    with gzip.open(gz, "wb") as f:
        f.write(full.read_bytes())
    full.unlink()
    gz.write_bytes(gz.read_bytes()[: gz.stat().st_size // 2])
    samples._resources.clear()
    c = samples.resolve("FAM")
    samples.ensure_index(c)
    samples._threads[c.key].join(timeout=30)
    st = samples.status(c)
    assert st["state"] == "error" and "truncated or corrupt" in st["error"] and "gzip -t" in st["error"]


# ---------------------------------------------------------------- API
def test_api_flow(env):
    with TestClient(app) as client:
        assert client.get("/api/structvar/samples").json() == {"samples": [{"key": "FAM", "label": "Synthetic family"}], "default": "FAM"}
        assert client.get("/api/structvar/svs", params={"sample": "nope"}).status_code == 404
        assert client.get("/api/structvar/status").json()["state"] in ("building", "ready")
        samples._threads["FAM"].join(timeout=60)
        assert client.get("/api/structvar/status").json() == {"state": "ready"}
        meta = client.get("/api/structvar/meta").json()
        assert meta["label"] == "Synthetic family" and meta["n_svs"] == 410 and meta["roles"]["proband"] == "KID"
        r = client.get("/api/structvar/svs", params={"inheritance": "de_novo", "impact": "exon", "types": "DEL,INS"}).json()
        assert [x["vid"] for x in r["rows"]] == ["a"]
        one = client.get(f"/api/structvar/svs/{r['rows'][0]['id']}").json()
        assert one["genes"][0]["symbol"] == "G1" and len(one["family"]) == 3
        assert client.get("/api/structvar/svs/999999").status_code == 404
        for bad in ({"zygosity": "x"}, {"inheritance": "x"}, {"impact": "x"}):
            assert client.get("/api/structvar/svs", params=bad).status_code == 422
        assert client.get("/api/structvar/region", params={"chrom": "chr1", "start": 1, "end": 110000}).json()["n_genes"] == 3
        assert client.get("/api/structvar/region", params={"chrom": "chr1", "start": 1, "end": 30_000_000}).status_code == 422
        assert client.get("/api/structvar/genes/G2").json()["transcript"] == "NM_000002.1"
        assert client.get("/api/structvar/genes/NOPE").status_code == 404
        assert client.get("/api/structvar/overview").json()["funnel"][0]["count"] == 410


def test_api_answers_503_until_the_index_is_built(env, monkeypatch):
    release = threading.Event()
    real = index_mod.build_index

    def slow(*a, **k):
        release.wait(timeout=30)
        return real(*a, **k)

    monkeypatch.setattr(index_mod, "build_index", slow)
    with TestClient(app) as client:
        assert client.get("/api/structvar/status").json()["state"] == "building"
        resp = client.get("/api/structvar/svs")
        assert resp.status_code == 503 and resp.json()["detail"]["state"] == "building"
        release.set()
        samples._threads["FAM"].join(timeout=60)
        assert client.get("/api/structvar/svs").status_code == 200
