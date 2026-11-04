#!/usr/bin/env python3
"""Build the gene-model table the structural-variant suite annotates against.

Inputs are the Ensembl 101 gene models and the gene-level files from PacBio's HiFi WDL resource
bundle (v1.0.2, ``static_resources/``), the same files its svpack annotation uses:

  GRCh38/ensembl.GRCh38.101.reformatted.gff3.gz   gene models. Despite the name this file is plain
                                                  text in the bundle; compression is detected from
                                                  the first bytes, not the extension.
  hpo/ensembl.inheritance.tsv                     Ensembl gene -> mode of inheritance
  hpo/ensembl.hpoPhenotype.tsv, hpo/hpoTerms.txt  Ensembl gene -> HPO phenotype terms
  slivar/lof_lookup.v2.1.1.txt                    gnomAD v2.1.1 pLI and observed/expected LoF
  slivar/clinvar_gene_desc.*.txt                  gene -> ClinVar condition names (as of 2022-12)

  ~/sv_data/ucsc/mane.bb                          UCSC's MANE track (https://hgdownload.soe.ucsc.edu/gbdb/hg38/mane/
                                                  mane.bb): the MANE Select transcript of each gene, from the NCBI /
                                                  EMBL-EBI MANE project, with its RefSeq NM_ accession

One row per gene, with one transcript whose exons are what an SV is compared against. That is the MANE
Select transcript wherever the gene has one (about 19,400 mostly protein-coding genes), so exon numbers match
the numbering clinical reports use. A gene with no MANE transcript keeps the Ensembl 101 transcript with the
longest coding sequence (non-coding genes: the longest transcript), and its ``model_source`` says so. The
gene's span is the Ensembl 101 gene locus (widened to cover the MANE transcript if needed). The gene-level
columns are context shown next to a hit.

Usage:
    python3 build_sv_genes.py --bundle ~/sv_data/pacbio_bundle/static_resources \
        --mane ~/sv_data/ucsc/mane.bb --out ~/sv_data/sv_genes.grch38.tsv.gz
"""

from __future__ import annotations

import argparse
import gzip
import sys
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from structvar.bigbed import BigBed  # noqa: E402
from structvar.normalize import canon_chrom  # noqa: E402

COLUMNS = (
    "gene_id", "symbol", "biotype", "chrom", "start", "end", "strand", "transcript", "ensembl_transcript", "model_source",
    "exons", "cds_start", "cds_end", "pli", "oe_lof", "inheritance", "hpo_n", "hpo_terms", "diseases",
)
GENE_TYPES = {"gene", "ncRNA_gene", "pseudogene"}
TRANSCRIPT_TYPES = {"mRNA", "lnc_RNA", "ncRNA", "pseudogenic_transcript", "transcript", "miRNA", "snRNA", "snoRNA", "rRNA", "tRNA"}
MAX_HPO = 25
MAX_DISEASES = 10


def open_text(path: Path):
    with open(path, "rb") as f:
        magic = f.read(2)
    return gzip.open(path, "rt", encoding="utf-8") if magic == b"\x1f\x8b" else open(path, "rt", encoding="utf-8")


def attrs(field: str) -> dict:
    out = {}
    for kv in field.rstrip("\n").split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = unquote(v)
    return out


def read_gff(path: Path):
    genes, transcripts = {}, {}
    exons = defaultdict(list)
    cds = defaultdict(list)
    with open_text(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t")
            if len(p) < 9:
                continue
            kind = p[2]
            if kind in GENE_TYPES:
                a = attrs(p[8])
                gid = a.get("gene_id") or a.get("ID", "").replace("gene:", "")
                genes[gid] = {"symbol": a.get("Name", gid), "biotype": a.get("biotype", ""), "chrom": canon_chrom(p[0]),
                              "start": int(p[3]), "end": int(p[4]), "strand": p[6]}
            elif kind in TRANSCRIPT_TYPES:
                a = attrs(p[8])
                tid = a.get("transcript_id") or a.get("ID", "").replace("transcript:", "")
                transcripts[tid] = {"gene": a.get("Parent", "").replace("gene:", ""), "biotype": a.get("biotype", ""), "basic": "basic" in a.get("tag", "")}
            elif kind == "exon":
                exons[attrs(p[8]).get("Parent", "").replace("transcript:", "")].append((int(p[3]), int(p[4])))
            elif kind == "CDS":
                cds[attrs(p[8]).get("Parent", "").replace("transcript:", "")].append((int(p[3]), int(p[4])))
    return genes, transcripts, exons, cds


def canonical(gene_id: str, by_gene, exons, cds):
    best, best_key = None, None
    for tid in by_gene.get(gene_id, ()):
        ex = exons.get(tid)
        if not ex:
            continue
        cds_len = sum(e - s + 1 for s, e in cds.get(tid, ()))
        exon_len = sum(e - s + 1 for s, e in ex)
        key = (cds_len, exon_len)
        if best_key is None or key > best_key:
            best, best_key = tid, key
    return best


def read_mane(path: Path) -> dict:
    """MANE Select transcripts from UCSC's bigBed, keyed by unversioned Ensembl gene id."""
    bb = BigBed(path.read_bytes())
    names = bb.field_names()
    out = {}
    for chrom, start0, _end, rest in bb.records():
        d = dict(zip(names[3:], rest))
        if d.get("maneStat") != "MANE Select":
            continue
        sizes = [int(x) for x in d["blockSizes"].split(",") if x]
        rel = [int(x) for x in d["chromStarts"].split(",") if x]
        exons = [(start0 + r + 1, start0 + r + z) for r, z in zip(rel, sizes)]
        thick_s, thick_e = int(d["thickStart"]), int(d["thickEnd"])
        gid = d["geneName"].split(".")[0]
        out[gid] = {
            "symbol": d["geneName2"], "biotype": d["geneType"], "chrom": canon_chrom(chrom), "strand": d["strand"],
            "refseq": d["ncbiId"], "enst": d["name"], "exons": exons,
            "cds": (thick_s + 1, thick_e) if thick_e > thick_s else None,
            "start": exons[0][0], "end": exons[-1][1],
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bundle", required=True, help="the bundle's static_resources directory")
    ap.add_argument("--mane", default="~/sv_data/ucsc/mane.bb", help="UCSC mane.bb (MANE Select transcripts)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    b = Path(args.bundle).expanduser()

    mane_path = Path(args.mane).expanduser()
    if not mane_path.is_file():
        raise SystemExit(f"{mane_path} not found; download https://hgdownload.soe.ucsc.edu/gbdb/hg38/mane/mane.bb")
    mane = read_mane(mane_path)
    genes, transcripts, exons, cds = read_gff(b / "GRCh38" / "ensembl.GRCh38.101.reformatted.gff3.gz")
    # genes in MANE that Ensembl 101 predates are added from MANE alone
    for gid, m in mane.items():
        if gid not in genes:
            genes[gid] = {"symbol": m["symbol"], "biotype": m["biotype"], "chrom": m["chrom"], "start": m["start"], "end": m["end"], "strand": m["strand"]}
    by_gene = defaultdict(list)
    for tid, t in transcripts.items():
        by_gene[t["gene"]].append(tid)

    inheritance = defaultdict(set)
    for line in open(b / "hpo" / "ensembl.inheritance.tsv", encoding="utf-8"):
        g, _, mode = line.rstrip("\n").partition("\t")
        if mode:
            inheritance[g].add(mode)
    terms = {}
    for line in open(b / "hpo" / "hpoTerms.txt", encoding="utf-8"):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            terms[p[0]] = p[1]
    hpo = defaultdict(list)
    for line in open(b / "hpo" / "ensembl.hpoPhenotype.tsv", encoding="utf-8"):
        p = line.rstrip("\n").split("\t")
        name = terms.get(p[1], "")
        if name and "inheritance" not in name.lower():
            hpo[p[0]].append(name)
    constraint = {}
    for line in open(b / "slivar" / "lof_lookup.v2.1.1.txt", encoding="utf-8"):
        sym, _, rest = line.rstrip("\n").partition("\t")
        kv = dict(x.split("=", 1) for x in rest.split(";") if "=" in x)
        constraint[sym] = (kv.get("pLI", ""), kv.get("oe_lof", ""))
    diseases = defaultdict(list)
    desc = next((b / "slivar").glob("clinvar_gene_desc.*.txt"))
    for line in open(desc, encoding="utf-8"):
        p = line.rstrip("\n").split("\t")
        if len(p) == 2 and p[0] != "AssociatedGenes":
            for sym in p[0].split(","):
                diseases[sym.strip()].append(p[1])

    out = Path(args.out).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    n = n_coding = n_mane = 0
    with gzip.open(out, "wt", encoding="utf-8") as f:
        f.write("# MANE Select (UCSC mane.bb) and Ensembl 101 gene models with gene-level context, built by build_sv_genes.py\n")
        f.write("\t".join(COLUMNS) + "\n")
        for gid, g in sorted(genes.items(), key=lambda kv: (kv[1]["chrom"], kv[1]["start"])):
            m = mane.get(gid)
            if m is not None:
                ex, cd = m["exons"], m["cds"]
                display, enst, source = m["refseq"], m["enst"], "MANE Select"
                start_, end_ = min(g["start"], m["start"]), max(g["end"], m["end"])
                n_mane += 1
            else:
                tid = canonical(gid, by_gene, exons, cds)
                ex = sorted(exons.get(tid, ())) if tid else []
                cdr = cds.get(tid, ()) if tid else ()
                cd = (min(s for s, _ in cdr), max(e for _, e in cdr)) if cdr else None
                display = enst = tid or ""
                source = ("Ensembl 101 (longest CDS)" if cdr else "Ensembl 101 (longest transcript)") if tid else ""
                start_, end_ = g["start"], g["end"]
            pli, oe = constraint.get(g["symbol"], ("", ""))
            names = sorted(set(hpo.get(gid, ())))
            f.write("\t".join([
                gid, g["symbol"], g["biotype"], g["chrom"], str(start_), str(end_), g["strand"], display, enst, source,
                ",".join(f"{s_}-{e_}" for s_, e_ in ex), str(cd[0]) if cd else "", str(cd[1]) if cd else "",
                pli, oe, ",".join(sorted(inheritance.get(gid, ()))), str(len(names)), "; ".join(names[:MAX_HPO]),
                "; ".join(list(dict.fromkeys(diseases.get(g["symbol"], ())))[:MAX_DISEASES]),
            ]) + "\n")
            n += 1
            n_coding += g["biotype"] == "protein_coding"
    print(f"{n:,} genes ({n_coding:,} protein-coding), {n_mane:,} with a MANE Select transcript -> {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
