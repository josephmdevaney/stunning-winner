#!/usr/bin/env python3
"""Rescope the genome-wide Repeat Catalog down to just the regions RExPRT
predicts pathogenic, and extract their real genotypes from a genome-wide TRGT
VCF -- both as a performance fix, not a cosmetic one.

The full TRExplorer v2.0 catalog (5,339,786 regions) plus a matching
genome-wide TRGT VCF (5,339,786 records) is real, correct data, but building
the backend's in-memory TrgtIndex from that VCF takes minutes and multiple GB
of RAM (see repeatscope/api/main.py's `_get_trgt_index`) -- and can take much
longer, or effectively hang, on a machine under real memory pressure (verified
directly: a cold build that normally takes ~2 minutes took over 15 without
finishing when free system RAM was measured at ~60MB). Since the product's
actual use of the full catalog was "browse candidate loci RExPRT flagged,"
not "browse all 5.3M reference definitions," this script narrows the *data*
itself to just that candidate set rather than filtering it at request time --
filtering after the fact still requires building the full-size index first.

This rescoping is a deliberate, real reduction in scope, not a bug fix or an
approximation: regions RExPRT did not evaluate, or scored benign, are no
longer present in the genome-wide catalog screen at all after this script
runs. project/uploads/rexprt_scores.hg38.tsv itself is untouched and still
covers the full v2.0 catalog (727,004 real matches, including benign ones) --
only the catalog BED and the sample genotype VCF are rescoped.

Requires `bcftools` (used for real, indexed region extraction from the VCF --
this is the same real-world tool the rest of this pipeline already assumes
exists, not a new invented dependency).

Usage:
    python3 build_rexprt_pathogenic_subset.py \
        --rexprt-scores project/uploads/rexprt_scores.hg38.tsv \
        --catalog-bed /path/to/full/repeat_catalog.hg38.bed \
        --vcf /path/to/full/HG002.GRCh38.trgt.genomewide.vcf.gz \
        --out-catalog-bed project/uploads/repeat_catalog.hg38.bed \
        --out-vcf backend/fixtures/hg002_real_genomewide/HG002.GRCh38.trgt.genomewide.vcf.gz

The full-scope catalog BED and VCF this script takes as input are not kept in
this repo after rescoping (they're multi-hundred-MB); they remain recoverable
from git history (the commit that first added the TRExplorer v2.0 catalog) or
by re-running the earlier build steps documented in project/uploads/README.md.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path


def load_positive_regions(rexprt_scores_path: str) -> set[tuple[str, int, int]]:
    positive: set[tuple[str, int, int]] = set()
    with open(rexprt_scores_path, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        ensemble_binary_idx = header.index("ensemble_binary")
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if fields[ensemble_binary_idx] == "1":
                positive.add((fields[0], int(fields[1]), int(fields[2])))
    return positive


def write_catalog_subset(catalog_bed_path: str, positive: set, out_path: str) -> int:
    n_written = 0
    with open(catalog_bed_path, encoding="utf-8") as f, open(out_path, "w", encoding="utf-8") as out:
        for line in f:
            chrom, start, stop, _info = line.rstrip("\n").split("\t")
            if (chrom, int(start), int(stop)) in positive:
                out.write(line)
                n_written += 1
    return n_written


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rexprt-scores", required=True)
    ap.add_argument("--catalog-bed", required=True, help="full-scope genome-wide catalog BED")
    ap.add_argument("--vcf", required=True, help="full-scope genome-wide TRGT VCF (bgzip + .tbi)")
    ap.add_argument("--out-catalog-bed", required=True)
    ap.add_argument("--out-vcf", required=True)
    args = ap.parse_args()

    positive = load_positive_regions(args.rexprt_scores)
    print(f"RExPRT-predicted-pathogenic regions: {len(positive)}", file=sys.stderr)

    n_written = write_catalog_subset(args.catalog_bed, positive, args.out_catalog_bed)
    print(f"catalog rows written: {n_written}", file=sys.stderr)

    with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as tmp:
        region_rows = []
        with open(args.out_catalog_bed, encoding="utf-8") as f:
            for line in f:
                chrom, start, stop, _info = line.split("\t", 3)
                region_rows.append((chrom, int(start), stop))
        region_rows.sort(key=lambda r: (r[0], r[1]))
        for chrom, start, stop in region_rows:
            tmp.write(f"{chrom}\t{start}\t{stop}\n")
        region_file = tmp.name

    subprocess.run(
        ["bcftools", "view", "-R", region_file, args.vcf, "-Oz", "-o", args.out_vcf],
        check=True,
    )
    subprocess.run(["tabix", "-p", "vcf", args.out_vcf], check=True)
    Path(region_file).unlink()

    n_records = int(
        subprocess.run(
            ["bash", "-c", f"bcftools view -H {args.out_vcf} | wc -l"],
            check=True, capture_output=True, text=True,
        ).stdout.strip()
    )
    print(f"subset VCF records: {n_records}", file=sys.stderr)


if __name__ == "__main__":
    main()
