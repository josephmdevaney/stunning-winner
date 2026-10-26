#!/usr/bin/env python3
"""Build the Pangenome Repeat Catalog and its HG002 genotypes.

The catalog is the set of TRExplorer v2.0 variation-cluster loci that the
pangenome SV paper's Supplementary Table 3 (Lin et al., medRxiv
10.64898/2026.08.21.26361050; tables 3a, 3b, 3d, 3e, 3f) found in TRExplorer.
Each found paper locus maps to one variation-cluster locus, because that is the
catalog TRGT genotypes against:

  1. the variation-cluster locus that contains the locus the TRExplorer portal matched
     to the paper locus, by ID (the usual case; a cluster's ID is the comma-joined IDs
     of the TRs merged into it);
  2. otherwise the variation-cluster locus that overlaps the paper locus's
     coordinates by the most bases (the portal's locus was merged into a larger
     cluster, or the paper row had no portal match);
  3. otherwise the paper locus is left out and reported -- a locus with no
     variation-cluster definition can't be genotyped, and no definition is
     invented for it.

Paper rows with no coordinates at all (gene name only) are left out and
reported too.

Input is the "found" rows as CSV with these columns (the "In Repeat Catalog" and
"Not in Repeat Catalog" tabs of the comparison sheet, concatenated):
    Source table, Paper locus_id, Chrom, Start (0-based), End,
    TRExplorer portal match: locus ID

The full variation-cluster catalog BED and the full HG002 genome-wide TRGT VCF
(bgzip + .tbi) are not kept in the repo; see project/uploads/README.md for where
they come from. Requires `bcftools` and `tabix`.

Usage:
    python3 build_pangenome_catalog.py \
        --found-loci-csv found_loci.csv \
        --catalog-bed /path/to/full/repeat_catalog.hg38.bed \
        --vcf /path/to/full/HG002.GRCh38.trgt.genomewide.vcf.gz \
        --out-catalog-bed project/uploads/pangenome_repeat_catalog.hg38.bed \
        --out-vcf backend/fixtures/hg002_real_genomewide/HG002.GRCh38.trgt.pangenome.vcf.gz
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

PORTAL_ID = "TRExplorer portal match: locus ID"


def chrom_key(chrom: str) -> tuple[int, str]:
    name = chrom.removeprefix("chr")
    if name.isdigit():
        return (int(name), "")
    return ({"X": 23, "Y": 24, "M": 25}.get(name, 99), name)


def locus_id(info: str) -> str:
    for field in info.split(";"):
        if field.startswith("ID="):
            return field[3:]
    return ""


def read_found_loci(path: str) -> tuple[list[dict], list[dict]]:
    """(rows with coordinates, rows without -- gene-name-only paper entries)."""
    with_coords, without = [], []
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            (with_coords if row["Start (0-based)"] and row["End"] else without).append(row)
    return with_coords, without


def resolve(found: list[dict], catalog_bed: str) -> tuple[dict[str, str], list[dict]]:
    """Single pass over the full catalog. Returns ({locus ID: BED line}, paper rows left out).

    A variation-cluster locus's ID is the comma-joined IDs of the TRs merged into it, so a
    portal ID is matched against each of those constituents, and the cluster is the locus."""
    wanted_ids = {r[PORTAL_ID] for r in found if r[PORTAL_ID]}
    by_id: dict[str, str] = {}
    matched_portal_ids: set[str] = set()

    with open(catalog_bed, encoding="utf-8") as f:
        for line in f:
            lid = locus_id(line.rstrip("\n").split("\t")[3])
            hits = wanted_ids.intersection(lid.split(","))
            if hits:
                matched_portal_ids |= hits
                by_id[lid] = line if line.endswith("\n") else line + "\n"

    # Paper rows whose portal ID isn't in any variation-cluster locus are resolved by
    # overlap, which needs every catalog line on their chromosomes.
    unmatched_by_chrom: dict[str, list[dict]] = defaultdict(list)
    for row in found:
        if row[PORTAL_ID] not in matched_portal_ids:
            unmatched_by_chrom[row["Chrom"]].append(
                {"row": row, "start": int(row["Start (0-based)"]), "end": int(row["End"]), "best": (0, None)}
            )

    if unmatched_by_chrom:
        with open(catalog_bed, encoding="utf-8") as f:
            for line in f:
                chrom, start, stop, _info = line.rstrip("\n").split("\t")
                pending = unmatched_by_chrom.get(chrom)
                if not pending:
                    continue
                start_i, stop_i = int(start), int(stop)
                for p in pending:
                    overlap = min(p["end"], stop_i) - max(p["start"], start_i)
                    if overlap > p["best"][0]:
                        p["best"] = (overlap, line if line.endswith("\n") else line + "\n")

    left_out = []
    for pending in unmatched_by_chrom.values():
        for p in pending:
            _overlap, line = p["best"]
            if line is None:
                left_out.append(p["row"])
            else:
                by_id[locus_id(line.rstrip("\n").split("\t")[3])] = line
    return by_id, left_out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--found-loci-csv", required=True)
    ap.add_argument("--catalog-bed", required=True, help="full TRExplorer v2.0 variation-cluster catalog BED")
    ap.add_argument("--vcf", required=True, help="full genome-wide TRGT VCF against that catalog (bgzip + .tbi)")
    ap.add_argument("--out-catalog-bed", required=True)
    ap.add_argument("--out-vcf", required=True)
    args = ap.parse_args()

    found, no_coords = read_found_loci(args.found_loci_csv)
    print(f"found paper rows with coordinates: {len(found)}; gene-name-only rows left out: {len(no_coords)}", file=sys.stderr)

    by_id, left_out = resolve(found, args.catalog_bed)
    for row in left_out:
        print(
            f"left out (no variation-cluster locus overlaps it): table {row['Source table']} "
            f"{row['Paper locus_id']} {row['Chrom']}:{row['Start (0-based)']}-{row['End']}",
            file=sys.stderr,
        )

    def sort_key(line: str):
        chrom, start, stop, _ = line.split("\t", 3)
        return (chrom_key(chrom), int(start), int(stop))

    lines = sorted(by_id.values(), key=sort_key)
    with open(args.out_catalog_bed, "w", encoding="utf-8") as out:
        out.writelines(lines)
    print(f"catalog loci written: {len(lines)}", file=sys.stderr)

    # Extract the HG002 genotypes for exactly these loci: -R selects by position,
    # the TRID filter drops any neighbouring locus the position window also catches.
    with tempfile.TemporaryDirectory() as tmp:
        regions = Path(tmp) / "regions.bed"
        ids = Path(tmp) / "ids.txt"
        regions.write_text("".join(line.split("\t", 3)[0] + "\t" + line.split("\t")[1] + "\t" + line.split("\t")[2] + "\n" for line in lines))
        ids.write_text("".join(locus_id(line.rstrip("\n").split("\t")[3]) + "\n" for line in lines))
        subprocess.run(
            f"bcftools view -R {regions} {args.vcf} | bcftools view -i 'INFO/TRID=@{ids}' -Oz -o {args.out_vcf}",
            shell=True, check=True,
        )
    subprocess.run(["tabix", "-f", "-p", "vcf", args.out_vcf], check=True)
    n_records = int(
        subprocess.run(
            f"bcftools view -H {args.out_vcf} | wc -l", shell=True, check=True, capture_output=True, text=True
        ).stdout.strip()
    )
    print(f"subset VCF records: {n_records}", file=sys.stderr)


if __name__ == "__main__":
    main()
