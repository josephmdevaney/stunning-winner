#!/usr/bin/env python3
"""Build a real trgt-format BED catalog for the 80 STRchive disease loci.

Merges STRchive's coordinates/IDs with PacBio's own official pathogenic-repeat
catalog (from the PacificBiosciences/trgt GitHub repo's repeats/pathogenic_repeats.hg38.bed)
for accurate MOTIFS/STRUC where the two match by gene + coordinates, falling back to a
simple single-motif STRUC built from STRchive's own motif field otherwise.

Requires a local clone of https://github.com/PacificBiosciences/trgt for the official
catalog file:

    git clone --depth 1 https://github.com/PacificBiosciences/trgt.git /tmp/trgt_probe

Usage:
    python3 build_strchive_trgt_catalog.py > strchive_trgt.bed
"""

import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
STRCHIVE_JSON = REPO_ROOT / "project" / "uploads" / "STRchive-loci.json"
PACBIO_OFFICIAL_BED = Path("/tmp/trgt_probe/repeats/pathogenic_repeats.hg38.bed")


def main():
    with open(STRCHIVE_JSON) as f:
        strchive = json.load(f)

    pacbio_official = {}
    with open(PACBIO_OFFICIAL_BED) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            chrom, start, stop, info = parts[0], int(parts[1]), int(parts[2]), parts[3]
            kvs = dict(kv.split("=", 1) for kv in info.split(";"))
            pacbio_official[kvs["ID"]] = (chrom, start, stop, kvs.get("MOTIFS", ""), kvs.get("STRUC", ""))

    rows = []
    matched = 0
    for locus in strchive:
        gene, chrom, start, stop = locus["gene"], locus["chrom"], locus["start_hg38"], locus["stop_hg38"]
        motif = None
        for key in ("pathogenic_motif_reference_orientation", "reference_motif_reference_orientation"):
            vals = locus.get(key) or []
            if vals:
                motif = vals[0]
                break
        motif = motif or "N"

        off = pacbio_official.get(gene)
        if off and off[0] == chrom and abs(off[1] - start) < 200 and abs(off[2] - stop) < 200:
            chrom, start, stop, motifs, struc = off
            matched += 1
        else:
            motifs, struc = motif, f"({motif})n"

        rows.append((chrom, start, stop, f"ID={locus['id']};MOTIFS={motifs};STRUC={struc}"))

    rows.sort(key=lambda r: (r[0], r[1]))
    for chrom, start, stop, info in rows:
        print(f"{chrom}\t{start}\t{stop}\t{info}")

    print(f"wrote {len(rows)} loci, {matched} matched PacBio's official catalog by gene+coords", file=sys.stderr)


if __name__ == "__main__":
    main()
