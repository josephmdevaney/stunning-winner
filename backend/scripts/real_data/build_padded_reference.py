#!/usr/bin/env python3
"""Build an N-padded per-chromosome reference FASTA, with real GRCh38 sequence
patched in only at the STRchive loci windows (+ padding). Chromosome lengths match
the real GCA_000001405.15_GRCh38_no_alt_analysis_set reference exactly, so trgt's
consistency checks against the BAM header pass, while we only need to fetch a tiny
fraction of the genome's real sequence over the network.
"""
import subprocess
import sys

REF_URL = "https://giab.s3.amazonaws.com/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz"
FAI = "ref.fasta.gz.fai"
BED = "strchive_trgt.bed"
PAD = 2000
OUT = "padded_reference.fasta"
LINE_WIDTH = 70


def load_lengths():
    lengths = {}
    with open(FAI) as f:
        for line in f:
            parts = line.split("\t")
            lengths[parts[0]] = int(parts[1])
    return lengths


def load_loci_by_chrom():
    by_chrom = {}
    with open(BED) as f:
        for line in f:
            chrom, start, stop = line.split("\t")[:3]
            by_chrom.setdefault(chrom, []).append((int(start), int(stop)))
    return by_chrom


def fetch_seq(chrom, start, end, retries=4):
    """1-based inclusive region fetch via remote-indexed faidx; returns uppercase seq, no newlines."""
    region = f"{chrom}:{start}-{end}"
    last_err = None
    for attempt in range(retries):
        proc = subprocess.run(
            ["samtools", "faidx", REF_URL, region], capture_output=True, text=True
        )
        if proc.returncode == 0:
            lines = proc.stdout.splitlines()[1:]  # drop ">header" line
            return "".join(lines).upper()
        last_err = proc.stderr
        print(f"  retry {attempt+1}/{retries} for {region}: {last_err.strip()}", file=sys.stderr)
    raise RuntimeError(f"failed to fetch {region} after {retries} attempts: {last_err}")


def main():
    lengths = load_lengths()
    by_chrom = load_loci_by_chrom()

    with open(OUT, "w") as out:
        for chrom in sorted(by_chrom, key=lambda c: list(lengths.keys()).index(c)):
            length = lengths[chrom]
            seq = bytearray(b"N" * length)

            for start, stop in by_chrom[chrom]:
                win_start = max(1, start - PAD)  # 1-based inclusive for faidx
                win_end = min(length, stop + PAD)
                real = fetch_seq(chrom, win_start, win_end)
                # seq is 0-based; win_start is 1-based -> offset win_start-1
                seq[win_start - 1 : win_start - 1 + len(real)] = real.encode()
                print(f"  patched {chrom}:{win_start}-{win_end} ({len(real)}bp)", file=sys.stderr)

            print(f"writing {chrom} ({length:,}bp, {len(by_chrom[chrom])} loci patched)", file=sys.stderr)
            out.write(f">{chrom}\n")
            seq_str = seq.decode()
            for i in range(0, len(seq_str), LINE_WIDTH):
                out.write(seq_str[i : i + LINE_WIDTH])
                out.write("\n")

    print("done", file=sys.stderr)


if __name__ == "__main__":
    main()
