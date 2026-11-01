#!/usr/bin/env python3
"""Reduce the four population SV VCFs in PacBio's HiFi WDL resource bundle to one table.

From ``static_resources/GRCh38/sv_pop_vcfs/`` (see the README there for each source):

  gnomad_sv   nstd166.GRCh38.variant_call.vcf.gz   gnomAD-SV (Collins 2020) as deposited in dbVar;
                                                   short-read, allele frequency overall and per ancestry
  hprc        hprc.GRCh38.pbsv.vcf.gz              HPRC release 1, pbsv long-read calls, about 100 samples;
                                                   allele frequency (non-PASS calls are skipped)
  decode      ont_sv_high_confidence_SVs.sorted.vcf.gz   DECODE (Beyter 2021), ONT long-read; presence only
  eee         EEE_SV-Pop_1.ALL.sites.20181204.vcf.gz     Audano 2019 long-read discovery set; the number
                                                         of discovery genomes carrying it, no frequency

Output is one row per call: source, chrom, start, end, svtype, svlen, ac, an, af, popmax. Breakends and
the DECODE INSDEL records (both uncommon, and not comparable by size) are left out. A source that has no
frequency leaves af empty; the match is then reported as presence only.

Usage:
    python3 build_sv_catalogs.py --bundle ~/sv_data/pacbio_bundle/static_resources \
        --out ~/sv_data/sv_catalogs.grch38.tsv.gz
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from structvar.normalize import canon_chrom  # noqa: E402

COLUMNS = ("source", "chrom", "start", "end", "svtype", "svlen", "ac", "an", "af", "popmax")
KEEP = {"DEL", "DUP", "INS", "INV", "CNV"}


def info(field: str) -> dict:
    out = {}
    for kv in field.split(";"):
        k, _, v = kv.partition("=")
        out[k] = v
    return out


def first_int(v: str):
    try:
        return int(float(v.split(",")[0]))
    except (ValueError, IndexError):
        return None


def rows(path: Path, source: str):
    with gzip.open(path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t", 8)
            chrom, pos, ref, alt, flt, inf = p[0], int(p[1]), p[3], p[4], p[6], p[7]
            i = info(inf)
            t = i.get("SVTYPE", "")
            if t not in KEEP:
                continue
            if source == "hprc" and flt not in ("PASS", "."):
                continue
            end = first_int(i.get("END", "")) or pos
            svlen = first_int(i.get("SVLEN", ""))
            if svlen is not None:
                svlen = abs(svlen)
            elif t == "INS":
                svlen = max(len(alt) - len(ref), 0) or None
            else:
                svlen = end - pos
            ac = an = af = popmax = ""
            if source == "gnomad_sv":
                ac, an, af = i.get("AC", ""), i.get("AN", ""), i.get("AF", "")
                anc = [float(i[k]) for k in ("AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF") if i.get(k, "") not in ("", ".")]
                popmax = f"{max(anc):.6g}" if anc else ""
            elif source == "hprc":
                ac, an, af = i.get("AC", ""), i.get("AN", ""), i.get("AF", "")
            elif source == "eee":
                ac = i.get("MERGE_AC", "")
            yield (source, canon_chrom(chrom), pos, end, t, svlen if svlen is not None else "", ac, an, af, popmax)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    d = Path(args.bundle).expanduser() / "GRCh38" / "sv_pop_vcfs"
    files = [
        ("gnomad_sv", d / "nstd166.GRCh38.variant_call.vcf.gz"),
        ("hprc", d / "hprc.GRCh38.pbsv.vcf.gz"),
        ("decode", d / "ont_sv_high_confidence_SVs.sorted.vcf.gz"),
        ("eee", d / "EEE_SV-Pop_1.ALL.sites.20181204.vcf.gz"),
    ]
    out = Path(args.out).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    counts = {}
    with gzip.open(out, "wt", encoding="utf-8") as f:
        f.write("# Population SV catalogs reduced by build_sv_catalogs.py from the PacBio HiFi WDL resource bundle\n")
        f.write("\t".join(COLUMNS) + "\n")
        for source, path in files:
            n = 0
            for r in rows(path, source):
                f.write("\t".join(map(str, r)) + "\n")
                n += 1
            counts[source] = n
    print(", ".join(f"{k} {v:,}" for k, v in counts.items()) + f" -> {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
