#!/usr/bin/env python3
"""Generate a small, spec-accurate synthetic TRGT v5.1.0 VCF for development and tests.

Covers a handful of real STRchive loci at benign/intermediate/pathogenic genotypes
(mirroring the original design's narrative overrides) plus a few generic genome-wide
catalog regions, so the annotation pipeline can be developed and tested end-to-end
before a real HG002 TRGT VCF is available. Run with:

    python backend/scripts/generate_synthetic_vcf.py > backend/fixtures/synthetic_trgt.vcf
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
STRCHIVE_JSON = REPO_ROOT / "project" / "uploads" / "STRchive-loci.json"

SAMPLE_NAME = "SYNTH-001"

# (locus_id, a1_units, a2_units, a1_depth, a2_depth, a1_methylation, a2_methylation)
NARRATIVE_CALLS = [
    ("HD_HTT", 18, 44, 16, 15, None, None),
    ("FXS_FMR1", 29, 29, 22, 20, 0.04, 0.05),
    ("FTDALS1_C9orf72", 2, 412, 14, 9, None, 0.71),
    ("SCA1_ATXN1", 29, 37, 18, 17, None, None),
    ("SCA3_ATXN3", 14, 22, 19, 18, None, None),
    ("DM1_DMPK", 5, 5, 24, 23, None, None),
    ("FRDA_FXN", 580, 612, 12, 11, None, None),
    ("SCA6_CACNA1A", 11, 13, 20, 19, None, None),
]

# A few generic genome-wide catalog regions (id matches repeat_catalog.hg38.bed's ID=),
# to exercise Repeat Catalog matching outside the STRchive panel.
GENERIC_CALLS = [
    # (chrom, start, stop, motif, a1_units, a2_units, a1_depth, a2_depth)
    ("chr1", 44835, 44867, "AAAT", 8, 9, 15, 14),
    ("chr1", 151101, 151105, "CT", 2, 2, 18, 17),
    ("chr1", 832736, 832781, "GTTTT", 9, 11, 13, 12),
]

HEADER = f"""##fileformat=VCFv4.3
##source=synthetic-trgt-v5.1.0-fixture
##INFO=<ID=TRID,Number=1,Type=String,Description="Identifier of the repeat region from the input catalog">
##INFO=<ID=END,Number=1,Type=Integer,Description="Ending position of the repeat region">
##INFO=<ID=MOTIFS,Number=.,Type=String,Description="Motif units observed at the locus">
##INFO=<ID=STRUC,Number=1,Type=String,Description="Expected locus structure (deprecated)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AL,Number=.,Type=Integer,Description="Allele length in base pairs">
##FORMAT=<ID=ALLR,Number=.,Type=String,Description="Length confidence interval per allele">
##FORMAT=<ID=SD,Number=.,Type=Integer,Description="Number of spanning reads supporting each allele">
##FORMAT=<ID=MC,Number=.,Type=String,Description="Motif copy count(s) per allele">
##FORMAT=<ID=MS,Number=.,Type=String,Description="Run-length segmentation string per allele">
##FORMAT=<ID=AP,Number=.,Type=Float,Description="Allele purity">
##FORMAT=<ID=AM,Number=.,Type=Float,Description="Mean 5mCpG methylation level per allele">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{SAMPLE_NAME}"""


def fmt(v) -> str:
    return "." if v is None else (f"{v:.3f}" if isinstance(v, float) else str(v))


def build_record(chrom, start, stop, motif, trid, a1_units, a2_units, a1_depth, a2_depth,
                  a1_meth=None, a2_meth=None):
    """Build one VCF record matching TRGT's real output format.

    Verified against an actual trgt v5.1.0 genotype call on HG002: AL/ALLR/SD/MC/AP/AM
    are NOT Number=R (one slot per declared REF+ALT allele) -- they're always exactly
    one slot per haplotype (2 for diploid), ordered by ascending *original* allele
    index of the two haplotypes actually called. A homozygous call still gets two
    independent per-haplotype-copy values (e.g. different real SD read counts), not
    one shared/deduplicated value.
    """
    motif_len = len(motif)
    ref_units = max(1, round((stop - start) / motif_len))
    ref_bp = ref_units * motif_len
    a1_bp, a2_bp = a1_units * motif_len, a2_units * motif_len

    # Assign each haplotype an allele index: 0 if it matches the reference length,
    # otherwise a new ALT index in order of first appearance.
    alt_bp: list[int] = []

    def allele_index(bp):
        if bp == ref_bp:
            return 0
        for i, existing in enumerate(alt_bp):
            if existing == bp:
                return i + 1
        alt_bp.append(bp)
        return len(alt_bp)

    idx1 = allele_index(a1_bp)
    idx2 = allele_index(a2_bp)
    alt_str = ",".join(f"<STR{i}n{bp}>" for i, bp in enumerate(alt_bp, start=1)) or "."

    # FORMAT arrays: exactly 2 slots, ordered by ascending original allele index.
    # hap1 -> slot 0 on ties (matches observed real behavior for homozygous calls).
    hap1_first = idx1 <= idx2
    slot0 = (a1_bp, a1_units, a1_depth, a1_meth) if hap1_first else (a2_bp, a2_units, a2_depth, a2_meth)
    slot1 = (a2_bp, a2_units, a2_depth, a2_meth) if hap1_first else (a1_bp, a1_units, a1_depth, a1_meth)

    al_full = [slot0[0], slot1[0]]
    allr_full = [(slot0[0] - 1, slot0[0] + 1), (slot1[0] - 1, slot1[0] + 1)]
    sd_full = [slot0[2], slot1[2]]
    mc_full = [slot0[1], slot1[1]]
    ap_full = [0.97, 0.97]
    am_full = [slot0[3], slot1[3]]

    gt = f"{idx1}|{idx2}"
    info = f"TRID={trid};END={stop};MOTIFS={motif};STRUC=({motif})n"
    sample_fmt = ":".join([
        gt,
        ",".join(str(v) for v in al_full),
        ",".join(f"{lo}-{hi}" for lo, hi in allr_full),
        ",".join(str(v) for v in sd_full),
        ",".join(str(v) for v in mc_full),
        ",".join(["."] * 2),
        ",".join(fmt(v) for v in ap_full),
        ",".join(fmt(v) for v in am_full),
        "1",
    ])
    ref_seq = "N"
    return "\t".join([
        chrom, str(start), trid, ref_seq, alt_str, ".", "PASS", info,
        "GT:AL:ALLR:SD:MC:MS:AP:AM:PS", sample_fmt,
    ])


def main():
    with open(STRCHIVE_JSON, encoding="utf-8") as f:
        loci = {d["id"]: d for d in json.load(f)}

    records = []
    for locus_id, a1, a2, d1, d2, m1, m2 in NARRATIVE_CALLS:
        d = loci[locus_id]
        motif = (d.get("pathogenic_motif_reference_orientation") or ["N"])[0]
        records.append(
            (d["chrom"], d["start_hg38"],
             build_record(d["chrom"], d["start_hg38"], d["stop_hg38"], motif, locus_id,
                           a1, a2, d1, d2, m1, m2))
        )

    for i, (chrom, start, stop, motif, a1, a2, d1, d2) in enumerate(GENERIC_CALLS):
        trid = f"{chrom}_{start}_{stop}"
        records.append((chrom, start, build_record(chrom, start, stop, motif, trid, a1, a2, d1, d2)))

    # VCF must be position-sorted per chromosome for our index's binary search to be valid;
    # sort by (chrom, start) using natural chromosome order.
    chrom_order = {f"chr{c}": i for i, c in enumerate(list(range(1, 23)) + ["X", "Y", "M"])}
    records.sort(key=lambda r: (chrom_order.get(r[0], 99), r[1]))

    print(HEADER)
    for _, _, line in records:
        print(line)


if __name__ == "__main__":
    main()
