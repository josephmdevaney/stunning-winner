#!/usr/bin/env python3
"""Generate synthetic TRGT v5.1.0 VCFs for the Carrier screen fixture panel.

Carrier mode's config (`config/modes.json`) has always pointed `vcf_dir` at
`fixtures/carrier_panel` and `panel_bed` at the real, hand-curated
`fixtures/panels/carrier_screen.repeat_definitions.hg38.bed` (5 STRchive repeat loci
carried on PacBio's own carrier-screening panel) plus
`carrier_screen.non_repeat_targets.json` (8 genes the panel also screens via
non-repeat callers this pipeline doesn't run) -- but no VCF fixtures existed yet, so
the mode always showed "no data". Unlike `gv_puretarget_panel` (real vendor-provided
`trgt` output), no real carrier-panel BAMs were available here, so -- exactly like
`fixtures/synthetic_trgt.vcf` -- this builds spec-accurate synthetic `trgt` v5.1.0
records by hand for just the 5 repeat loci this panel actually calls, one
single-sample VCF per fixture sample, covering a realistic carrier-screening spread
(negative, intermediate/premutation carriers, full-expansion carriers, and a no-call).

Run with:

    python backend/scripts/generate_carrier_panel_vcfs.py
"""

from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
STRCHIVE_JSON = REPO_ROOT / "project" / "uploads" / "STRchive-loci.json"
OUT_DIR = REPO_ROOT / "backend" / "fixtures" / "carrier_panel"

LOCUS_IDS = ["FRDA_FXN", "PRTS_ARX", "EIEE1_ARX", "FXS_FMR1", "FRAXE_AFF2"]

# Per sample, per locus: (a1_units, a2_units, a1_depth, a2_depth) or None for no-call.
# Values are chosen against each locus's real STRchive benign/intermediate/pathogenic
# thresholds to exercise the full range a carrier panel is meant to catch.
SAMPLES: dict[str, dict[str, tuple[int, int, int, int] | None]] = {
    # Clean negative screen: every locus within its benign range.
    "CR-01": {
        "FRDA_FXN": (9, 11, 22, 21),
        "PRTS_ARX": (12, 12, 19, 18),
        "EIEE1_ARX": (11, 13, 20, 19),
        "FXS_FMR1": (29, 30, 24, 23),
        "FRAXE_AFF2": (18, 22, 21, 20),
    },
    # FRDA_FXN carrier: one benign allele, one past the pathogenic threshold (>=56) --
    # the classic asymptomatic heterozygous carrier of an autosomal-recessive repeat
    # expansion.
    "CR-02": {
        "FRDA_FXN": (13, 92, 18, 15),
        "PRTS_ARX": (12, 12, 20, 19),
        "EIEE1_ARX": (12, 14, 21, 20),
        "FXS_FMR1": (28, 31, 23, 22),
        "FRAXE_AFF2": (20, 24, 22, 21),
    },
    # FXS_FMR1 premutation carrier: one allele in the 45-200 intermediate/premutation
    # range -- unaffected but at risk of expansion in offspring (the core FXS carrier
    # scenario).
    "CR-03": {
        "FRDA_FXN": (10, 15, 20, 19),
        "PRTS_ARX": (12, 12, 18, 17),
        "EIEE1_ARX": (10, 12, 19, 18),
        "FXS_FMR1": (29, 82, 21, 12),
        "FRAXE_AFF2": (19, 23, 20, 19),
    },
    # FRAXE_AFF2 premutation carrier (40-200 range).
    "CR-04": {
        "FRDA_FXN": (11, 14, 19, 18),
        "PRTS_ARX": (12, 12, 19, 18),
        "EIEE1_ARX": (11, 11, 20, 20),
        "FXS_FMR1": (27, 30, 22, 21),
        "FRAXE_AFF2": (21, 96, 20, 11),
    },
    # PRTS_ARX carrier: heterozygous for the pathogenic 20-repeat polyalanine
    # expansion (X-linked recessive -- unaffected female carrier).
    "CR-05": {
        "FRDA_FXN": (8, 12, 21, 20),
        "PRTS_ARX": (12, 20, 17, 14),
        "EIEE1_ARX": (12, 13, 18, 17),
        "FXS_FMR1": (30, 30, 20, 19),
        "FRAXE_AFF2": (17, 25, 21, 20),
    },
    # EIEE1_ARX carrier: heterozygous for a 17-27 pathogenic-range allele.
    "CR-06": {
        "FRDA_FXN": (9, 13, 20, 19),
        "PRTS_ARX": (12, 12, 21, 20),
        "EIEE1_ARX": (13, 21, 16, 13),
        "FXS_FMR1": (29, 29, 23, 22),
        "FRAXE_AFF2": (16, 20, 19, 18),
    },
    # FRDA_FXN intermediate allele (34-55) only -- borderline finding warranting
    # follow-up rather than a clear-cut carrier call.
    "CR-07": {
        "FRDA_FXN": (14, 40, 17, 16),
        "PRTS_ARX": (12, 12, 18, 17),
        "EIEE1_ARX": (12, 12, 19, 18),
        "FXS_FMR1": (31, 33, 21, 20),
        "FRAXE_AFF2": (19, 21, 20, 19),
    },
    # Second clean negative screen.
    "CR-08": {
        "FRDA_FXN": (10, 12, 23, 22),
        "PRTS_ARX": (12, 12, 20, 19),
        "EIEE1_ARX": (10, 11, 21, 20),
        "FXS_FMR1": (28, 29, 22, 21),
        "FRAXE_AFF2": (18, 20, 21, 20),
    },
    # FXS_FMR1 full mutation (>=201): the panel also catches true affected cases, not
    # just carriers -- both alleles expanded here (homozygous full mutation).
    "CR-09": {
        "FRDA_FXN": (9, 14, 18, 17),
        "PRTS_ARX": (12, 12, 17, 16),
        "EIEE1_ARX": (11, 13, 18, 17),
        "FXS_FMR1": (240, 260, 9, 8),
        "FRAXE_AFF2": (20, 24, 19, 18),
    },
    # No-call at one locus (insufficient spanning reads) alongside an otherwise clean
    # screen -- QC should flag FRAXE_AFF2 here rather than silently showing benign.
    "CR-10": {
        "FRDA_FXN": (11, 13, 19, 18),
        "PRTS_ARX": (12, 12, 19, 18),
        "EIEE1_ARX": (12, 14, 20, 19),
        "FXS_FMR1": (29, 31, 21, 20),
        "FRAXE_AFF2": None,
    },
}

HEADER_TEMPLATE = """##fileformat=VCFv4.3
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
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}"""


def fmt(v) -> str:
    return "." if v is None else (f"{v:.3f}" if isinstance(v, float) else str(v))


def build_record(chrom, start, stop, motif, trid, call):
    """Build one TRGT-format VCF record; `call` is None for a no-call site."""
    motif_len = len(motif)
    ref_units = max(1, round((stop - start) / motif_len))
    ref_bp = ref_units * motif_len

    if call is None:
        info = f"TRID={trid};END={stop};MOTIFS={motif};STRUC=({motif})n"
        sample_fmt = ":".join(["./.", ".", ".", ".", ".", ".", ".", ".", "."])
        return "\t".join(
            [chrom, str(start), trid, "N", ".", ".", "PASS", info, "GT:AL:ALLR:SD:MC:MS:AP:AM:PS", sample_fmt]
        )

    a1_units, a2_units, a1_depth, a2_depth = call
    a1_bp, a2_bp = a1_units * motif_len, a2_units * motif_len

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

    hap1_first = idx1 <= idx2
    slot0 = (a1_bp, a1_units, a1_depth) if hap1_first else (a2_bp, a2_units, a2_depth)
    slot1 = (a2_bp, a2_units, a2_depth) if hap1_first else (a1_bp, a1_units, a1_depth)

    al_full = [slot0[0], slot1[0]]
    allr_full = [(slot0[0] - 1, slot0[0] + 1), (slot1[0] - 1, slot1[0] + 1)]
    sd_full = [slot0[2], slot1[2]]
    mc_full = [slot0[1], slot1[1]]
    ap_full = [0.97, 0.97]
    am_full = [None, None]

    gt = f"{idx1}|{idx2}"
    info = f"TRID={trid};END={stop};MOTIFS={motif};STRUC=({motif})n"
    sample_fmt = ":".join(
        [
            gt,
            ",".join(str(v) for v in al_full),
            ",".join(f"{lo}-{hi}" for lo, hi in allr_full),
            ",".join(str(v) for v in sd_full),
            ",".join(str(v) for v in mc_full),
            ",".join(["."] * 2),
            ",".join(fmt(v) for v in ap_full),
            ",".join(fmt(v) for v in am_full),
            "1",
        ]
    )
    return "\t".join([chrom, str(start), trid, "N", alt_str, ".", "PASS", info, "GT:AL:ALLR:SD:MC:MS:AP:AM:PS", sample_fmt])


def main():
    with open(STRCHIVE_JSON, encoding="utf-8") as f:
        loci = {d["id"]: d for d in json.load(f)}

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    for sample_name, calls in SAMPLES.items():
        records = []
        for locus_id in LOCUS_IDS:
            d = loci[locus_id]
            motif = (d.get("pathogenic_motif_reference_orientation") or ["N"])[0]
            line = build_record(
                d["chrom"], d["start_hg38"], d["stop_hg38"], motif, locus_id, calls[locus_id]
            )
            records.append((d["chrom"], d["start_hg38"], line))

        chrom_order = {f"chr{c}": i for i, c in enumerate(list(range(1, 23)) + ["X", "Y", "M"])}
        records.sort(key=lambda r: (chrom_order.get(r[0], 99), r[1]))

        out_path = OUT_DIR / f"{sample_name}.trgt.vcf"
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(HEADER_TEMPLATE.format(sample=sample_name) + "\n")
            for _, _, line in records:
                f.write(line + "\n")
        print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
