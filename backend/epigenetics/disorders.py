"""Imprinting-disorder patterns from per-region calls.

Each disorder lists the DMRs whose methylation change is its usual epigenetic signature and the direction of
that change, from two open-access consensus papers:
  * Multi-locus imprinting disturbance (MLID): interim joint statement, Clin Epigenet 2024 (PMC11295890)
  * First step towards a consensus strategy for multi-locus diagnostic testing, Clin Epigenet 2022 (PMC9641836)

A methylation pattern is *consistent with* a disorder; it never establishes one. Methylation alone cannot
tell an imprinting defect from a deletion or uniparental disomy (which need copy-number, ROH or trio data),
and mosaicism means a blood-derived result can miss a real disturbance. Every disorder summary carries those
limits, and results are for research use.
"""

from __future__ import annotations

from typing import Optional

from .analysis import ABERRANT, STATUS_GAIN, STATUS_IMPRINTED, STATUS_LOSS

# trigger: (dmr id, direction, primary?). A primary trigger is the DMR the literature treats as the main
# diagnostic locus for that disorder; supportive ones add weight but do not establish a pattern alone.
DISORDER_RULES: list[dict] = [
    {
        "id": "pws", "name": "Prader-Willi syndrome", "locus": "15q11.2",
        "triggers": [("snurf-tss-dmr", STATUS_GAIN, True)],
        "meaning": "Both alleles carry the maternal (methylated) pattern at the SNURF/SNRPN imprinting centre.",
        "mechanisms": "paternal 15q11-q13 deletion, maternal UPD15, or an imprinting defect",
    },
    {
        "id": "as", "name": "Angelman syndrome", "locus": "15q11.2",
        "triggers": [("snurf-tss-dmr", STATUS_LOSS, True)],
        "meaning": "Both alleles carry the paternal (unmethylated) pattern at the SNURF/SNRPN imprinting centre.",
        "mechanisms": "maternal 15q11-q13 deletion, paternal UPD15, or an imprinting defect (UBE3A variants are not epigenetic)",
    },
    {
        "id": "bws", "name": "Beckwith-Wiedemann syndrome", "locus": "11p15.5",
        "triggers": [("kcnq1ot1-tss-dmr", STATUS_LOSS, True), ("h19-igf2-ig-dmr-tss-dmr", STATUS_GAIN, True)],
        "meaning": "Loss of methylation at the KCNQ1OT1 TSS-DMR (IC2), or gain at the H19/IGF2 IG-DMR (IC1).",
        "mechanisms": "isolated IC2 loss (sometimes part of a multi-locus disturbance), IC1 gain (can accompany 11p15 CNVs), "
                      "or mosaic paternal UPD11 (IC1 gain with IC2 loss); CDKN1C variants are not epigenetic",
    },
    {
        "id": "srs", "name": "Silver-Russell syndrome", "locus": "11p15.5 / 7",
        "triggers": [("h19-igf2-ig-dmr-tss-dmr", STATUS_LOSS, True), ("mest-alt-tss-dmr", STATUS_GAIN, False),
                     ("grb10-alt-tss-dmr", STATUS_GAIN, False)],
        "meaning": "Loss of methylation at the H19/IGF2 IG-DMR (IC1), or gain at MEST and GRB10 (the maternal UPD7 pattern).",
        "mechanisms": "IC1 loss (the most frequent cause), maternal UPD7 (MEST/GRB10 gain), or 11p15 copy-number changes",
    },
    {
        "id": "tndm", "name": "Transient neonatal diabetes mellitus (6q24)", "locus": "6q24.2",
        "triggers": [("plagl1-alt-tss-dmr", STATUS_LOSS, True)],
        "meaning": "Loss of methylation at the PLAGL1 alt-TSS-DMR.",
        "mechanisms": "paternal UPD6, an imprinting defect (sometimes ZFP57-related multi-locus), or a paternal 6q24 duplication, "
                      "which does not change methylation",
    },
    {
        "id": "php1b", "name": "Pseudohypoparathyroidism type 1B (GNAS)", "locus": "20q13.32",
        "triggers": [("gnas-a-b-tss-dmr", STATUS_LOSS, True), ("gnas-xl-ex1-dmr", STATUS_LOSS, False),
                     ("gnas-as1-tss-dmr", STATUS_LOSS, False), ("gnas-nesp-tss-dmr", STATUS_GAIN, False)],
        "meaning": "Loss of maternal methylation at GNAS A/B, often with loss at XL/AS and gain at NESP.",
        "mechanisms": "sporadic imprinting defect, STX16 deletion (isolated A/B loss), or paternal UPD20",
    },
    {
        "id": "ts14", "name": "Temple syndrome", "locus": "14q32.2",
        "triggers": [("meg3-tss-dmr", STATUS_LOSS, True), ("meg3-dlk1-ig-dmr", STATUS_LOSS, False)],
        "meaning": "Loss of methylation at the MEG3 TSS-DMR and/or the MEG3/DLK1 IG-DMR (the maternal pattern on both alleles).",
        "mechanisms": "maternal UPD14, paternal 14q32 deletion, or an imprinting defect",
    },
    {
        "id": "kos14", "name": "Kagami-Ogata syndrome", "locus": "14q32.2",
        "triggers": [("meg3-tss-dmr", STATUS_GAIN, True), ("meg3-dlk1-ig-dmr", STATUS_GAIN, False)],
        "meaning": "Gain of methylation at the MEG3 TSS-DMR and/or the MEG3/DLK1 IG-DMR (the paternal pattern on both alleles).",
        "mechanisms": "paternal UPD14, maternal 14q32 deletion, or an imprinting defect",
    },
]

CAVEATS = (
    "A methylation pattern is consistent with a disorder, not a diagnosis. Methylation alone cannot separate an "
    "imprinting defect from a deletion or uniparental disomy (use copy number, ROH or trio data), and mosaicism "
    "can hide a real disturbance in blood. Confirm with a clinically validated test."
)


def evaluate_disorders(results: dict[str, dict]) -> list[dict]:
    """`results` maps a DMR id to its region_result. Returns one summary per disorder."""
    out = []
    for rule in DISORDER_RULES:
        triggers, primary_assessed, primary_hit, supportive_hit = [], 0, False, False
        for dmr_id, direction, primary in rule["triggers"]:
            r: Optional[dict] = results.get(dmr_id)
            status = r["status"] if r else "not_assessed"
            callable_ = r is not None and status not in ("low_coverage", "no_data", "not_assessed")
            if primary and callable_:
                primary_assessed += 1
            hit = callable_ and status == direction
            if hit and primary:
                primary_hit = True
            if hit and not primary:
                supportive_hit = True
            triggers.append({
                "dmr": dmr_id, "name": r["name"] if r else dmr_id, "expected": direction, "primary": primary,
                "observed": status, "hit": bool(hit),
                "methylation": (r["combined"]["methylation"] if r and r["combined"] else None),
            })
        n_primary = sum(1 for t in rule["triggers"] if t[2])
        if primary_hit:
            status = "pattern_present"
        elif supportive_hit:
            status = "supportive_only"
        elif primary_assessed == n_primary:
            status = "not_observed"
        else:
            status = "insufficient_data"
        out.append({
            "id": rule["id"], "name": rule["name"], "locus": rule["locus"], "status": status,
            "meaning": rule["meaning"], "mechanisms": rule["mechanisms"], "triggers": triggers,
        })
    return out


def evaluate_mlid(results: dict[str, dict], clinical_only: bool = True) -> dict:
    """Multi-locus imprinting disturbance: aberrant calls on two or more chromosomes.

    Counting chromosomes rather than DMRs keeps a single event that moves several DMRs on one chromosome
    (a UPD, or IC1 plus IC2 in paternal UPD11) from reading as multi-locus. The 2024 joint statement treats
    comprehensive MLID analysis, which includes non-clinical DMRs, as research; `clinical_only` restricts to
    the clinically associated (CA) DMRs."""
    considered = [r for r in results.values() if r["category"] == "CA" or (not clinical_only and r["category"] == "NC")]
    aberrant = [r for r in considered if r["status"] in ABERRANT]
    chromosomes = sorted({r["chrom"] for r in aberrant})
    assessed = [r for r in considered if r["status"] not in ("low_coverage", "no_data")]
    return {
        "scope": "clinically associated DMRs" if clinical_only else "clinically associated and non-clinical DMRs (research)",
        "assessed": len(assessed), "considered": len(considered),
        "aberrant": [{"id": r["id"], "name": r["name"], "chrom": r["chrom"], "status": r["status"]} for r in aberrant],
        "chromosomes": chromosomes,
        "possible": len(chromosomes) >= 2,
    }
