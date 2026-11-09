"""FSHD2 research-evidence tiers: a permissive D4Z4 allele, D4Z4 hypomethylation, and a variant in an FSHD2 gene.

Pure functions over already-shaped data (no filesystem), so they are unit-testable with small synthetic inputs.

FSHD2 needs a DUX4-permissive D4Z4 allele (chr4, qAIntactPolyA, not 4A166) *and* a second hit that lowers D4Z4
methylation: a pathogenic variant in SMCHD1, DNMT3B or LRIF1. The tiers count how many of the three findings a sample
shows. They are research evidence, not a diagnosis: Kivvi defines no FSHD2 rule, the hypomethylation cut-off is this
project's (`parse.HYPOMETHYLATION_CUTOFF`), and the variants are supplied by the user rather than called here.

The permissive allele is a gate. A sample without one cannot express DUX4 from D4Z4, so no tier is given however
low its methylation or however damaging its variant; the tier then reads "not supported". "No permissive allele" means
none among the alleles Kivvi assembled: with fewer than four resolved, or an allele it could not place on a
chromosome, a permissive allele may be missing from the data, and the notes say so.
"""

from __future__ import annotations

import re
from typing import Optional

from .parse import EXPECTED_D4Z4_ALLELES

# The genes whose pathogenic variants cause FSHD2 (SMCHD1 most often; DNMT3B and LRIF1 rarely).
FSHD2_GENES = ("SMCHD1", "DNMT3B", "LRIF1")

# Consequence classes that remove or truncate the protein, which count as a hit without further evidence.
LOSS_OF_FUNCTION = ("frameshift", "nonsense", "splice_site")
CONSEQUENCE_LABELS = {
    "frameshift": "Frameshift",
    "nonsense": "Nonsense",
    "splice_site": "Canonical splice site",
    "missense": "Missense",
    "other": "Other / unclassified",
}
# ClinVar classifications that count as a hit for a variant that is not loss-of-function.
CLINVAR_COUNTS = ("pathogenic", "likely pathogenic", "pathogenic/likely pathogenic")

TIER_LABELS = {
    "all_three": "Tier 1 · all three findings",
    "two_of_three": "Tier 2 · permissive allele and one other finding",
    "permissive_only": "Tier 3 · permissive allele only",
    "no_permissive": "Not supported · no permissive allele found",
}

_SPLICE = re.compile(r"c\.-?\*?\d+[+-][12](?!\d)")
_FRAMESHIFT_P = re.compile(r"fs")
_NONSENSE_P = re.compile(r"p\.\(?[A-Z][a-z]{0,2}\d+(\*|X|Ter)\)?$")
_MISSENSE_P = re.compile(r"p\.\(?[A-Z][a-z]{0,2}\d+[A-Z][a-z]{0,2}\)?$")
_INDEL_C = re.compile(r"c\.(\d+)(?:_(\d+))?(del|dup)")


def consequence(hgvs_c: Optional[str], hgvs_p: Optional[str]) -> str:
    """The consequence class of a variant from its HGVS notation: frameshift, nonsense, splice_site (canonical +-1/2),
    missense, or other. The protein change decides when there is one; the coding change is the fallback."""
    c = (hgvs_c or "").replace(" ", "")
    p = (hgvs_p or "").replace(" ", "")
    if _SPLICE.search(c):
        return "splice_site"
    if p:
        if _FRAMESHIFT_P.search(p):
            return "frameshift"
        if _NONSENSE_P.search(p):
            return "nonsense"
        if _MISSENSE_P.search(p):
            return "missense"
        return "other"
    indel = _INDEL_C.search(c)
    if indel:
        start, end = int(indel.group(1)), int(indel.group(2) or indel.group(1))
        return "frameshift" if (end - start + 1) % 3 else "other"
    return "other"


def annotate_variant(variant: dict) -> dict:
    """A supplied variant with its consequence class and whether it counts as a hit.

    A variant counts when it is in an FSHD2 gene and is either loss-of-function (frameshift, nonsense, canonical
    splice) or classified pathogenic / likely pathogenic in ClinVar. A missense variant of uncertain significance
    is shown but does not count: it would overstate the evidence.
    """
    gene = (variant.get("gene") or "").upper()
    cons = consequence(variant.get("hgvs_c"), variant.get("hgvs_p"))
    clinvar = (variant.get("clinvar_classification") or "").strip().lower()
    in_gene = gene in FSHD2_GENES
    if not in_gene:
        counts, why = False, f"{gene or 'The gene'} is not an FSHD2 gene"
    elif cons in LOSS_OF_FUNCTION:
        counts, why = True, f"{CONSEQUENCE_LABELS[cons].lower()} variant"
    elif clinvar in CLINVAR_COUNTS:
        counts, why = True, f"ClinVar: {variant['clinvar_classification']}"
    else:
        counts, why = False, "not loss-of-function and not ClinVar pathogenic / likely pathogenic"
    return {
        **variant,
        "gene": gene,
        "consequence": cons,
        "consequence_label": CONSEQUENCE_LABELS[cons],
        "in_fshd2_gene": in_gene,
        "counts": counts,
        "counts_reason": why,
    }


def _variant_criterion(variants: Optional[list[dict]]) -> dict:
    if variants is None:
        return {
            "state": "not_assessed",
            "met": False,
            "summary": "Not assessed: no variant calls were supplied for this sample.",
            "variants": [],
        }
    annotated = [annotate_variant(v) for v in variants]
    counting = [v for v in annotated if v["counts"]]
    if counting:
        names = ", ".join(f"{v['gene']} {v['hgvs_c']}" for v in counting)
        return {"state": "met", "met": True, "summary": f"{len(counting)} reported: {names}.", "variants": annotated}
    if annotated:
        return {
            "state": "uncertain",
            "met": False,
            "summary": "Reported, but none counts as a hit (not loss-of-function, not ClinVar pathogenic / likely pathogenic).",
            "variants": annotated,
        }
    return {"state": "none_reported", "met": False, "summary": "None reported in SMCHD1, DNMT3B or LRIF1.", "variants": []}


def evidence(fshd2: dict, variants: Optional[list[dict]]) -> dict:
    """The tier for one sample.

    `fshd2` is the block from `parse.shape_fshd2`; `variants` is the list supplied for the sample, an empty list when
    the sample was checked and none was found, and None when it was never assessed.
    """
    summary = fshd2["summary"]
    alleles = fshd2["alleles"]
    permissive = [a for a in alleles if a["permissive"]]
    hypo = fshd2["hypomethylation"]["positive"]
    variant = _variant_criterion(variants)

    criteria = {
        "permissive_allele": {
            "met": bool(permissive),
            "summary": (
                f"{len(permissive)} permissive allele{'s' if len(permissive) != 1 else ''}."
                if permissive
                else "No permissive allele (chr4, qAIntactPolyA, not 4A166)."
            ),
        },
        "hypomethylation": {
            "met": bool(hypo),
            "summary": (
                "Kivvi reported no sample methylation."
                if hypo is None
                else f"Sample methylation is {'below' if hypo else 'at or above'} the {round(fshd2['hypomethylation']['cutoff'] * 100)}% cut-off."
            ),
        },
        "variant": variant,
    }

    others = int(bool(hypo)) + int(variant["met"])
    if not permissive:
        status = "no_permissive"
    elif others == 2:
        status = "all_three"
    elif others == 1:
        status = "two_of_three"
    else:
        status = "permissive_only"
    tier = {"all_three": 1, "two_of_three": 2, "permissive_only": 3}.get(status)

    notes = []
    if any(a["size_is_lower_bound"] for a in permissive):
        notes.append("A permissive allele is only partly assembled, so its size is a lower bound; whether it is permissive rests on its haplotype.")
    unplaced = [a for a in alleles if a["distal_haplotype"] == "qAIntactPolyA" and a["chromosome"] == "unknown"]
    if not permissive and len(alleles) < EXPECTED_D4Z4_ALLELES:
        notes.append(
            f"Kivvi resolved {len(alleles)} of {EXPECTED_D4Z4_ALLELES} expected D4Z4 alleles, so a permissive allele may be among those not assembled."
        )
    if unplaced and not permissive:
        notes.append("An allele with an intact qA haplotype could not be placed on a chromosome; if it is chr4 it would be permissive.")
    if variant["state"] == "uncertain":
        notes.append("A reported variant is shown but does not count towards the tier.")
    if status == "no_permissive" and (hypo or variant["met"]):
        notes.append("Hypomethylation or a variant is present, but without a permissive allele the evidence does not support FSHD2.")
    if hypo is None:
        notes.append("Hypomethylation could not be scored.")
    # A pending variant can still move the tier, so say so; it cannot when there is no permissive allele.
    pending = variant["state"] == "not_assessed" and bool(permissive)
    if pending:
        notes.append("The variant criterion has not been assessed, so the tier may rise once variant calls are supplied.")

    return {
        "status": status,
        "tier": tier,
        "label": TIER_LABELS[status],
        "criteria": criteria,
        "variants": variant["variants"],
        "pending_variant": pending,
        "permissive_alleles": summary["permissive_alleles"],
        "notes": notes,
    }
