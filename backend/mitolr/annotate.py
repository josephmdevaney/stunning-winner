"""Annotate real mtDNA variant calls against the rCRS feature map and the
curated pathogenic panel.

A VCF alone gives position/ref/alt/depth/heteroplasmy but not functional
consequence (synonymous vs missense needs the codon, hence the reference
sequence) or population frequency. This module annotates what is honestly
derivable — locus, functional class, and any known pathogenic association — and
leaves the rest marked as not evaluated rather than fabricating it.
"""

from __future__ import annotations

from typing import Optional

from . import gnomad as gnomad_mod
from .reference import (
    CONSEQUENCE_BY_TYPE,
    PATHOGENIC_PANEL,
    STATUS_COLORS,
    feature_for_pos,
)

# Known pathogenic variants keyed by (pos, ref, alt) for exact-call matching.
_PATHO_BY_KEY: dict[tuple[int, str, str], dict] = {}
for _v in PATHOGENIC_PANEL:
    _m = None
    _pos_str = _v["position"]
    # Parse "m.<pos><ref>><alt>" for SNVs; deletions are matched by pos only.
    import re as _re

    _snv = _re.match(r"m\.(\d+)([ACGT])>([ACGT])", _pos_str)
    if _snv:
        _PATHO_BY_KEY[(int(_snv.group(1)), _snv.group(2), _snv.group(3))] = _v


def match_pathogenic(pos: int, ref: str, alt: str) -> Optional[dict]:
    return _PATHO_BY_KEY.get((pos, ref, alt))


def locus_and_consequence(pos: int) -> tuple[str, str]:
    """Return (locus_name, consequence_class) for a genomic position."""
    feat = feature_for_pos(pos)
    if feat is None:
        return ("Intergenic", "Noncoding")
    name = "D-loop" if feat["name"] == "MT-CR" else feat["name"]
    return (name, CONSEQUENCE_BY_TYPE.get(feat["type"], "Noncoding"))


def variant_notation(pos: int, ref: str, alt: str) -> tuple[str, str, bool]:
    """Build (display position, change, is_large_deletion) in rCRS `m.` notation.

    Handles SNVs and the small insertions/deletions MitoRSaw emits (anchored,
    e.g. A>AC insertion, CN>C deletion)."""
    if len(ref) == 1 and len(alt) == 1:
        return (f"m.{pos}{ref}>{alt}", f"{ref}>{alt}", False)
    if len(alt) > len(ref):  # insertion (anchored on ref[0])
        inserted = alt[len(ref):]
        return (f"m.{pos}_{pos + 1}ins{inserted}", f"ins{inserted}", False)
    if len(ref) > len(alt):  # deletion (anchored on ref[0])
        del_start = pos + len(alt)
        del_stop = pos + len(ref) - 1
        span = del_stop - del_start + 1
        label = f"m.{del_start}del" if span == 1 else f"m.{del_start}_{del_stop}del"
        # Treat only large (>100 bp) deletions as the circular-arc "isDeletion" kind.
        return (label, "del", span > 100)
    # Same-length MNV / substitution.
    return (f"m.{pos}{ref}>{alt}", f"{ref}>{alt}", False)


def classify(
    pos: int,
    ref: str,
    alt: str,
    het: float,
    gnomad_hit: Optional[gnomad_mod.GnomadHit] = None,
    gnomad_checked: bool = False,
) -> dict:
    """Classify a called variant into the display status taxonomy.

    Known-pathogenic calls use the heteroplasmy bands, unaffected by
    population frequency. Novel (non-pathogenic-panel) calls are reported by
    heteroplasmy state; when ``gnomad_checked`` is True (real VCF-derived
    calls only -- see vcf.py), a homoplasmic novel call is further split into
    Common vs Rare using its real gnomAD v3.1 AF_hom, rather than assuming
    "homoplasmic implies common." ``gnomad_hit`` is None either because the
    variant was checked and not found in the callset (gnomad_checked=True) or
    because gnomAD wasn't consulted at all (gnomad_checked=False, e.g. the
    synthetic demo sample) -- the note text distinguishes the two."""
    patho = match_pathogenic(pos, ref, alt)
    if patho is not None:
        if het >= 30:
            status, note = (
                "Pathogenic",
                f"Elevated heteroplasmy ({het}%) above the clinical significance threshold "
                f"(>30%); consistent with {patho['disease']}.",
            )
        elif het >= 5:
            status, note = (
                "Low-level heteroplasmy",
                f"Known pathogenic variant at low-level heteroplasmy ({het}%); below typical "
                f"clinical significance threshold but may warrant follow-up.",
            )
        else:
            status, note = (
                "Not detected",
                f"Known pathogenic position below reporting threshold (<5%).",
            )
        return {
            "status": status,
            "statusColor": STATUS_COLORS[status],
            "note": note,
            "disease": patho["disease"],
            "evidence": patho["evidence"],
            "inheritance": patho["inheritance"],
            "origin": "pathogenic",
        }

    # No known disease association.
    if het < 5:
        status = "Not detected"
        note = f"Sub-threshold call ({het}%); below the 5% reporting threshold."
    elif het < 95:
        status = "VUS"
        if gnomad_hit is not None:
            note = (
                f"Heteroplasmic call ({het}%) with no known disease association — variant of "
                f"uncertain significance. Also observed heteroplasmic in {gnomad_mod.SOURCE_LABEL} "
                f"(AF_het={gnomad_hit['af_het']:.2%}, n={gnomad_hit['an']:,}); population frequency "
                f"alone doesn't resolve significance for a heteroplasmic call."
            )
        elif gnomad_checked:
            note = (
                f"Heteroplasmic call ({het}%) with no known disease association — variant of "
                f"uncertain significance. Not found in {gnomad_mod.SOURCE_LABEL} "
                f"(n={gnomad_mod.COHORT_SIZE:,})."
            )
        else:
            note = (
                f"Heteroplasmic call ({het}%) with no known disease association — variant of "
                f"uncertain significance. Population frequency not evaluated."
            )
    elif gnomad_hit is not None:
        af = gnomad_hit["af_hom"]
        if af >= 0.01:
            status = "Common variant"
            note = (
                f"Homoplasmic call vs rCRS ({het}%); no known disease association. Common in "
                f"{gnomad_mod.SOURCE_LABEL} (AF_hom={af:.2%}, n={gnomad_hit['an']:,}) — likely "
                f"haplogroup-defining or benign."
            )
        else:
            status = "Rare variant"
            note = (
                f"Homoplasmic call vs rCRS ({het}%); no known disease association. Rare in "
                f"{gnomad_mod.SOURCE_LABEL} (AF_hom={af:.2%}, n={gnomad_hit['an']:,})."
            )
        if not gnomad_hit["pass_qc"]:
            note += f" gnomAD flags this site: {', '.join(gnomad_hit['filters'])}."
    elif gnomad_checked:
        status = "Rare variant"
        note = (
            f"Homoplasmic call vs rCRS ({het}%); no known disease association. Not found in "
            f"{gnomad_mod.SOURCE_LABEL} (n={gnomad_mod.COHORT_SIZE:,}) — genuinely rare or "
            f"private, not merely uncommon."
        )
    else:
        status = "Common variant"
        note = (
            f"Homoplasmic call vs rCRS ({het}%); no known disease association. Likely "
            f"haplogroup-defining or benign (population frequency not evaluated)."
        )
    return {
        "status": status,
        "statusColor": STATUS_COLORS[status],
        "note": note,
        "disease": None,
        "evidence": "—",
        "inheritance": "Maternal",
        "origin": "novel",
    }
