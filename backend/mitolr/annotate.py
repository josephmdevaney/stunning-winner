"""Annotate real mtDNA variant calls against the rCRS feature map and the
curated pathogenic panel.

A VCF alone gives position/ref/alt/depth/heteroplasmy but not functional
consequence (synonymous vs missense needs the codon, hence the reference
sequence) or population frequency. This module annotates what is honestly
derivable — locus, functional class, and any known pathogenic association — and
leaves the rest marked as not evaluated rather than fabricating it.
"""

from __future__ import annotations

import re
from typing import Optional

from . import gnomad as gnomad_mod
from .reference import (
    CONSEQUENCE_BY_TYPE,
    GENOME_LENGTH,
    PATHOGENIC_PANEL,
    STATUS_COLORS,
    features_for_pos,
)

# Known pathogenic SNVs keyed by (pos, ref, alt) for exact-call matching.
_PATHO_BY_KEY: dict[tuple[int, str, str], dict] = {}
for _v in PATHOGENIC_PANEL:
    _m = None
    _pos_str = _v["position"]
    # Parse "m.<pos><ref>><alt>" for SNVs; deletions are indexed separately below.
    _snv = re.match(r"m\.(\d+)([ACGT])>([ACGT])", _pos_str)
    if _snv:
        _PATHO_BY_KEY[(int(_snv.group(1)), _snv.group(2), _snv.group(3))] = _v

# Known pathogenic large-scale deletions keyed by exact (del_start, del_stop)
# breakpoint span -- e.g. the classic KSS/CPEO common deletion, m.8483_13459del4977
# (pos=8483, stop=13459 already ARE that span). Matched by EXACT breakpoint only:
# real samples can carry deletions of similar size at slightly different
# breakpoints, which are related but distinct events, not curated here, so no
# breakpoint-tolerance heuristic is applied without a real basis for one.
_PATHO_DELETIONS: dict[tuple[int, int], dict] = {
    (_v["pos"], _v["stop"]): _v for _v in PATHOGENIC_PANEL if _v.get("isDeletion") and "stop" in _v
}


def deletion_span(pos: int, ref: str, alt: str) -> Optional[tuple[int, int]]:
    """(del_start, del_stop), 1-based, for a MitoRSaw-style anchored deletion
    call (REF longer than ALT) -- the same anchoring convention variant_notation()
    uses for display. None if this call isn't a deletion.

    Wraps del_stop back into 1..GENOME_LENGTH if the deletion spans the circular
    genome's origin (position GENOME_LENGTH+1 is position 1 again) -- mirrors how
    features_for_pos() already treats the D-loop as wrapping, rather than
    silently returning a position beyond the genome's real length."""
    if len(ref) <= len(alt):
        return None
    del_start = pos + len(alt)
    del_stop = pos + len(ref) - 1
    if del_stop > GENOME_LENGTH:
        del_stop -= GENOME_LENGTH
    return (del_start, del_stop)


def _deletion_span_length(del_start: int, del_stop: int) -> int:
    """Length of a deletion span, accounting for del_stop having wrapped past
    the circular genome's origin (so del_stop can be numerically less than
    del_start for a real, non-empty deletion)."""
    if del_stop >= del_start:
        return del_stop - del_start + 1
    return GENOME_LENGTH - del_start + 1 + del_stop


def match_pathogenic(pos: int, ref: str, alt: str) -> Optional[dict]:
    hit = _PATHO_BY_KEY.get((pos, ref, alt))
    if hit is not None:
        return hit
    span = deletion_span(pos, ref, alt)
    return _PATHO_DELETIONS.get(span) if span is not None else None


# Position -> panel entries, restricted to tRNA/rRNA variants, for PM5_supporting.
_PATHO_BY_POS_RNA: dict[int, list[dict]] = {}
for _v in PATHOGENIC_PANEL:
    if _v["type"] in ("tRNA point mutation", "rRNA variant"):
        _PATHO_BY_POS_RNA.setdefault(_v["pos"], []).append(_v)


def pm5_supporting_candidate(pos: int, alt: str) -> Optional[dict]:
    """PM5_supporting (McCormick et al. 2020): a *different* nucleotide change at
    the same position as a curated tRNA/rRNA pathogenic variant raises suspicion
    for pathogenicity, even without a protein-coding amino-acid match -- tRNA/rRNA
    genes aren't translated, so position identity substitutes for the amino-acid
    comparison the original (protein-coding) PM5 rule uses. Scoped to tRNA/rRNA
    panel entries only: the mRNA form of this criterion needs codon-level
    consequence annotation this pipeline doesn't perform.

    Returns the panel entry it matched (for its disease/evidence context), or None.
    Deliberately excludes an exact allele match -- that's already `match_pathogenic`,
    a much stronger line of evidence, not this supporting one."""
    for entry in _PATHO_BY_POS_RNA.get(pos, ()):
        m = re.match(r"m\.(\d+)([ACGT])>([ACGT])", entry["position"])
        if m and m.group(3) != alt:
            return entry
    return None


# Generic heteroplasmy-band reporting convention applied uniformly to every
# curated pathogenic variant. This is NOT a variant-specific phenotype-expression
# threshold -- real mtDNA disease expression varies substantially by variant (e.g.
# m.8993T>G NARP/Leigh typically needs ~70-90%+ heteroplasmy for the severe
# phenotype, while m.3243A>G can produce a MIDD phenotype at much lower levels).
# This build doesn't curate a per-variant expression threshold, so a single
# generic convention is applied to every panel entry and labeled as such below,
# rather than implying variant-specific precision it doesn't have.
PATHOGENIC_HET_THRESHOLD = 30.0
LOW_LEVEL_HET_THRESHOLD = 5.0


def classify_heteroplasmy_band(het: float, disease: str) -> tuple[str, str]:
    """(status, note) from the generic heteroplasmy-band convention above,
    shared by the real-VCF path (classify(), below) and the synthetic demo
    (generate.py) so both apply -- and disclose -- the identical convention."""
    if het >= PATHOGENIC_HET_THRESHOLD:
        return (
            "Pathogenic",
            f"Elevated heteroplasmy ({het}%) above the generic {PATHOGENIC_HET_THRESHOLD:.0f}% "
            f"reporting threshold; consistent with {disease}. This threshold is a uniform "
            f"screening convention applied to every curated variant here, not a "
            f"variant-specific phenotype-expression threshold -- real disease-expression "
            f"heteroplasmy varies substantially by variant.",
        )
    if het >= LOW_LEVEL_HET_THRESHOLD:
        return (
            "Low-level heteroplasmy",
            f"Known pathogenic variant at low-level heteroplasmy ({het}%); below the generic "
            f"{PATHOGENIC_HET_THRESHOLD:.0f}% reporting threshold but may warrant follow-up.",
        )
    return (
        "Not detected",
        f"Known pathogenic position below the generic {LOW_LEVEL_HET_THRESHOLD:.0f}% reporting threshold.",
    )


# Real, citable per-variant heteroplasmy thresholds -- verified directly against
# McCormick et al. 2020 (Hum Mutat 41(12):2028-2057) and the primary papers it
# cites, not a secondary summary. Of the 21 curated panel variants, only these
# have a specific quantitative (or, for LHON, penetrance-based) relationship with
# a real citation; every other panel variant still falls through to the generic
# classify_heteroplasmy_band() convention above because no per-variant threshold
# has been sourced for it -- this is deliberately not an exhaustive system.

def _classify_m8993(het: float, disease: str) -> tuple[str, str]:
    """m.8993T>C/G (NARP/Leigh): disease manifestations typically require blood
    heteroplasmy >90% (Thorburn, Rahman & Rahman, 1993) -- a real, much higher
    threshold than the generic 30% convention."""
    if het >= 90:
        return (
            "Pathogenic",
            f"Heteroplasmy ({het}%) exceeds the ~90% threshold at which m.8993T>C/G "
            f"typically causes disease manifestations (Thorburn, Rahman & Rahman, "
            f"1993); consistent with {disease}.",
        )
    if het >= LOW_LEVEL_HET_THRESHOLD:
        return (
            "Low-level heteroplasmy",
            f"Known pathogenic variant at {het}% heteroplasmy -- below the ~90% "
            f"threshold at which NARP/Leigh manifestations typically occur "
            f"(Thorburn, Rahman & Rahman, 1993), but may still warrant follow-up.",
        )
    return (
        "Not detected",
        f"Known pathogenic position below the generic {LOW_LEVEL_HET_THRESHOLD:.0f}% reporting threshold.",
    )


def _classify_m3243(het: float, disease: str) -> tuple[str, str]:
    """m.3243A>G: real disease-specific heteroplasmy bands (McCormick et al. 2020,
    citing Kaufmann et al. 2009/2011) -- MELAS at 50-70%, Leigh above 90%, MIDD
    below 40%. Laloi-Michelin et al. 2009 report MIDD occurring even below 10%
    heteroplasmy in some carriers, so a low or undetectable reading here does
    NOT rule out MIDD the way it would for the other phenotypes -- that
    exception is disclosed explicitly rather than silently reported as clean."""
    if het >= 90:
        return (
            "Pathogenic",
            f"Heteroplasmy ({het}%) is in the range typically associated with "
            f"Leigh syndrome for m.3243A>G (McCormick et al., 2020).",
        )
    if het >= 50:
        return (
            "Pathogenic",
            f"Heteroplasmy ({het}%) is in the range typically associated with "
            f"MELAS for m.3243A>G (McCormick et al., 2020).",
        )
    if het >= LOW_LEVEL_HET_THRESHOLD:
        return (
            "Low-level heteroplasmy",
            f"Heteroplasmy ({het}%) is below the typical MELAS range (50-70%), but "
            f"m.3243A>G can still cause MIDD (maternally-inherited diabetes and "
            f"deafness) well below that -- reported even below 10% heteroplasmy in "
            f"blood or urine in some carriers (Laloi-Michelin et al., 2009). A low "
            f"reading does not rule out MIDD specifically.",
        )
    return (
        "Not detected",
        f"Below the generic {LOW_LEVEL_HET_THRESHOLD:.0f}% reporting threshold. "
        f"Note: m.3243A>G-associated MIDD has been reported even below 10% "
        f"heteroplasmy in blood or urine (Laloi-Michelin et al., 2009) -- a low or "
        f"undetectable reading here does not fully exclude MIDD.",
    )


def _classify_lhon(het: float, disease: str) -> tuple[str, str]:
    """LHON variants (m.11778G>A, m.3460G>A, m.14484T>C) are frequently
    homoplasmic already -- the clinically operative variable is penetrance
    (~50% of male, ~10% of female carriers ever develop vision loss), not
    heteroplasmy level (Yu-Wai-Man & Chinnery, 1993). The generic band
    convention still runs (a homoplasmic pathogenic call is expected to read
    Pathogenic), but the note discloses that heteroplasmy level itself isn't
    what determines whether disease actually manifests for this variant."""
    status, note = classify_heteroplasmy_band(het, disease)
    if status == "Pathogenic":
        note += (
            " LHON-associated variants are frequently homoplasmic; penetrance "
            "(~50% of male, ~10% of female carriers ever develop vision loss), not "
            "heteroplasmy level, is what actually determines whether disease "
            "manifests (Yu-Wai-Man & Chinnery, 1993)."
        )
    return status, note


# Variant id -> override classifier, for the handful of curated variants above
# with a real per-variant heteroplasmy relationship.
_HETEROPLASMY_OVERRIDES = {
    "m8993TG": _classify_m8993,
    "m3243AG": _classify_m3243,
    "m11778GA": _classify_lhon,
    "m3460GA": _classify_lhon,
    "m14484TC": _classify_lhon,
}


def classify_pathogenic_heteroplasmy(variant_id: str, het: float, disease: str) -> tuple[str, str]:
    """(status, note) for a curated pathogenic-panel variant -- using a real
    per-variant threshold where the literature gives one (_HETEROPLASMY_OVERRIDES
    above), falling back to the generic classify_heteroplasmy_band() convention
    otherwise. Shared by the real-VCF path (classify(), below) and the synthetic
    demo (generate.py) so both apply the identical, real classification rather
    than two copies that could drift apart."""
    override = _HETEROPLASMY_OVERRIDES.get(variant_id)
    if override is not None:
        return override(het, disease)
    return classify_heteroplasmy_band(het, disease)


def tissue_caveat(status: str, variant_id: str, tissue: Optional[str]) -> Optional[str]:
    """A note fragment per the Mitochondrial Medicine Society consensus statement
    (Parikh et al., 2015): a negative or low reading in blood or buccal samples
    does not rule out clinically significant heteroplasmy in another tissue.
    Only relevant when the finding is itself negative/low -- a false negative is
    exactly what the guidance warns about -- and the tested tissue is one of the
    two it specifically flags as insufficient alone. None if there's nothing to
    add (unknown tissue, or a tissue the guidance doesn't address, e.g. urine or
    muscle, or a finding that's already positive)."""
    if tissue not in ("blood", "buccal"):
        return None
    if status not in ("Not detected", "Low-level heteroplasmy"):
        return None
    if variant_id == "m3243AG":
        return (
            f" Sample tissue is {tissue}; per Mitochondrial Medicine Society guidance (Parikh "
            f"et al., 2015), urine heteroplasmy analysis can be more informative and accurate "
            f"than blood or buccal alone for m.3243A>G specifically (Parikh et al., 2015; de "
            f"Laat et al., 2012) -- consider testing another tissue if clinical suspicion "
            f"remains high."
        )
    return (
        f" Sample tissue is {tissue}; per Mitochondrial Medicine Society guidance (Parikh et "
        f"al., 2015), a negative or low reading in blood or buccal samples does not rule out "
        f"clinically significant heteroplasmy in another tissue -- consider testing another "
        f"tissue if clinical suspicion remains high."
    )


def locus_and_consequence(pos: int) -> tuple[str, str, list[str]]:
    """Return (display locus name, consequence class, underlying gene symbols)
    for a genomic position. A position inside one of mtDNA's overlapping
    reading frames (MT-ATP8/MT-ATP6, MT-ND4L/MT-ND4) returns both genes --
    joined for display, and both symbols in the gene list for constraint
    lookup -- rather than only the first-listed gene."""
    feats = features_for_pos(pos)
    if not feats:
        return ("Intergenic", "Noncoding", [])
    names = ["D-loop" if f["name"] == "MT-CR" else f["name"] for f in feats]
    consequence = CONSEQUENCE_BY_TYPE.get(feats[0]["type"], "Noncoding")
    return ("/".join(names), consequence, names)


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
        del_start, del_stop = deletion_span(pos, ref, alt)
        span = _deletion_span_length(del_start, del_stop)
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
    tissue: Optional[str] = None,
) -> dict:
    """Classify a called variant into the display status taxonomy.

    Known-pathogenic calls use the heteroplasmy bands (real per-variant
    thresholds where the literature gives one -- see
    classify_pathogenic_heteroplasmy() -- otherwise the generic convention),
    unaffected by population frequency. When ``tissue`` names the sampled
    tissue (e.g. "blood", "buccal") and the finding is negative or low-level,
    a Parikh et al. (2015) caveat is appended noting that another tissue may
    be needed to rule out clinically significant heteroplasmy -- omitted
    entirely when tissue is unknown, since there's nothing honest to say
    about a tissue that was never recorded.

    Novel (non-pathogenic-panel) calls are reported by heteroplasmy state;
    when ``gnomad_checked`` is True (real VCF-derived calls only -- see
    vcf.py), a homoplasmic novel call is further split into Common vs Rare
    using its real gnomAD v3.1 AF_hom, rather than assuming "homoplasmic
    implies common." ``gnomad_hit`` is None either because the variant was
    checked and not found in the callset (gnomad_checked=True) or because
    gnomAD wasn't consulted at all (gnomad_checked=False, e.g. the synthetic
    demo sample) -- the note text distinguishes the two."""
    patho = match_pathogenic(pos, ref, alt)
    if patho is not None:
        status, note = classify_pathogenic_heteroplasmy(patho["id"], het, patho["disease"])
        caveat = tissue_caveat(status, patho["id"], tissue)
        if caveat is not None:
            note += caveat
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
