"""Static mtDNA reference data for MitoLR.

This is seed/reference data the real annotation pipeline would also draw from
(a MITOMAP-style pathogenic-variant knowledge base, an rCRS gene map, and a
haplogroup classifier). In this build it stands in for those services so the
frontend has realistic, stable data to render; it is NOT sequencing output.

Coordinates are 1-based against the rCRS (NC_012920.1, 16,569 bp).
"""

from __future__ import annotations

from typing import Optional, TypedDict


class GeneRegion(TypedDict):
    name: str
    start: int
    stop: int
    type: str  # 'control' | 'rRNA' | 'protein'


class PathogenicVariant(TypedDict, total=False):
    id: str
    position: str
    locus: str
    pos: int
    stop: int
    type: str
    disease: str
    inheritance: str
    evidence: str
    isDeletion: bool


GENOME_LENGTH = 16569

# Ordered annotated regions used for the circular ideogram and for placing
# synthetic/novel calls in plausible loci.
GENE_MAP: list[GeneRegion] = [
    {"name": "D-loop", "start": 16024, "stop": 16569, "type": "control"},
    {"name": "D-loop", "start": 1, "stop": 576, "type": "control"},
    {"name": "MT-RNR1", "start": 648, "stop": 1601, "type": "rRNA"},
    {"name": "MT-RNR2", "start": 1671, "stop": 3229, "type": "rRNA"},
    {"name": "MT-ND1", "start": 3307, "stop": 4262, "type": "protein"},
    {"name": "MT-ND2", "start": 4470, "stop": 5511, "type": "protein"},
    {"name": "MT-CO1", "start": 5904, "stop": 7445, "type": "protein"},
    {"name": "MT-CO2", "start": 7586, "stop": 8269, "type": "protein"},
    {"name": "MT-ATP8/6", "start": 8366, "stop": 9207, "type": "protein"},
    {"name": "MT-CO3", "start": 9207, "stop": 9990, "type": "protein"},
    {"name": "MT-ND3", "start": 10059, "stop": 10404, "type": "protein"},
    {"name": "MT-ND4L/4", "start": 10470, "stop": 12137, "type": "protein"},
    {"name": "MT-ND5", "start": 12337, "stop": 14148, "type": "protein"},
    {"name": "MT-ND6", "start": 14149, "stop": 14673, "type": "protein"},
    {"name": "MT-CYB", "start": 14747, "stop": 15887, "type": "protein"},
]

# Curated well-established pathogenic mtDNA variants (the "known-pathogenic"
# origin panel). Each carries real disease-association metadata.
PATHOGENIC_PANEL: list[PathogenicVariant] = [
    {"id": "m3243AG", "position": "m.3243A>G", "locus": "MT-TL1", "pos": 3243, "type": "tRNA point mutation", "disease": "MELAS / MIDD (maternally inherited diabetes and deafness)", "inheritance": "Maternal", "evidence": "Definitive"},
    {"id": "m3271TC", "position": "m.3271T>C", "locus": "MT-TL1", "pos": 3271, "type": "tRNA point mutation", "disease": "MELAS", "inheritance": "Maternal", "evidence": "Strong"},
    {"id": "m3260AG", "position": "m.3260A>G", "locus": "MT-TL1", "pos": 3260, "type": "tRNA point mutation", "disease": "Maternally inherited cardiomyopathy / myopathy", "inheritance": "Maternal", "evidence": "Moderate"},
    {"id": "m8344AG", "position": "m.8344A>G", "locus": "MT-TK", "pos": 8344, "type": "tRNA point mutation", "disease": "MERRF (myoclonic epilepsy with ragged-red fibers)", "inheritance": "Maternal", "evidence": "Definitive"},
    {"id": "m8356TC", "position": "m.8356T>C", "locus": "MT-TK", "pos": 8356, "type": "tRNA point mutation", "disease": "MERRF/MELAS overlap syndrome", "inheritance": "Maternal", "evidence": "Moderate"},
    {"id": "m8993TG", "position": "m.8993T>G", "locus": "MT-ATP6", "pos": 8993, "type": "Missense", "disease": "NARP / Leigh syndrome", "inheritance": "Maternal", "evidence": "Definitive"},
    {"id": "m9176TC", "position": "m.9176T>C", "locus": "MT-ATP6", "pos": 9176, "type": "Missense", "disease": "Leigh syndrome", "inheritance": "Maternal", "evidence": "Strong"},
    {"id": "m9185TC", "position": "m.9185T>C", "locus": "MT-ATP6", "pos": 9185, "type": "Missense", "disease": "Leigh syndrome / NARP", "inheritance": "Maternal", "evidence": "Moderate"},
    {"id": "m11778GA", "position": "m.11778G>A", "locus": "MT-ND4", "pos": 11778, "type": "Missense", "disease": "Leber hereditary optic neuropathy (LHON)", "inheritance": "Maternal", "evidence": "Definitive"},
    {"id": "m3460GA", "position": "m.3460G>A", "locus": "MT-ND1", "pos": 3460, "type": "Missense", "disease": "Leber hereditary optic neuropathy (LHON)", "inheritance": "Maternal", "evidence": "Definitive"},
    {"id": "m14484TC", "position": "m.14484T>C", "locus": "MT-ND6", "pos": 14484, "type": "Missense", "disease": "Leber hereditary optic neuropathy (LHON)", "inheritance": "Maternal", "evidence": "Definitive"},
    {"id": "m1555AG", "position": "m.1555A>G", "locus": "MT-RNR1", "pos": 1555, "type": "rRNA variant", "disease": "Aminoglycoside-induced & nonsyndromic hearing loss", "inheritance": "Maternal", "evidence": "Definitive"},
    {"id": "m1494CT", "position": "m.1494C>T", "locus": "MT-RNR1", "pos": 1494, "type": "rRNA variant", "disease": "Aminoglycoside-induced hearing loss", "inheritance": "Maternal", "evidence": "Strong"},
    {"id": "m7445AG", "position": "m.7445A>G", "locus": "MT-TS1", "pos": 7445, "type": "tRNA point mutation", "disease": "Sensorineural hearing loss", "inheritance": "Maternal", "evidence": "Moderate"},
    {"id": "m7512TC", "position": "m.7512T>C", "locus": "MT-TS1", "pos": 7512, "type": "tRNA point mutation", "disease": "Myoclonic epilepsy, ataxia, hearing loss", "inheritance": "Maternal", "evidence": "Moderate"},
    {"id": "m4269AG", "position": "m.4269A>G", "locus": "MT-TI", "pos": 4269, "type": "tRNA point mutation", "disease": "Hypertrophic cardiomyopathy", "inheritance": "Maternal", "evidence": "Moderate"},
    {"id": "m5523TC", "position": "m.5523T>C", "locus": "MT-TW", "pos": 5523, "type": "tRNA point mutation", "disease": "Mitochondrial encephalomyopathy", "inheritance": "Maternal", "evidence": "Limited"},
    {"id": "m10197GA", "position": "m.10197G>A", "locus": "MT-ND3", "pos": 10197, "type": "Missense", "disease": "Leigh syndrome / MELAS overlap", "inheritance": "Maternal", "evidence": "Strong"},
    {"id": "m15990CT", "position": "m.15990C>T", "locus": "MT-TP", "pos": 15990, "type": "tRNA point mutation", "disease": "Mitochondrial myopathy", "inheritance": "Maternal", "evidence": "Limited"},
    {"id": "m14709TC", "position": "m.14709T>C", "locus": "MT-TE", "pos": 14709, "type": "tRNA point mutation", "disease": "Myopathy with diabetes mellitus", "inheritance": "Maternal", "evidence": "Moderate"},
    {"id": "kss_common_del", "position": "m.8483_13459del4977", "locus": "—", "pos": 8483, "stop": 13459, "type": "Large-scale deletion", "disease": "Kearns-Sayre syndrome / CPEO", "inheritance": "Sporadic (somatic)", "evidence": "Definitive", "isDeletion": True},
]


# Full rCRS feature map (1-based, NC_012920.1) for precise locus + consequence
# annotation of real variant calls: 2 rRNAs, 22 tRNAs, 13 protein genes, and the
# control region. Coarser than would be needed for codon-level consequence, but
# enough to name the locus and its functional class. tRNA/ND6 strand is ignored
# here (position containment only). The control region wraps the origin.
MT_FEATURES: list[GeneRegion] = [
    {"name": "MT-TF", "start": 577, "stop": 647, "type": "tRNA"},
    {"name": "MT-RNR1", "start": 648, "stop": 1601, "type": "rRNA"},
    {"name": "MT-TV", "start": 1602, "stop": 1670, "type": "tRNA"},
    {"name": "MT-RNR2", "start": 1671, "stop": 3229, "type": "rRNA"},
    {"name": "MT-TL1", "start": 3230, "stop": 3304, "type": "tRNA"},
    {"name": "MT-ND1", "start": 3307, "stop": 4262, "type": "protein"},
    {"name": "MT-TI", "start": 4263, "stop": 4331, "type": "tRNA"},
    {"name": "MT-TQ", "start": 4329, "stop": 4400, "type": "tRNA"},
    {"name": "MT-TM", "start": 4402, "stop": 4469, "type": "tRNA"},
    {"name": "MT-ND2", "start": 4470, "stop": 5511, "type": "protein"},
    {"name": "MT-TW", "start": 5512, "stop": 5579, "type": "tRNA"},
    {"name": "MT-TA", "start": 5587, "stop": 5655, "type": "tRNA"},
    {"name": "MT-TN", "start": 5657, "stop": 5729, "type": "tRNA"},
    {"name": "MT-TC", "start": 5761, "stop": 5826, "type": "tRNA"},
    {"name": "MT-TY", "start": 5826, "stop": 5891, "type": "tRNA"},
    {"name": "MT-CO1", "start": 5904, "stop": 7445, "type": "protein"},
    {"name": "MT-TS1", "start": 7446, "stop": 7514, "type": "tRNA"},
    {"name": "MT-TD", "start": 7518, "stop": 7585, "type": "tRNA"},
    {"name": "MT-CO2", "start": 7586, "stop": 8269, "type": "protein"},
    {"name": "MT-TK", "start": 8295, "stop": 8364, "type": "tRNA"},
    {"name": "MT-ATP8", "start": 8366, "stop": 8572, "type": "protein"},
    {"name": "MT-ATP6", "start": 8527, "stop": 9207, "type": "protein"},
    {"name": "MT-CO3", "start": 9207, "stop": 9990, "type": "protein"},
    {"name": "MT-TG", "start": 9991, "stop": 10058, "type": "tRNA"},
    {"name": "MT-ND3", "start": 10059, "stop": 10404, "type": "protein"},
    {"name": "MT-TR", "start": 10405, "stop": 10469, "type": "tRNA"},
    {"name": "MT-ND4L", "start": 10470, "stop": 10766, "type": "protein"},
    {"name": "MT-ND4", "start": 10760, "stop": 12137, "type": "protein"},
    {"name": "MT-TH", "start": 12138, "stop": 12206, "type": "tRNA"},
    {"name": "MT-TS2", "start": 12207, "stop": 12265, "type": "tRNA"},
    {"name": "MT-TL2", "start": 12266, "stop": 12336, "type": "tRNA"},
    {"name": "MT-ND5", "start": 12337, "stop": 14148, "type": "protein"},
    {"name": "MT-ND6", "start": 14149, "stop": 14673, "type": "protein"},
    {"name": "MT-TE", "start": 14674, "stop": 14742, "type": "tRNA"},
    {"name": "MT-CYB", "start": 14747, "stop": 15887, "type": "protein"},
    {"name": "MT-TT", "start": 15888, "stop": 15953, "type": "tRNA"},
    {"name": "MT-TP", "start": 15956, "stop": 16023, "type": "tRNA"},
]

# Human-readable functional class per feature type, used as the fallback
# "consequence" when codon-level annotation is not available from a VCF alone.
CONSEQUENCE_BY_TYPE = {
    "control": "Control region variant",
    "rRNA": "rRNA variant",
    "tRNA": "tRNA variant",
    "protein": "Protein-coding variant",
}


def feature_for_pos(pos: int) -> Optional[GeneRegion]:
    """Precise locus lookup against MT_FEATURES; control region wraps the origin."""
    if pos >= 16024 or pos <= 576:
        return {"name": "MT-CR", "start": 16024, "stop": 576, "type": "control"}
    for f in MT_FEATURES:
        if f["start"] <= pos <= f["stop"]:
            return f
    return None


class Haplogroup(TypedDict):
    name: str
    confidence: float
    region: str
    markers: list[str]


HAPLOGROUP: Haplogroup = {
    "name": "H2a5b1",
    "confidence": 99.6,
    "region": "Predominantly Western Eurasia",
    "markers": ["A2706G", "T7028C", "A4769G", "A8860G", "A750G", "G4048A"],
}


# Status classification -> badge color (OKLCH), shared with the frontend theme.
STATUS_COLORS = {
    "Pathogenic": "oklch(0.55 0.20 25)",
    "Low-level heteroplasmy": "oklch(0.62 0.15 75)",
    "Not detected": "oklch(0.55 0.14 150)",
    "VUS": "oklch(0.55 0.1 280)",
    "Rare variant": "oklch(0.55 0.12 230)",
    "Common variant": "oklch(0.6 0.02 260)",
}


def region_for_pos(pos: int) -> Optional[GeneRegion]:
    for g in GENE_MAP:
        if g["start"] <= pos <= g["stop"]:
            return g
    return None
