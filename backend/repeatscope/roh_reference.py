"""Static reference data for the Runs-of-Homozygosity (ROH) tool.

Chromosome lengths and centromere midpoints are standard GRCh38 assembly
values (matching UCSC/NCBI GRCh38 chromosome sizes and cytoBand centromere
ranges); they are public reference facts, not sequencing output.

The OMIM disease-gene model below is user-curated (GRCh38 coordinates
verified against NCBI Gene / OMIM by the user), used for gene-overlap /
autozygosity-mapping against real ROH regions.

The UPD-screening length thresholds follow a commonly cited SNP-microarray
ROH/UPD screening convention (telomeric >=5 Mb, interstitial >=10 Mb on an
imprinted chromosome, interstitial >=15 Mb elsewhere -- often attributed to
Hoppman-Chaney et al.'s work on unmasking autosomal recessive disease via
SNP array ROH). That attribution has NOT been independently verified against
a primary source in this codebase, so it is surfaced as a labeled screening
heuristic, not a cited/authoritative threshold.
"""

from __future__ import annotations

from typing import Optional, TypedDict

# GRCh38 chromosome lengths, bp. Standard assembly values (chr1-22, chrX).
# chrY and chrM are handled separately -- see CHROM_ORDER note below.
CHROM_LENGTH: dict[str, int] = {
    "chr1": 248_956_422, "chr2": 242_193_529, "chr3": 198_295_559,
    "chr4": 190_214_555, "chr5": 181_538_259, "chr6": 170_805_979,
    "chr7": 159_345_973, "chr8": 145_138_636, "chr9": 138_394_717,
    "chr10": 133_797_422, "chr11": 135_086_622, "chr12": 133_275_309,
    "chr13": 114_364_328, "chr14": 107_043_718, "chr15": 101_991_189,
    "chr16": 90_338_345, "chr17": 83_257_441, "chr18": 80_373_285,
    "chr19": 58_617_616, "chr20": 64_444_167, "chr21": 46_709_983,
    "chr22": 50_818_468, "chrX": 156_040_895, "chrY": 57_227_415,
}

# Approximate GRCh38 centromere midpoints, bp (from UCSC cytoBand centromere
# ranges). Used only to classify a region as spanning the centromere
# (pericentromeric), not as an exact acrocentric/metacentric boundary.
CENTROMERE_MID: dict[str, int] = {
    "chr1": 123_400_000, "chr2": 93_900_000, "chr3": 90_900_000,
    "chr4": 50_000_000, "chr5": 48_800_000, "chr6": 59_800_000,
    "chr7": 60_100_000, "chr8": 45_200_000, "chr9": 43_000_000,
    "chr10": 39_800_000, "chr11": 53_400_000, "chr12": 35_500_000,
    "chr13": 17_700_000, "chr14": 17_200_000, "chr15": 19_000_000,
    "chr16": 36_800_000, "chr17": 25_100_000, "chr18": 18_500_000,
    "chr19": 26_200_000, "chr20": 28_100_000, "chr21": 12_000_000,
    "chr22": 15_000_000, "chrX": 61_000_000,
}

# Autosomes in karyotype order, then chrX. chrY/chrM are reported (real
# bcftools calls exist for them) but excluded from autosomal F-value/
# consanguinity math -- see roh.py.
AUTOSOMES: list[str] = [f"chr{i}" for i in range(1, 23)]
CHROM_ORDER: list[str] = AUTOSOMES + ["chrX", "chrY", "chrM"]

# Imprinted autosomes (chromosomes carrying at least one well-established
# imprinted locus), used by the UPD screening heuristic below.
IMPRINTED_CHROMS: set[str] = {"chr6", "chr7", "chr11", "chr14", "chr15", "chr20"}

# Total autosomal genome length, bp -- denominator for the ROH-based
# inbreeding-coefficient (F) estimate.
AUTOSOMAL_GENOME_LENGTH = sum(CHROM_LENGTH[c] for c in AUTOSOMES)

# Region near a chromosome terminus counted as "telomeric" for position
# classification, bp.
TELOMERE_PROXIMITY_BP = 1_000_000

# ROH length thresholds for the UPD screening heuristic, bp. See module
# docstring: commonly cited SNP-microarray convention, not independently
# verified against a primary source in this codebase.
UPD_TELOMERIC_IMPRINTED_BP = 5_000_000
UPD_INTERSTITIAL_IMPRINTED_BP = 10_000_000
UPD_INTERSTITIAL_NONIMPRINTED_BP = 15_000_000

# Minimum ROH length counted toward the autosomal F-value estimate, bp.
# Short ROH are common in outbred individuals from background LD/drift, not
# just IBD; ROH analyses commonly restrict "IBD-informative" length to a
# multi-Mb floor. 3 Mb is used here as a conventional floor for that
# distinction (not tied to a specific cited threshold).
F_VALUE_MIN_LENGTH_BP = 3_000_000

# Theoretical inbreeding coefficient (F) by pedigree relationship degree,
# derived from the standard F = (1/2)^n formula for n meioses separating
# the parents' common ancestor. Textbook population-genetics values, not
# from a disease-specific guideline.
RELATIONSHIP_F_VALUES: list[tuple[str, float]] = [
    ("Parent-child / full siblings", 0.25),
    ("Uncle-niece / half-siblings / double first cousins", 0.125),
    ("First cousins", 0.0625),
    ("First cousins once removed", 0.03125),
    ("Second cousins", 0.015625),
    ("Second cousins once removed", 0.0078125),
    ("Third cousins", 0.00390625),
]


class GeneModelEntry(TypedDict):
    sym: str
    chrom: str
    start: int
    end: int
    inh: str  # 'AR' | 'AD'
    disease: str
    omim: str


# Curated OMIM disease-gene model, GRCh38 coordinates. User-provided;
# coordinates verified against NCBI Gene / OMIM.
GENE_MODEL: list[GeneModelEntry] = [
    {"sym": "ZAP70", "chrom": "chr2", "start": 97_713_576, "end": 97_756_364,
     "inh": "AR", "disease": "Immunodeficiency 48 (SCID)", "omim": "176947"},
    {"sym": "CNGA3", "chrom": "chr2", "start": 98_346_188, "end": 98_399_797,
     "inh": "AR", "disease": "Achromatopsia 2", "omim": "600053"},
    {"sym": "HFE", "chrom": "chr6", "start": 26_087_429, "end": 26_098_343,
     "inh": "AR", "disease": "Hemochromatosis, type 1", "omim": "613609"},
    {"sym": "SHH", "chrom": "chr7", "start": 155_799_980, "end": 155_812_463,
     "inh": "AD", "disease": "Holoprosencephaly 3", "omim": "600725"},
    {"sym": "MNX1", "chrom": "chr7", "start": 156_994_915, "end": 157_004_509,
     "inh": "AD", "disease": "Currarino syndrome", "omim": "142994"},
    {"sym": "TJP2", "chrom": "chr9", "start": 69_121_264, "end": 69_255_208,
     "inh": "AR", "disease": "Cholestasis, progressive familial intrahepatic 4", "omim": "607709"},
    {"sym": "TMC1", "chrom": "chr9", "start": 72_521_608, "end": 72_838_297,
     "inh": "AR", "disease": "Deafness, DFNB7/11 (AR); DFNA36 (AD)", "omim": "606706"},
    {"sym": "RB1", "chrom": "chr13", "start": 48_303_751, "end": 48_481_890,
     "inh": "AD", "disease": "Retinoblastoma", "omim": "614041"},
    {"sym": "OTOA", "chrom": "chr16", "start": 21_689_000, "end": 21_819_000,
     "inh": "AR", "disease": "Deafness, AR 22", "omim": "607038"},
    {"sym": "COG7", "chrom": "chr16", "start": 23_368_000, "end": 23_432_000,
     "inh": "AR", "disease": "Congenital disorder of glycosylation IIe", "omim": "606978"},
    {"sym": "CLN3", "chrom": "chr16", "start": 28_466_000, "end": 28_490_000,
     "inh": "AR", "disease": "Neuronal ceroid lipofuscinosis 3 (Batten)", "omim": "607042"},
    {"sym": "TUFM", "chrom": "chr16", "start": 28_850_000, "end": 28_857_000,
     "inh": "AR", "disease": "Combined oxidative phosphorylation deficiency 4", "omim": "602389"},
]
