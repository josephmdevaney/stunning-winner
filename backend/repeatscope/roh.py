"""Runs-of-Homozygosity (ROH) analysis.

Parses real `bcftools roh` RG (region-summary) output and derives, from those
real regions only:

* position classification (telomeric / pericentromeric / interstitial) against
  real GRCh38 chromosome lengths and centromere midpoints
* gene overlap against a curated OMIM disease-gene model (autozygosity
  mapping / AR-candidate-gene flagging)
* an autosomal inbreeding-coefficient (F) estimate from total ROH length
* a UPD (uniparental disomy) screening flag using a commonly cited but
  NOT independently verified SNP-microarray ROH heuristic -- see
  roh_reference.py's module docstring

chrX/chrY/chrM regions are reported (bcftools genuinely calls "ROH" there)
but excluded from the autosomal F-value: chrX/Y homozygosity in a male
sample is hemizygosity, not autozygosity, and chrM has no meaningful
diploid homozygosity concept at all.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Optional

from . import roh_reference as ref
from . import sex_inference

BACKEND_ROOT = Path(__file__).resolve().parents[1]

# sample -> real bcftools roh RG-output. HG002 is the public GIAB fixture
# (committed); CDS27601 is real per-participant hackathon data, kept outside
# the repo (never committed) -- see backend/config/modes.local.json's own
# genome_wide entry for the matching TRGT data for this same sample.
ROH_SAMPLES: dict[str, str] = {
    "HG002": str(BACKEND_ROOT / "fixtures" / "roh" / "HG002.GRCh38.roh.txt"),
    "CDS27601": "/Users/joedevaney/hackathon_cds2760x/CDS27601.GRCh38.roh.txt",
    "CDS48001": "/Users/joedevaney/hackathon_cds48001/CDS48001.GRCh38.roh.txt",
    "CDS57401": "/Users/joedevaney/hackathon_cds57401/CDS57401.GRCh38.roh.txt",
}

# sample -> real mosdepth `*.mosdepth.summary.txt`, used to infer biological
# sex from chrX/chrY coverage relative to the autosomal mean (see
# sex_inference.py). This is what lets chrX/chrY ROH calls below be read as
# hemizygosity rather than autozygosity for a male sample.
MOSDEPTH_SUMMARIES: dict[str, str] = {
    "CDS27601": "/Users/joedevaney/hackathon_cds2760x/CDS27601.GRCh38.mosdepth.summary.txt",
}

# Fallback for samples with no mosdepth coverage in this build. HG002/NA24385
# is the GIAB trio's father -- documented public pedigree fact, not something
# this pipeline infers from its own data.
KNOWN_SEX: dict[str, str] = {
    "HG002": "M",
}


@lru_cache(maxsize=None)
def sex_for_sample(sample: str) -> dict:
    """Best-available sex call for `sample`: real coverage-based inference
    when a mosdepth summary is on file, else known external provenance, else
    unknown (never guessed)."""
    mosdepth_path = MOSDEPTH_SUMMARIES.get(sample)
    if mosdepth_path and Path(mosdepth_path).exists():
        result = sex_inference.infer_sex_from_mosdepth(mosdepth_path)
        if result["sex"]:
            return result
    known = KNOWN_SEX.get(sample)
    if known:
        return {
            "sex": known,
            "method": "known sample provenance (not inferred from coverage in this build)",
        }
    return {"sex": None, "method": None}


@dataclass
class RohRegion:
    sample: str
    chrom: str
    start: int
    end: int
    length_bp: int
    n_markers: int
    quality: float


def parse_roh_file(path: str | Path) -> list[RohRegion]:
    """Parse RG lines from real `bcftools roh` output. Ignores ST (per-marker)
    lines and comment/header lines -- only RG region summaries are needed."""
    regions: list[RohRegion] = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.startswith("RG\t"):
                continue
            # A stray empty tab-separated field has been observed in this fixture
            # (paste artifact, not a real bcftools output shape) -- drop empty
            # fields rather than fail, since a genuine RG row never has one.
            fields = [x for x in line.rstrip("\n").split("\t") if x != ""]
            _, sample, chrom, start, end, length_bp, n_markers, quality = fields
            regions.append(
                RohRegion(
                    sample=sample,
                    chrom=chrom,
                    start=int(start),
                    end=int(end),
                    length_bp=int(length_bp),
                    n_markers=int(n_markers),
                    quality=float(quality),
                )
            )
    return regions


def classify_position(chrom: str, start: int, end: int) -> Optional[str]:
    """Telomeric / pericentromeric / interstitial, using real GRCh38 chromosome
    length + centromere midpoint reference data. None for chromosomes outside
    that reference (chrY, chrM)."""
    length = ref.CHROM_LENGTH.get(chrom)
    if length is None:
        return None
    near_start = start <= ref.TELOMERE_PROXIMITY_BP
    near_end = end >= length - ref.TELOMERE_PROXIMITY_BP
    if near_start or near_end:
        return "telomeric"
    cen = ref.CENTROMERE_MID.get(chrom)
    if cen is not None and start <= cen <= end:
        return "pericentromeric"
    return "interstitial"


def length_band(length_bp: int) -> str:
    mb = length_bp / 1_000_000
    if mb >= 5:
        return "≥5 Mb"
    if mb >= 3:
        return "3–5 Mb"
    if mb >= 1:
        return "1–3 Mb"
    return "<1 Mb"


def gene_overlaps(chrom: str, start: int, end: int) -> list[ref.GeneModelEntry]:
    return [
        g for g in ref.GENE_MODEL
        if g["chrom"] == chrom and g["start"] <= end and g["end"] >= start
    ]


def upd_flag(chrom: str, position: Optional[str], length_bp: int) -> Optional[str]:
    """Screening-only UPD flag. See roh_reference.py docstring: thresholds
    follow a commonly cited SNP-microarray convention not independently
    verified against a primary source in this codebase."""
    if chrom not in ref.AUTOSOMES:
        return None
    imprinted = chrom in ref.IMPRINTED_CHROMS
    if position == "telomeric" and length_bp >= ref.UPD_TELOMERIC_IMPRINTED_BP:
        return (
            f"Possible UPD -- telomeric ≥{ref.UPD_TELOMERIC_IMPRINTED_BP // 1_000_000} Mb"
            f"{' (imprinted chromosome)' if imprinted else ''}"
        )
    if position == "interstitial" and imprinted and length_bp >= ref.UPD_INTERSTITIAL_IMPRINTED_BP:
        return f"Possible UPD -- interstitial ≥{ref.UPD_INTERSTITIAL_IMPRINTED_BP // 1_000_000} Mb (imprinted chromosome)"
    if position == "interstitial" and not imprinted and length_bp >= ref.UPD_INTERSTITIAL_NONIMPRINTED_BP:
        return f"Possible UPD -- interstitial ≥{ref.UPD_INTERSTITIAL_NONIMPRINTED_BP // 1_000_000} Mb"
    return None


def _non_autosomal_note(sex_info: dict) -> str:
    sex = sex_info["sex"]
    if sex == "M":
        if "chrX_mean_depth" in sex_info:
            evidence = (
                f"Inferred sex: male (XY), from real mosdepth coverage -- chrX mean depth "
                f"{sex_info['chrX_mean_depth']:.2f}x is {sex_info['chrX_autosomal_ratio']:.0%} of the "
                f"autosomal mean ({sex_info['autosomal_mean_depth']:.2f}x, consistent with one copy), "
                f"and chrY mean depth {sex_info['chrY_mean_depth']:.2f}x confirms a Y chromosome is "
                f"present."
            )
        else:
            evidence = "Sample sex: male (XY), from known sample provenance (no mosdepth coverage on file)."
        return (
            f"{evidence} chrX/chrY regions below are real bcftools roh calls, but bcftools genotypes "
            "chrX/chrY as if diploid -- for a male sample this makes the entire non-PAR region trivially "
            "\"homozygous,\" reflecting hemizygosity (one copy), not autozygosity. They are excluded from "
            "the autosomal F-value and AR-candidate-gene search below. Mitochondrial DNA (chrM) has no "
            "diploid homozygosity concept either."
        )
    if sex == "F":
        evidence = "Inferred sex: female (XX), from real mosdepth coverage."
        return (
            f"{evidence} Unlike a male sample, chrX here is genuinely diploid, so chrX ROH regions may "
            "reflect real autozygosity -- however they are still conservatively excluded from the "
            "autosomal F-value and AR-candidate-gene search below (this build's reference gene model and "
            "F-value baseline are autosome-only; extending them to chrX would need dedicated X-chromosome "
            "reference data). chrY coverage is negligible as expected. Mitochondrial DNA (chrM) has no "
            "diploid homozygosity concept."
        )
    return (
        "chrX, chrY, and chrM regions are real bcftools roh calls but are excluded from the autosomal "
        "F-value: this sample's sex could not be determined (no mosdepth coverage on file and no known "
        "provenance), and chrX/chrY homozygosity means something different depending on sex -- "
        "hemizygosity (not autozygosity) for a male, potentially real autozygosity for a female. "
        "Mitochondrial DNA has no diploid homozygosity concept regardless."
    )


def _region_dict(r: RohRegion) -> dict:
    position = classify_position(r.chrom, r.start, r.end)
    genes = gene_overlaps(r.chrom, r.start, r.end)
    sex = sex_for_sample(r.sample)["sex"]
    hemizygous = sex == "M" and r.chrom in ("chrX", "chrY")
    return {
        "sample": r.sample,
        "chrom": r.chrom,
        "start": r.start,
        "end": r.end,
        "length_bp": r.length_bp,
        "length_mb": round(r.length_bp / 1_000_000, 2),
        "n_markers": r.n_markers,
        "quality": r.quality,
        "position": position,
        "length_band": length_band(r.length_bp),
        "imprinted_chrom": r.chrom in ref.IMPRINTED_CHROMS,
        "genes": [
            {
                "symbol": g["sym"], "inheritance": g["inh"], "disease": g["disease"], "omim": g["omim"],
                "start": g["start"], "end": g["end"],
            }
            for g in genes
        ],
        "upd_flag": upd_flag(r.chrom, position, r.length_bp),
        "autosomal": r.chrom in ref.AUTOSOMES,
        "f_value_eligible": r.chrom in ref.AUTOSOMES and r.length_bp >= ref.F_VALUE_MIN_LENGTH_BP,
        "hemizygous": hemizygous,
    }


def build_regions(sample: str, path: str | Path) -> list[dict]:
    return [_region_dict(r) for r in parse_roh_file(path)]


@lru_cache(maxsize=None)
def regions_for_sample(sample: str) -> list[dict]:
    path = ROH_SAMPLES[sample]
    return build_regions(sample, path)


def build_overview(sample: str, path: str | Path) -> dict:
    regions = build_regions(sample, path)
    sex_info = sex_for_sample(sample)
    sex = sex_info["sex"]

    autosomal = [r for r in regions if r["autosomal"]]
    f_eligible = [r for r in autosomal if r["f_value_eligible"]]

    total_f_length_bp = sum(r["length_bp"] for r in f_eligible)
    f_value = total_f_length_bp / ref.AUTOSOMAL_GENOME_LENGTH

    upd_flags = [r for r in regions if r["upd_flag"]]

    # Surface AR genes overlapped by an ROH >= AR_CANDIDATE_MIN_LENGTH_BP --
    # deliberately a lower, separate floor than the F-value's own >=3 Mb
    # (see AR_CANDIDATE_MIN_LENGTH_BP's docstring): short background ROH are
    # ubiquitous genome-wide and overlapping one by chance isn't a meaningful
    # autozygosity-mapping candidate, but requiring the F-value's stricter
    # floor here would silently drop real candidates in a smaller genuinely
    # -autozygous ROH. Keep the largest overlapping region per gene if more
    # than one qualifies.
    ar_eligible = [r for r in autosomal if r["length_bp"] >= ref.AR_CANDIDATE_MIN_LENGTH_BP]
    ar_genes: dict[str, dict] = {}
    for r in ar_eligible:
        for g in r["genes"]:
            if g["inheritance"] != "AR":
                continue
            key = g["symbol"]
            if key not in ar_genes or r["length_mb"] > ar_genes[key]["roh_length_mb"]:
                ar_genes[key] = {
                    "symbol": g["symbol"], "disease": g["disease"], "omim": g["omim"],
                    "chrom": r["chrom"], "roh_length_mb": r["length_mb"],
                }

    by_chrom: dict[str, dict] = {}
    for chrom in ref.CHROM_ORDER:
        chrom_regions = [r for r in regions if r["chrom"] == chrom]
        by_chrom[chrom] = {
            "chrom": chrom,
            "chrom_length_bp": ref.CHROM_LENGTH.get(chrom),
            "centromere_mid_bp": ref.CENTROMERE_MID.get(chrom),
            "n_regions": len(chrom_regions),
            "total_length_bp": sum(r["length_bp"] for r in chrom_regions),
            "imprinted": chrom in ref.IMPRINTED_CHROMS,
            "hemizygous": sex == "M" and chrom in ("chrX", "chrY"),
        }

    return {
        "sample": sample,
        "genome_build": "GRCh38",
        "caller": "bcftools roh",
        "by_chromosome": [by_chrom[c] for c in ref.CHROM_ORDER if by_chrom[c]["n_regions"] or c in ref.AUTOSOMES],
        "stats": {
            "n_regions_total": len(regions),
            "n_regions_autosomal": len(autosomal),
            "n_regions_chrX": len([r for r in regions if r["chrom"] == "chrX"]),
            "n_regions_chrY": len([r for r in regions if r["chrom"] == "chrY"]),
            "n_regions_chrM": len([r for r in regions if r["chrom"] == "chrM"]),
            "n_regions_f_eligible": len(f_eligible),
            "total_autosomal_roh_bp": sum(r["length_bp"] for r in autosomal),
            "total_f_eligible_roh_bp": total_f_length_bp,
            "f_value_min_length_mb": ref.F_VALUE_MIN_LENGTH_BP // 1_000_000,
            "ar_candidate_min_length_mb": ref.AR_CANDIDATE_MIN_LENGTH_BP // 1_000_000,
            "autosomal_genome_length_bp": ref.AUTOSOMAL_GENOME_LENGTH,
            "f_value": round(f_value, 5),
            "f_value_pct": round(f_value * 100, 3),
            "n_upd_flags": len(upd_flags),
            "n_ar_candidate_genes": len(ar_genes),
        },
        "upd_flags": upd_flags,
        "ar_candidate_genes": sorted(ar_genes.values(), key=lambda g: -g["roh_length_mb"]),
        "relationship_f_values": [
            {"relationship": name, "f_value": val, "f_value_pct": round(val * 100, 4)}
            for name, val in ref.RELATIONSHIP_F_VALUES
        ],
        "sex_inference": sex_info,
        "non_autosomal_note": _non_autosomal_note(sex_info),
        "upd_methodology_note": (
            "UPD flags are a screening heuristic (telomeric ≥5 Mb / interstitial ≥10 Mb on an "
            "imprinted chromosome / interstitial ≥15 Mb elsewhere), commonly attributed to "
            "SNP-microarray ROH/UPD screening conventions (Hoppman-Chaney-style criteria). This "
            "attribution has not been independently verified against a primary source in this build -- "
            "treat flags as candidates for orthogonal confirmation (methylation-specific analysis, "
            "microsatellite/trio testing), not as a diagnostic call."
        ),
    }


def sample_keys() -> list[str]:
    return list(ROH_SAMPLES)


@lru_cache(maxsize=None)
def overview_for_sample(sample: str) -> dict:
    return build_overview(sample, ROH_SAMPLES[sample])


def query_regions(
    sample: str,
    chrom: str = "all",
    min_len_mb: float = 0,
    page: int = 0,
    page_size: int = 100,
) -> dict:
    all_regions = regions_for_sample(sample)
    filtered = all_regions
    if chrom != "all":
        filtered = [r for r in filtered if r["chrom"] == chrom]
    if min_len_mb:
        filtered = [r for r in filtered if r["length_mb"] >= min_len_mb]

    total = len(filtered)
    start = page * page_size
    page_rows = filtered[start : start + page_size]
    chroms = [c for c in ref.CHROM_ORDER if any(r["chrom"] == c for r in all_regions)]

    return {
        "total_all": len(all_regions),
        "total_matching": total,
        "chromosomes": ["all"] + chroms,
        "page": page,
        "page_size": page_size,
        "rows": page_rows,
    }
