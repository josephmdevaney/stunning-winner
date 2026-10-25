"""Infers biological sex from real mosdepth per-chromosome coverage.

A male (XY) sample has one copy of chrX, so its mean depth sits at roughly
half the autosomal mean, and a real chrY shows substantial non-zero coverage.
A female (XX) sample's chrX mean depth matches the autosomal mean (two
copies), and chrY coverage is near zero (only stray cross-mapping noise).

This is a coverage-ratio heuristic, not a clinical karyotype call -- it exists
so chrX/chrY `bcftools roh` output (which genotypes chrX as if diploid
regardless of true copy number) can be correctly read as hemizygosity rather
than autozygosity for a male sample.
"""

from __future__ import annotations

from pathlib import Path

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]

# A real male's chrX ratio clusters near 0.5; a real female's clusters near
# 1.0. A real male's chrY ratio clusters near 0.3-0.5 (lower than chrX's
# due to chrY's poor mappability); a real female's chrY ratio is noise, well
# under 0.1. These thresholds sit in the (wide) gap between those clusters.
CHRX_MALE_MAX_RATIO = 0.75
CHRY_PRESENT_MIN_RATIO = 0.1


def parse_mosdepth_summary(path: str | Path) -> dict[str, dict]:
    """Parse a mosdepth `*.mosdepth.summary.txt` file into {chrom: row}.
    Drops the `<chrom>_region` duplicate rows (identical here since mosdepth
    was run without a separate --by regions file)."""
    rows: dict[str, dict] = {}
    with open(path, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            fields = line.rstrip("\n").split("\t")
            row = dict(zip(header, fields))
            chrom = row["chrom"]
            if chrom.endswith("_region"):
                continue
            rows[chrom] = row
    return rows


def infer_sex_from_mosdepth(path: str | Path) -> dict:
    """Returns sex ('M'/'F'/None) plus the coverage ratios behind the call.
    None means the ratios didn't clearly fall into either cluster -- callers
    should not guess past that, e.g. low-coverage or aneuploid samples."""
    rows = parse_mosdepth_summary(path)
    auto_rows = [rows[c] for c in AUTOSOMES if c in rows]
    if not auto_rows or "chrX" not in rows:
        return {
            "sex": None,
            "method": "mosdepth coverage ratio",
            "reason": "missing autosomal or chrX coverage rows",
        }

    autosomal_len = sum(int(r["length"]) for r in auto_rows)
    autosomal_bases = sum(int(r["bases"]) for r in auto_rows)
    autosomal_mean = autosomal_bases / autosomal_len

    chrx_mean = float(rows["chrX"]["mean"])
    chry_mean = float(rows["chrY"]["mean"]) if "chrY" in rows else 0.0
    chrx_ratio = chrx_mean / autosomal_mean
    chry_ratio = chry_mean / autosomal_mean

    if chrx_ratio <= CHRX_MALE_MAX_RATIO and chry_ratio >= CHRY_PRESENT_MIN_RATIO:
        sex = "M"
    elif chrx_ratio > CHRX_MALE_MAX_RATIO and chry_ratio < CHRY_PRESENT_MIN_RATIO:
        sex = "F"
    else:
        sex = None

    return {
        "sex": sex,
        "method": "mosdepth coverage ratio (chrX/chrY mean depth vs autosomal mean)",
        "autosomal_mean_depth": round(autosomal_mean, 3),
        "chrX_mean_depth": chrx_mean,
        "chrY_mean_depth": chry_mean,
        "chrX_autosomal_ratio": round(chrx_ratio, 3),
        "chrY_autosomal_ratio": round(chry_ratio, 3),
    }
