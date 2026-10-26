"""Parse a MitoRSaw mitochondrial VCF into MitoLR's variant shape.

MitoRSaw (github.com/PacificBiosciences/mitorsaw) calls homoplasmic and
heteroplasmic mtDNA variants from PacBio HiFi reads, filtering NUMT-derived
reads, against the reference (GRCh38 chrM == rCRS / NC_012920.1). Its FORMAT
fields:

    GT   per-haplotype genotype (phased index per haplotype)
    DP   total reads spanning the locus (circular re-spans counted once)
    AD   [ref_reads, alt_reads] definitively assignable (sum <= DP)
    VAF  AD.alt / (AD.ref + AD.alt) -- the alternate-allele (heteroplasmy) fraction

Heteroplasmy is taken directly from VAF; the read pileup uses AD so the
mutant:wild-type ratio reflects real per-allele read support.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Optional

from . import annotate, constraint, gnomad


def _open(path: str):
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, "rt", encoding="utf-8")
    return open(p, "rt", encoding="utf-8")


def _fmt_pct(vaf: float) -> float:
    return round(vaf * 100, 1)


def parse_records(path: str) -> list[dict]:
    """Read raw variant records (chrM/MT only) from a VCF path."""
    records: list[dict] = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10:
                continue
            chrom, pos, _id, ref, alt, qual, filt, _info, fmt, sample = cols[:10]
            if chrom not in ("chrM", "MT", "chrMT", "M"):
                continue
            # First ALT allele only (MitoRSaw emits biallelic records here).
            alt = alt.split(",")[0]
            keys = fmt.split(":")
            vals = sample.split(":")
            f = dict(zip(keys, vals))

            dp = int(f["DP"]) if f.get("DP", ".") not in (".", "") else None
            ad = [int(x) for x in f["AD"].split(",")] if f.get("AD", ".") not in (".", "") else []
            vaf = float(f["VAF"].split(",")[0]) if f.get("VAF", ".") not in (".", "") else None

            records.append(
                {
                    "pos": int(pos),
                    "ref": ref,
                    "alt": alt,
                    "qual": None if qual in (".", "") else float(qual),
                    # '.' FILTER means the caller applied no filter -> the call passed.
                    "filter": "PASS" if filt in (".", "", "PASS") else filt,
                    "gt": f.get("GT", "."),
                    "dp": dp,
                    "ad": ad,
                    "vaf": vaf,
                }
            )
    return records


def _format_gnomad_freq(hit: dict | None) -> str:
    """Compact display string for the Overview table's gnomAD column."""
    if hit is None:
        return "Not in gnomAD"
    af_hom, af_het = hit["af_hom"], hit["af_het"]
    if af_hom >= af_het:
        return f"{af_hom:.2%} (hom)"
    return f"{af_het:.2%} (het)"


def build_variants(path: str, tissue: Optional[str] = None) -> list[dict]:
    """Parse + annotate a MitoRSaw VCF into sorted MitoLR variant dicts.

    ``tissue`` (e.g. "blood", "buccal"), when known for the sample, is passed
    through to annotate.classify() so a negative/low pathogenic-panel finding
    can carry the Parikh et al. (2015) tissue-discordance caveat."""
    out: list[dict] = []
    for r in parse_records(path):
        pos, ref, alt = r["pos"], r["ref"], r["alt"]
        het = _fmt_pct(r["vaf"]) if r["vaf"] is not None else 0.0
        position, change, is_del = annotate.variant_notation(pos, ref, alt)
        locus, consequence, genes = annotate.locus_and_consequence(pos)
        gnomad_hit = gnomad.frequency_for(pos, ref, alt)
        cls = annotate.classify(
            pos, ref, alt, het, gnomad_hit=gnomad_hit, gnomad_checked=True, tissue=tissue
        )
        gene_constraint = constraint.constraint_for_variant(pos, genes)

        # ACMG/AMP evidence tags (McCormick et al. 2020) only for calls that aren't
        # already an exact curated-pathogenic match -- BA1/BS1 benign evidence would
        # contradict an established-pathogenic classification, and PM2/PM5_supporting
        # exist to raise or lower suspicion for variants that aren't already resolved.
        acmg_evidence: list[dict] = []
        if cls["origin"] == "novel":
            acmg_evidence = gnomad.acmg_frequency_evidence(gnomad_hit)
            pm5 = annotate.pm5_supporting_candidate(pos, alt)
            if pm5 is not None:
                acmg_evidence.append(
                    {
                        "code": "PM5_supporting",
                        "strength": "Supporting (pathogenic)",
                        "note": (
                            f"Same position as the curated pathogenic variant {pm5['position']} "
                            f"({pm5['disease']}), different allele."
                        ),
                    }
                )
        # AD gives real per-allele read support for the pileup.
        ad_ref = r["ad"][0] if len(r["ad"]) >= 1 else None
        ad_alt = r["ad"][1] if len(r["ad"]) >= 2 else None

        out.append(
            {
                "id": f"chrM_{pos}_{ref}_{alt}",
                "position": position,
                "change": change,
                "locus": locus,
                "pos": pos,
                "stop": None,
                "type": consequence,
                "consequence": consequence,
                "isDeletion": is_del,
                "qual": r["qual"],
                "filter": r["filter"],
                "gnomadFreq": _format_gnomad_freq(gnomad_hit),
                "gnomad": gnomad_hit,
                "constraint": gene_constraint,
                "acmg_evidence": acmg_evidence,
                "heteroplasmy": het,
                "depth": r["dp"],
                "ad_ref": ad_ref,
                "ad_alt": ad_alt,
                **cls,
            }
        )
    out.sort(key=lambda v: v["pos"])
    return out


def simulate_pileup(variant: dict, display_reads: int = 48) -> list[dict]:
    """IGV-style pileup from real allele depths (AD). Subsamples to display_reads
    while preserving the mutant:wild-type ratio; falls back to heteroplasmy if AD
    is missing. Deterministic layout (no RNG) so it is stable across requests."""
    ad_alt = variant.get("ad_alt")
    ad_ref = variant.get("ad_ref")
    if ad_alt is not None and ad_ref is not None and (ad_alt + ad_ref) > 0:
        frac = ad_alt / (ad_alt + ad_ref)
    else:
        frac = (variant.get("heteroplasmy") or 0) / 100

    n_shown = display_reads
    n_mut = max(0, min(n_shown, round(n_shown * frac)))

    reads = []
    for idx in range(n_shown):
        # Even (Bresenham-style) spread of mutant rows through wild-type rows,
        # so the tick pattern reads as heteroplasmy rather than a solid block.
        is_mut = (idx * n_mut) // n_shown != ((idx + 1) * n_mut) // n_shown
        # Deterministic pseudo-varied span/offset from the row index.
        length = 200 + (idx * 37) % 120
        x = (idx * 13) % 20
        reads.append({"kind": "mut" if is_mut else "wt", "x": float(x), "len": float(length)})
    return reads
