"""Deterministic, seeded generation of a MitoLR sample's variant call set.

Ports the design prototype's seeded PRNG and generation logic to Python so the
data is produced server-side and is stable across requests. In a production
build this module would be replaced by a mitochondrial variant caller's VCF
output joined against a pathogenic-variant knowledge base; the shapes returned
here are what the frontend consumes.

The PRNG (a linear congruential generator) and the sequence of draws mirror the
prototype exactly so the demo tells the same clinical story it was designed to.
"""

from __future__ import annotations

import math
import re
from typing import Callable

from .annotate import classify_pathogenic_heteroplasmy
from .reference import (
    GENE_MAP,
    GENOME_LENGTH,
    HAPLOGROUP,
    PATHOGENIC_PANEL,
    STATUS_COLORS,
)

# Narrative overrides so the sample tells a clear clinical story. Every other
# curated variant gets a plausible low / undetected heteroplasmy.
_OVERRIDES = {
    "m3243AG": 58.0,       # MELAS -- within the real 50-70% MELAS heteroplasmy band
                           # (McCormick et al. 2020); see annotate._classify_m3243
    "m1555AG": 97.0,       # aminoglycoside deafness risk allele -- homoplasmic
    "kss_common_del": 11.0,  # low-level common deletion, below significance
}

_BASES = ("A", "C", "G", "T")
_CHANGE_RE = re.compile(r"(\d+)([ACGT])>([ACGT])")
_ALT_RE = re.compile(r">([ACGT])")


def seeded_rand(seed: int) -> Callable[[], float]:
    """LCG matching the prototype: s = (s*9301 + 49297) % 233280; return s/233280."""
    state = {"s": seed}

    def nxt() -> float:
        state["s"] = (state["s"] * 9301 + 49297) % 233280
        return state["s"] / 233280

    return nxt


def _js_round(x: float) -> int:
    """JavaScript Math.round: round half toward +infinity."""
    return math.floor(x + 0.5)


def _round1(x: float) -> float:
    """Round to 1 decimal, JS-style (Math.round(x*10)/10)."""
    return _js_round(x * 10) / 10


def _classify_pathogenic(het: float, disease: str, variant_id: str) -> tuple[str, str, str]:
    """Thin wrapper around annotate.classify_pathogenic_heteroplasmy() so the demo
    sample and the real-VCF path (annotate.classify()) apply -- and disclose -- the
    identical classification (real per-variant threshold where the literature
    gives one, generic band convention otherwise) rather than two copies that
    could silently drift apart."""
    status, note = classify_pathogenic_heteroplasmy(variant_id, het, disease)
    return status, STATUS_COLORS[status], note


def build_pathogenic_panel() -> list[dict]:
    out: list[dict] = []
    for v in PATHOGENIC_PANEL:
        rng = seeded_rand(7 + sum(ord(c) for c in v["id"]))
        het = _round1(rng() * 3)  # draw regardless, then override
        if v["id"] in _OVERRIDES:
            het = _OVERRIDES[v["id"]]

        status, status_color, note = _classify_pathogenic(het, v.get("disease", ""), v["id"])
        m = _CHANGE_RE.search(v["position"])
        if m:
            change = f"{m.group(2)}>{m.group(3)}"
        else:
            change = "del" if v.get("isDeletion") else "—"

        qual = 60 + math.floor(rng() * 39)
        depth = 3800 + math.floor(rng() * 1800)

        out.append(
            {
                **v,
                "change": change,
                "consequence": v["type"],
                "qual": qual,
                "filter": "PASS",
                "gnomadFreq": "<0.0001",
                "origin": "pathogenic",
                "heteroplasmy": het,
                "status": status,
                "statusColor": status_color,
                "note": note,
                "depth": depth,
            }
        )
    return out


def build_novel_variants(count: int) -> list[dict]:
    """Simulate the broader benign/VUS calls a real mtDNA VCF would carry vs rCRS."""
    patho_positions = {v["pos"] for v in PATHOGENIC_PANEL}
    rng = seeded_rand(1337)
    cons_by_type = {
        "control": ["Control region variant", "Noncoding"],
        "rRNA": ["rRNA variant"],
        "protein": ["Synonymous", "Synonymous", "Missense VUS"],
    }

    out: list[dict] = []
    tries = 0
    while len(out) < count and tries < count * 8:
        tries += 1
        region = GENE_MAP[math.floor(rng() * len(GENE_MAP))]
        span = max(1, region["stop"] - region["start"])
        pos = region["start"] + math.floor(rng() * span)
        if pos in patho_positions or any(o["pos"] == pos for o in out):
            continue
        ref = _BASES[math.floor(rng() * 4)]
        alt = _BASES[math.floor(rng() * 4)]
        while alt == ref:
            alt = _BASES[math.floor(rng() * 4)]

        cons_options = cons_by_type.get(region["type"], ["Noncoding"])
        consequence = cons_options[math.floor(rng() * len(cons_options))]
        qual = 8 + math.floor(rng() * 11) if rng() < 0.08 else 35 + math.floor(rng() * 64)
        filter_ = "PASS" if qual >= 20 else "low_qual"
        rare = rng() < 0.12
        gnomad_freq = f"{0.0001 + rng() * 0.009:.4f}" if rare else f"{0.05 + rng() * 0.94:.2f}"
        is_het = rng() < 0.07
        heteroplasmy = _round1(10 + rng() * 55) if is_het else _round1(95 + rng() * 5)

        if is_het:
            status = "VUS"
            note = (
                f"Heteroplasmic call ({heteroplasmy}%) with no established disease association "
                f"— variant of uncertain significance."
            )
        elif rare:
            status = "Rare variant"
            note = (
                f"Homoplasmic call, low population frequency ({gnomad_freq}); "
                f"no known disease association."
            )
        else:
            status = "Common variant"
            note = (
                f"Homoplasmic, common population polymorphism (freq {gnomad_freq}); "
                f"likely haplogroup-defining or benign."
            )

        out.append(
            {
                "id": f"novel_{pos}_{ref}{alt}",
                "position": f"m.{pos}{ref}>{alt}",
                "change": f"{ref}>{alt}",
                "locus": region["name"],
                "pos": pos,
                "stop": None,
                "type": consequence,
                "consequence": consequence,
                "disease": None,
                "inheritance": "Maternal",
                "evidence": "—",
                "isDeletion": False,
                "qual": qual,
                "filter": filter_,
                "gnomadFreq": gnomad_freq,
                "origin": "novel",
                "heteroplasmy": heteroplasmy,
                "status": status,
                "statusColor": STATUS_COLORS[status],
                "note": note,
                "depth": 3600 + math.floor(rng() * 2200),
            }
        )
    return out


def build_variants() -> list[dict]:
    variants = build_pathogenic_panel() + build_novel_variants(58)
    variants.sort(key=lambda v: v["pos"])
    return variants


# A real locus is covered by thousands of reads; an IGV-style pileup shows a
# legible subsample. Draw this many rows, preserving the heteroplasmy ratio.
PILEUP_DISPLAY_READS = 48


def simulate_pileup(variant: dict) -> list[dict]:
    """IGV-style per-read layout: one row per (subsampled) read with a
    deterministic span and offset, tagged mutant vs wild-type in proportion to
    heteroplasmy. Subsampled to PILEUP_DISPLAY_READS so the track stays legible
    even for high-depth (~thousands of reads) calls."""
    het = variant.get("heteroplasmy") or 0
    n_shown = min(PILEUP_DISPLAY_READS, variant.get("depth") or PILEUP_DISPLAY_READS)
    n_mut = max(0, _js_round(n_shown * het / 100))
    rows = ["mut"] * n_mut + ["wt"] * (n_shown - n_mut)

    seed = 3 + sum(ord(c) for c in (variant.get("id") or "x"))
    rng = seeded_rand(seed)
    # Deterministic shuffle (Fisher-Yates), matching the prototype.
    for i in range(len(rows) - 1, 0, -1):
        j = math.floor(rng() * (i + 1))
        rows[i], rows[j] = rows[j], rows[i]

    reads = []
    for kind in rows:
        length = 200 + rng() * 120
        x = rng() * 20
        reads.append({"kind": kind, "x": round(x, 2), "len": round(length, 2)})
    return reads


def alt_base(variant: dict) -> str:
    m = _ALT_RE.search(variant.get("position", ""))
    return m.group(1) if m else "N"


def haplogroup() -> dict:
    return dict(HAPLOGROUP)


def genome_length() -> int:
    return GENOME_LENGTH
