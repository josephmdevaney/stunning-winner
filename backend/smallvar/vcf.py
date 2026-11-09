"""Streaming reader for a single-sample small-variant VCF.

Reads what PacBio's HiFi WGS pipeline produces (DeepVariant, phased by HiPhase:
``GT:GQ:DP:AD:VAF:PL:PS``) and anything shaped like it, such as a GIAB truth set (``GT``
only). A genome has 4-5 million calls, so records are yielded one at a time and nothing
is kept: callers aggregate as they go.

One record can carry several ALT alleles. Each ALT that the genotype actually carries is
yielded as its own ``Call`` with the number of copies the sample has of it, which is how a
``1/2`` site becomes two heterozygous calls in trans.
"""

from __future__ import annotations

import gzip
import re
from typing import Dict, Iterator, List, NamedTuple, Optional

from .normalize import canon_chrom, trim_alleles, variant_kind

PASS_FILTERS = ("PASS", ".")
_GT_SPLIT = re.compile(r"[/|]")
_SYMBOLIC_ALT = re.compile(r"^(<.*>|\*|.*[\[\]].*|\.)$")


class Call(NamedTuple):
    chrom: str  # canonical, chrN
    pos: int  # as written in the VCF
    ref: str
    alt: str
    kind: str  # snv | mnv | ins | del | complex, from the trimmed alleles
    copies: int  # copies of this ALT in the genotype: 1 het, 2 hom
    gt: str  # genotype as written
    phased: bool
    ps: Optional[str]  # phase set
    haplotypes: str  # haplotypes carrying this ALT ("1", "2" or "12"); "" if unphased
    qual: Optional[float]
    filter: str
    dp: Optional[int]
    gq: Optional[int]
    vaf: Optional[float]
    norm_pos: int  # trimmed position/alleles: the form every database join uses
    norm_ref: str
    norm_alt: str


class ReadStats:
    """What the reader saw besides the calls it yielded."""

    def __init__(self) -> None:
        self.records = 0
        self.no_alt_call = 0  # hom-ref / no-call records
        self.symbolic_or_star = 0  # ALT alleles that are not plain sequence
        self.samples = 0  # sample columns in the file (only the first is read)
        self.sample_name: Optional[str] = None


def open_text(path: str):
    p = str(path)
    if p.endswith((".gz", ".bgz")):
        return gzip.open(p, "rt", encoding="utf-8", newline="")
    return open(p, "rt", encoding="utf-8", newline="")


def _int(value: str) -> Optional[int]:
    try:
        return int(value)
    except ValueError:
        return None


def _float(value: str) -> Optional[float]:
    try:
        return float(value)
    except ValueError:
        return None


def _nth(value: Optional[str], i: int) -> Optional[str]:
    """The i-th (0-based) comma-separated element, or the value itself if it has none."""
    if value is None:
        return None
    parts = value.split(",")
    return parts[i] if i < len(parts) else None


def iter_calls(path: str, stats: Optional[ReadStats] = None) -> Iterator[Call]:
    stats = stats if stats is not None else ReadStats()
    fmt_cache: Dict[str, Dict[str, int]] = {}
    with open_text(path) as f:
        for line in f:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    cols = line.rstrip("\n").split("\t")
                    stats.samples = max(0, len(cols) - 9)
                    stats.sample_name = cols[9] if len(cols) > 9 else None
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            stats.records += 1
            chrom, pos_s, _id, ref, alt_s, qual_s, filt, _info, fmt, sample = fields[:10]
            keys = fmt_cache.get(fmt)
            if keys is None:
                keys = fmt_cache[fmt] = {k: i for i, k in enumerate(fmt.split(":"))}
            if "GT" not in keys:
                stats.no_alt_call += 1
                continue
            values = sample.split(":")
            gt = values[keys["GT"]] if keys["GT"] < len(values) else "."
            alleles = _GT_SPLIT.split(gt)
            if "." in alleles or all(a == "0" for a in alleles):
                stats.no_alt_call += 1
                continue
            phased = "|" in gt

            def field(name: str) -> Optional[str]:
                i = keys.get(name)
                return values[i] if i is not None and i < len(values) and values[i] != "." else None

            ps = field("PS")
            dp = _int(field("DP") or "")
            gq = _int(field("GQ") or "")
            vaf_field = field("VAF")
            qual = _float(qual_s) if qual_s != "." else None
            pos = int(pos_s)
            chrom = canon_chrom(chrom)
            alts: List[str] = alt_s.split(",")
            for idx in sorted({int(a) for a in alleles if a != "0"}):
                if idx > len(alts):
                    continue
                alt = alts[idx - 1]
                if _SYMBOLIC_ALT.match(alt):
                    stats.symbolic_or_star += 1
                    continue
                copies = sum(1 for a in alleles if a == str(idx))
                haps = "".join(str(h + 1) for h, a in enumerate(alleles) if a == str(idx)) if phased else ""
                npos, nref, nalt = trim_alleles(pos, ref.upper(), alt.upper())
                vaf = _float(_nth(vaf_field, idx - 1) or "")
                yield Call(
                    chrom, pos, ref, alt, variant_kind(nref, nalt), copies, gt, phased, ps, haps,
                    qual, filt, dp, gq, vaf, npos, nref, nalt,
                )
