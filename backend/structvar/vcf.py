"""Streaming reader for structural-variant VCFs (sawfish, pbsv, sniffles and the like).

One record is one call. Each caller writes the types a little differently, so this reader
normalises what the review screens need:

  * INS: the full inserted sequence in ALT (sawfish, pbsv) or a symbolic <INS>; size from SVLEN or from
    the allele lengths.
  * DEL: a long REF with a one-base ALT (sawfish, pbsv) or a symbolic <DEL> (imprecise depth-based events).
  * DUP / DUP:TANDEM / INV / CNV: symbolic; size from SVLEN, else END - POS.
  * BND: a breakend; the mate chromosome and position are read from ALT. How a breakend relates to an
    inversion or translocation is left to the caller's EVENT/MATEID fields, which are kept.

A genotype is kept per sample. ``copies`` is the number of ALT alleles the genotype carries: 0 reference,
1 heterozygous, 2 homozygous, and -1 for a no-call. A half-called genotype such as ``./1`` still carries
the variant (copies 1, ``partial`` set). A record whose REF or ALT is not valid VCF (one pbsv insertion at a contig end in the data this was built on
has a NUL byte for REF) is skipped and counted, not read as if it were sound. Support comes from whichever fields the caller wrote: AD (ref,alt),
DV/DR (sniffles), DP, GQ, and CN for copy-number events.
"""

from __future__ import annotations

import gzip
import re
from typing import Dict, Iterator, List, NamedTuple, Optional, Tuple

from .normalize import canon_chrom

PASS_FILTERS = ("PASS", ".")
_GT_SPLIT = re.compile(r"[/|]")
_BND = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")
_SYMBOLIC = re.compile(r"^<([A-Z]+)(?::[A-Za-z]+)?>$")
_REF_OK = re.compile(r"^[ACGTNRYSWKMBDHVacgtnryswkmbdhv]+$")
_ALT_OK = re.compile(r"^(?:[ACGTNRYSWKMBDHVacgtnryswkmbdhv]+|<[^<>\s]+>|\*|.*[\[\]].*)$")


class Genotype(NamedTuple):
    gt: str
    copies: int
    partial: bool
    gq: Optional[int]
    ad_ref: Optional[int]
    ad_alt: Optional[int]
    dp: Optional[int]
    cn: Optional[float]


class SvCall(NamedTuple):
    vid: str
    chrom: str
    pos: int
    end: int
    svtype: str  # DEL | DUP | INS | INV | CNV | BND
    size: Optional[int]
    filter: str
    passed: bool
    qual: Optional[float]
    imprecise: bool
    mate_chrom: Optional[str]
    mate_pos: Optional[int]
    event: Optional[str]
    symbolic: bool  # ALT is <TYPE> rather than a sequence
    samples: Tuple[Genotype, ...]


class ReadStats:
    def __init__(self) -> None:
        self.records = 0
        self.skipped = 0  # records of a type this reader does not handle
        self.invalid = 0  # records with an allele that is not valid VCF (e.g. a stray NUL byte); skipped
        self.sample_names: List[str] = []
        self.source: Optional[str] = None  # the caller, from the ##source header


def open_text(path: str):
    """Plain or gzip/bgzip, decided by the first bytes (a name ending in .gz is not proof)."""
    with open(path, "rb") as f:
        magic = f.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "rt", encoding="utf-8", newline="")


def _int(v: Optional[str]) -> Optional[int]:
    try:
        return int(float(v)) if v not in (None, ".", "") else None
    except ValueError:
        return None


def _float(v: Optional[str]) -> Optional[float]:
    try:
        return float(v) if v not in (None, ".", "") else None
    except ValueError:
        return None


def _info(field: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for kv in field.split(";"):
        k, _, v = kv.partition("=")
        out[k] = v
    return out


def parse_genotype(sample: str, keys: Dict[str, int]) -> Genotype:
    vals = sample.split(":")

    def get(name: str) -> Optional[str]:
        i = keys.get(name)
        return vals[i] if i is not None and i < len(vals) and vals[i] not in (".", "") else None

    gt = get("GT") or "./."
    alleles = _GT_SPLIT.split(gt)
    alt = sum(1 for a in alleles if a not in ("0", "."))
    called = [a for a in alleles if a != "."]
    if not called:
        copies, partial = -1, False
    else:
        copies, partial = alt, "." in alleles and alt > 0
    ad_ref = ad_alt = None
    ad = get("AD")
    if ad:
        parts = ad.split(",")
        ad_ref = _int(parts[0])
        ad_alt = _int(parts[1]) if len(parts) > 1 else None
    if ad_alt is None and get("DV") is not None:  # sniffles
        ad_alt, ad_ref = _int(get("DV")), _int(get("DR"))
    dp = _int(get("DP"))
    if dp is None and ad_ref is not None and ad_alt is not None:
        dp = ad_ref + ad_alt
    return Genotype(gt, copies, partial, _int(get("GQ")) if get("GQ") is not None else _int(get("CNQ")), ad_ref, ad_alt, dp, _float(get("CN")))


def iter_calls(path: str, stats: Optional[ReadStats] = None) -> Iterator[SvCall]:
    stats = stats if stats is not None else ReadStats()
    fmt_cache: Dict[str, Dict[str, int]] = {}
    with open_text(path) as f:
        for line in f:
            if line.startswith("#"):
                if line.startswith("##source="):
                    stats.source = line.split("=", 1)[1].strip().strip('"')
                elif line.startswith("#CHROM"):
                    stats.sample_names = line.rstrip("\n").split("\t")[9:]
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 8:
                continue
            stats.records += 1
            chrom, pos_s, vid, ref, alt, qual_s, flt, inf = p[:8]
            if not _REF_OK.match(ref) or not _ALT_OK.match(alt.split(",")[0]):
                stats.invalid += 1
                continue
            i = _info(inf)
            pos = int(pos_s)
            m = _SYMBOLIC.match(alt)
            svtype = i.get("SVTYPE") or (m.group(1) if m else "")
            if alt.startswith("<") and ":" in alt and not svtype:
                svtype = alt[1:].split(":")[0]
            mate_chrom = mate_pos = None
            bnd = _BND.search(alt)
            if bnd:
                svtype, mate_chrom, mate_pos = "BND", canon_chrom(bnd.group(1)), int(bnd.group(2))
            if svtype not in ("DEL", "DUP", "INS", "INV", "CNV", "BND"):
                stats.skipped += 1
                continue
            svlen = abs(_int(i.get("SVLEN", "").split(",")[0]) or 0) or None
            end = _int(i.get("END"))
            if svtype == "INS":
                end = pos
                size = svlen or (len(alt) - len(ref) if not m and len(alt) > len(ref) else None)
            elif svtype == "BND":
                end = pos
                size = None
            else:
                if end is None:
                    end = pos + (svlen or max(len(ref) - 1, 0))
                size = svlen or max(end - pos, 1)
            samples: Tuple[Genotype, ...] = ()
            if len(p) > 9:
                keys = fmt_cache.get(p[8])
                if keys is None:
                    keys = fmt_cache[p[8]] = {k: n for n, k in enumerate(p[8].split(":"))}
                samples = tuple(parse_genotype(s, keys) for s in p[9:])
            yield SvCall(
                vid, canon_chrom(chrom), pos, end, svtype, size, flt, flt in PASS_FILTERS, _float(qual_s), "IMPRECISE" in i,
                mate_chrom, mate_pos, i.get("EVENT") or None, bool(m), samples,
            )
