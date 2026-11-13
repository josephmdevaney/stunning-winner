"""Population SV catalogs and the matcher that compares a call with them.

The catalogs are the reduced table from ``scripts/real_data/build_sv_catalogs.py``. SVs are not matched
by exact allele, because two callers (or two people) place the same event a few bases apart, so a match is
an approximation by breakpoints and size, the way svpack and Truvari match SV sets:

  * DEL / DUP / INV / CNV: reciprocal overlap >= 50% (overlap divided by the longer of the two).
  * INS: within 500 bp and the shorter insertion at least half the length of the longer.

Insertion sequences are not compared, so two different insertions at one locus can match. DEL and DUP
can also match a catalog CNV (whose direction the source does not state); that is reported as the
matched type. A source with no frequency (DECODE, EEE) reports presence only.
"""

from __future__ import annotations

import bisect
import gzip
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Tuple

SOURCES = ("gnomad_sv", "hprc", "decode", "eee")
SOURCE_LABEL = {
    "gnomad_sv": "gnomAD-SV (short-read)", "hprc": "HPRC (long-read, pbsv)", "decode": "DECODE (long-read, ONT)",
    "eee": "EEE (long-read discovery)",
}
MIN_RO = 0.5
INS_WINDOW = 500
MIN_SIZE_RATIO = 0.5
_COMPATIBLE = {"DEL": ("DEL", "CNV"), "DUP": ("DUP", "CNV"), "INV": ("INV",), "INS": ("INS",), "CNV": ("DEL", "DUP", "CNV")}


class Match(NamedTuple):
    svtype: str  # the catalog record's type
    similarity: float  # reciprocal overlap, or the insertion size ratio
    n: int  # how many catalog records matched
    af: Optional[float]
    ac: Optional[int]
    an: Optional[int]
    popmax: Optional[float]


def _num(v: str, cast):
    try:
        return cast(v) if v != "" else None
    except ValueError:
        return None


class Catalogs:
    def __init__(self, rows: Dict[Tuple[str, str, str], List[tuple]], counts: Dict[str, int], stamp: str):
        self.stamp = stamp
        self.counts = counts
        # (source, chrom, "int"|"ins") -> parallel lists sorted by start
        self._idx: Dict[Tuple[str, str, str], Tuple[List[int], List[tuple]]] = {}
        for key, items in rows.items():
            items.sort(key=lambda r: r[0])
            self._idx[key] = ([r[0] for r in items], items)

    @classmethod
    def load(cls, path: str | Path) -> "Catalogs":
        p = Path(path)
        rows: Dict[Tuple[str, str, str], List[tuple]] = defaultdict(list)
        counts: Dict[str, int] = defaultdict(int)
        with gzip.open(p, "rt", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#") or line.startswith("source"):
                    continue
                v = line.rstrip("\n").split("\t")
                source, chrom, start, end, t = v[0], v[1], int(v[2]), int(v[3]), v[4]
                svlen = _num(v[5], int)
                rec = (start, end, t, svlen if svlen is not None else max(end - start, 1), _num(v[8], float), _num(v[6], int), _num(v[7], int), _num(v[9], float))
                rows[(source, chrom, "ins" if t == "INS" else "int")].append(rec)
                counts[source] += 1
        st = p.stat()
        return cls(rows, dict(counts), f"{p.name}:{st.st_size}:{st.st_mtime_ns}")

    def match(self, svtype: str, chrom: str, start: int, end: int, size: Optional[int]) -> Dict[str, Match]:
        """Best match per source for a DEL, DUP, INV, INS or CNV call; BND is not matched."""
        want = _COMPATIBLE.get(svtype)
        if want is None or size is None:
            return {}
        out: Dict[str, Match] = {}
        for source in SOURCES:
            if svtype == "INS":
                entry = self._idx.get((source, chrom, "ins"))
                if entry is None:
                    continue
                starts, items = entry
                lo, hi = bisect.bisect_left(starts, start - INS_WINDOW), bisect.bisect_right(starts, start + INS_WINDOW)
                cands = [(min(size, r[3]) / max(size, r[3]), r) for r in items[lo:hi] if max(size, r[3]) and min(size, r[3]) / max(size, r[3]) >= MIN_SIZE_RATIO]
            else:
                entry = self._idx.get((source, chrom, "int"))
                if entry is None:
                    continue
                starts, items = entry
                length = max(end - start + 1, 1)
                lo, hi = bisect.bisect_left(starts, start - length), bisect.bisect_right(starts, start + length)
                cands = []
                for r in items[lo:hi]:
                    if r[2] not in want:
                        continue
                    overlap = min(end, r[1]) - max(start, r[0]) + 1
                    if overlap <= 0:
                        continue
                    ro = overlap / max(length, r[1] - r[0] + 1)
                    if ro >= MIN_RO:
                        cands.append((ro, r))
            if cands:
                sim, best = max(cands, key=lambda c: (c[0], c[1][4] or 0))
                af_values = [c[1][4] for c in cands if c[1][4] is not None]
                out[source] = Match(best[2], round(sim, 3), len(cands), best[4], best[5], best[6], best[7])
                if af_values and best[4] is not None:  # report the highest frequency among the matches
                    top = max(cands, key=lambda c: c[1][4] if c[1][4] is not None else -1)[1]
                    out[source] = Match(top[2], round(sim, 3), len(cands), top[4], top[5], top[6], top[7])
        return out
