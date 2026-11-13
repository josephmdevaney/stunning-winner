"""Genome tracks that mark regions where SV calls are less trustworthy or are expected to be common:
segmental duplications, repeats (RepeatMasker and tandem repeats), and "odd" regions (gaps, telomeres,
heterochromatin, PARs, the immunoglobulin V(D)J loci), from PacBio's HiFi WDL resource bundle
(``static_resources/GRCh38/annotation/``, UCSC tables).

Repeats are queried with ``tabix`` (the repeat file is far too large to hold), and only for SVs up to
100 kb, since a repeat label means little for a megabase event. If ``tabix`` is missing the repeat
columns stay empty and the suite says so."""

from __future__ import annotations

import bisect
import gzip
import os
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from .normalize import canon_chrom

MAX_REPEAT_SV = 100_000


def _merge(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    out: List[Tuple[int, int]] = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


class Regions:
    """Disjoint intervals per chromosome, with the fraction of a span they cover."""

    def __init__(self, by_chrom: Dict[str, List[Tuple[int, int]]]):
        self._c = {}
        for chrom, iv in by_chrom.items():
            m = _merge(iv)
            self._c[chrom] = ([s for s, _ in m], [e for _, e in m])

    def fraction(self, chrom: str, start: int, end: int) -> float:
        entry = self._c.get(chrom)
        if entry is None:
            return 0.0
        starts, ends = entry
        length = max(end - start + 1, 1)
        j = max(bisect.bisect_right(starts, end) - 1, 0)
        covered = 0
        while j >= 0 and ends[j] >= start:
            covered += max(0, min(end, ends[j]) - max(start, starts[j]) + 1)
            j -= 1
        return min(covered / length, 1.0)


def _read_bed(path: Path, named: bool = False):
    with gzip.open(path, "rt", encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3:
                yield canon_chrom(p[0]), int(p[1]) + 1, int(p[2]), (p[3] if len(p) > 3 else "")


class Tracks:
    def __init__(self, bundle: str | Path):
        d = Path(bundle).expanduser() / "GRCh38" / "annotation"
        self.dir = d
        self.segdups: Optional[Regions] = None
        self.odd: Dict[str, Regions] = {}
        seg, odd = d / "GRCh38.segdups.bed.gz", d / "GRCh38.oddRegions.bed.gz"
        if seg.is_file():
            by = defaultdict(list)
            for c, s, e, _ in _read_bed(seg):
                by[c].append((s, e))
            self.segdups = Regions(by)
        if odd.is_file():
            names = defaultdict(lambda: defaultdict(list))
            for c, s, e, n in _read_bed(odd, True):
                names[n or "odd"][c].append((s, e))
            self.odd = {n: Regions(by) for n, by in names.items()}
        self.repeats = d / "GRCh38.repeats.bed.gz"
        self.has_repeats = self.repeats.is_file() and shutil.which("tabix") is not None
        self.stamp = self.stamp_for(bundle)

    @staticmethod
    def stamp_for(bundle: str | Path) -> str:
        """What a rebuild is keyed on: which track files exist (with sizes), and whether tabix does."""
        d = Path(bundle).expanduser() / "GRCh38" / "annotation"
        files = [d / "GRCh38.segdups.bed.gz", d / "GRCh38.oddRegions.bed.gz", d / "GRCh38.repeats.bed.gz"]
        parts = [p for p in files if p.is_file()]
        has_repeats = files[2].is_file() and shutil.which("tabix") is not None
        return ";".join(f"{p.name}:{p.stat().st_size}" for p in parts) + (";tabix" if has_repeats else "")

    def odd_regions(self, chrom: str, start: int, end: int) -> List[str]:
        return sorted(n for n, r in self.odd.items() if r.fraction(chrom, start, end) > 0)

    def repeat_labels(self, spans: Iterable[Tuple[int, str, int, int]]) -> Dict[int, Tuple[str, float]]:
        """For (key, chrom, start, end) spans, the repeat class covering most of each and how much of it."""
        spans = [s for s in spans if s[3] - s[2] + 1 <= MAX_REPEAT_SV]
        if not spans or not self.has_repeats:
            return {}
        regions = sorted((c, s, e) for _, c, s, e in spans)
        with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as t:
            for c, s, e in regions:
                t.write(f"{c}\t{max(s - 1, 0)}\t{e + 1}\n")
            bed = t.name
        try:
            proc = subprocess.run(["tabix", "-R", bed, str(self.repeats)], capture_output=True, text=True, timeout=600)
        finally:
            os.unlink(bed)
        by: Dict[str, Dict[str, List[Tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
        for line in proc.stdout.splitlines():
            p = line.split("\t")
            if len(p) < 4:
                continue
            by[canon_chrom(p[0])][_repeat_class(p[3])].append((int(p[1]) + 1, int(p[2])))
        merged = {c: {k: Regions({c: iv}) for k, iv in d.items()} for c, d in by.items()}
        out: Dict[int, Tuple[str, float]] = {}
        for key, c, s, e in spans:
            best = max(((lab, r.fraction(c, s, e)) for lab, r in merged.get(c, {}).items()), key=lambda x: x[1], default=None)
            if best and best[1] > 0:
                out[key] = (best[0], round(best[1], 3))
        return out


def _repeat_class(name: str) -> str:
    """A short label from a repeats.bed name: ``rmsk:<class>:<name>`` or ``trf:<unit>``."""
    parts = name.split(":")
    if parts[0] == "trf":
        return "Tandem repeat"
    cls, rep = (parts[1] if len(parts) > 1 else ""), (parts[2] if len(parts) > 2 else "")
    if cls.startswith("SINE") and rep.startswith("Alu"):
        return "Alu"
    if cls.startswith("LINE") and rep.startswith("L1"):
        return "L1"
    if cls.startswith("LTR"):
        return "LTR"
    if cls.startswith(("Simple_repeat", "Low_complexity")):
        return "Tandem repeat"
    if cls.startswith("Satellite"):
        return "Satellite"
    return "Other repeat"
