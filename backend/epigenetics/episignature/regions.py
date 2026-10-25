"""Region-averaged methylation for an episignature model's DMR panel.

A model's panel is a bed of single CpG positions, the fourth column naming the DMR each belongs to (many CpGs
share a name). A DMR's value is EpiSignaLR's own: the coverage-weighted pooled fraction

    sum(est_mod_count) / sum(cov)    over the DMR's listed CpGs that have at least `min_coverage` reads

read from the pb-CpG-tools count columns, exactly as `region_average_methylation.py --platform pacbio
--min-coverage 2 --weighted` does. A DMR with no qualifying CpG has no value (the model then substitutes its
training mean).
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from .. import cpg

MIN_COVERAGE = 2        # reads per CpG, as in EpiSignaLR's PacBio command


@dataclass(frozen=True)
class Dmr:
    name: str
    sites: dict[str, tuple[int, ...]]       # chrom -> 0-based CpG starts listed for this DMR

    @property
    def n_listed(self) -> int:
        return sum(len(v) for v in self.sites.values())

    @property
    def span(self) -> tuple[str, int, int]:
        """(chrom, start, end) around the listed CpGs, for display. A DMR is on one chromosome."""
        chrom = next(iter(self.sites))
        positions = self.sites[chrom]
        return chrom, min(positions), max(positions) + 1


@dataclass(frozen=True)
class RegionValue:
    n_cpgs: int                 # listed CpGs that were covered well enough to count
    mean_coverage: float
    fraction: Optional[float]   # 0-1, or None when no CpG qualified


def parse_panel(path: Path) -> dict[str, Dmr]:
    """DMRs in order of first appearance, from a model's `*_full_panel_DMPs.bed`."""
    sites: dict[str, dict[str, list[int]]] = {}
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end, name = line.rstrip("\n").split("\t")[:4]
            for pos in range(int(start), int(end)):
                sites.setdefault(name, {}).setdefault(chrom, []).append(pos)
    return {name: Dmr(name, {c: tuple(p) for c, p in by_chrom.items()}) for name, by_chrom in sites.items()}


def region_values(bed: Path, panel: dict[str, Dmr], min_coverage: int = MIN_COVERAGE) -> dict[str, RegionValue]:
    """Each DMR's pooled methylation from the combined pileup `bed` (one tabix query per DMR and chromosome)."""
    out: dict[str, RegionValue] = {}
    for name, dmr in panel.items():
        mod = cov = n = 0
        for chrom, positions in dmr.sites.items():
            wanted = Counter(positions)     # a position listed twice counts twice, as in the authors' script
            for c in cpg.fetch(bed, chrom, min(positions), max(positions) + 1):
                times = wanted.get(c.pos)
                if times and c.coverage >= min_coverage:
                    mod += c.est_mod * times
                    cov += c.coverage * times
                    n += times
        out[name] = RegionValue(n, cov / n if n else 0.0, mod / cov if cov else None)
    return out
