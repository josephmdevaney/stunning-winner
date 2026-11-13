"""Which gene body a variant falls in, from the same RefSeq gene-body BED RepeatScope uses
(``project/uploads/genes.hg38.bed``: chrom, start, stop, symbol, strand; one row per gene
locus). Overlap, not consequence: a variant in a gene body can be intronic or intergenic
to the coding sequence, and this index does not tell them apart."""

from __future__ import annotations

import bisect
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

from .normalize import canon_chrom


class GeneIndex:
    def __init__(self, rows: List[Tuple[str, int, int, str]]):
        by: Dict[str, List[Tuple[int, int, str]]] = defaultdict(list)
        for chrom, start, stop, symbol in rows:
            by[canon_chrom(chrom)].append((start, stop, symbol))
        self._c: Dict[str, Tuple[List[int], List[int], List[str], List[int]]] = {}
        for chrom, items in by.items():
            items.sort(key=lambda x: x[0])
            ends = [x[1] for x in items]
            prefix_max, m = [], -1
            for e in ends:
                m = max(m, e)
                prefix_max.append(m)
            self._c[chrom] = ([x[0] for x in items], ends, [x[2] for x in items], prefix_max)
        self.n_genes = len(rows)

    @classmethod
    def from_bed(cls, path: str | Path) -> "GeneIndex":
        rows = []
        with open(path, encoding="utf-8") as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 4:
                    rows.append((p[0], int(p[1]), int(p[2]), p[3]))
        return cls(rows)

    def genes_at(self, chrom: str, start0: int, end0: int) -> Tuple[str, ...]:
        """Symbols of every gene body overlapping [start0, end0), in genomic order."""
        entry = self._c.get(chrom)
        if entry is None:
            return ()
        starts, ends, names, prefix_max = entry
        end0 = max(end0, start0 + 1)
        j = bisect.bisect_left(starts, end0) - 1
        found = []
        while j >= 0 and prefix_max[j] > start0:
            if ends[j] > start0:
                found.append(names[j])
            j -= 1
        return tuple(reversed(found))
