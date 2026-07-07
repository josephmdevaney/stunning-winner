"""Loader for the genome-wide repeat catalog BED (TRGT's `--repeats` input format).

Format: tab-separated `chrom  start  stop  ID=...;MOTIFS=...;STRUC=...`, one row per
repeat region -- the same catalog file used to genotype the sample with TRGT, and the
same file the original design's Repeat Catalog screen parsed client-side (see
`project/RepeatCatalogScreen.dc.html`, `project/uploads/repeat_catalog.hg38.bed`).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator


@dataclass(frozen=True)
class CatalogRegion:
    chrom: str
    start: int
    stop: int
    id: str
    motifs: str
    struc: str

    @property
    def length(self) -> int:
        return self.stop - self.start


def iter_catalog_regions(bed_path: str | Path) -> Iterator[CatalogRegion]:
    with open(bed_path, encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, start, stop, info = parts[0], parts[1], parts[2], parts[3]
            try:
                start_i, stop_i = int(start), int(stop)
            except ValueError:
                continue
            region_id = ""
            motifs = ""
            struc = ""
            for kv in info.split(";"):
                if "=" not in kv:
                    continue
                k, v = kv.split("=", 1)
                if k == "ID":
                    region_id = v
                elif k == "MOTIFS":
                    motifs = v
                elif k == "STRUC":
                    struc = v
            yield CatalogRegion(
                chrom=chrom, start=start_i, stop=stop_i, id=region_id, motifs=motifs, struc=struc
            )
