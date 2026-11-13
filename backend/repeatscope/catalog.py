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


def chrom_sort_key(chrom: str):
    order = {f"chr{i}": i for i in range(1, 23)}
    order["chrX"] = 23
    order["chrY"] = 24
    order["chrM"] = 25
    return order.get(chrom, 99)


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


class CatalogIndex:
    """Chromosome-bucketed, start-sorted index of the *entire* genome-wide
    catalog, held as lightweight CatalogRegion values only -- no genotype/gene/
    pathogenicity annotation. Building the index this way, and deferring the
    comparatively expensive per-region annotation (TRGT genotype lookup, gene
    overlap, RExPRT score) to only the rows a request actually needs (e.g. one
    paginated page of ~50, not all of them), is what makes serving a
    multi-million-row catalog (TRExplorer v2.0: 5,599,658 regions) practical --
    eagerly annotating every row up front measured at ~4 minutes and 15GB+ RAM
    at that scale, vs. ~14s/1.8GB just to build this raw index.

    Callers needing a genuine full-catalog scan (e.g. Run QC's aggregate
    genotyped/depth stats) still pay an O(n) annotation cost, same as before --
    this index only removes the *redundant* re-parsing of the BED file between
    that use and the paginated /api/catalog route, and lets the latter skip
    annotating anything outside the requested page.
    """

    def __init__(self, regions: list[CatalogRegion]):
        self._by_chrom: dict[str, list[CatalogRegion]] = {}
        motifs_seen: set[str] = set()
        for r in regions:
            self._by_chrom.setdefault(r.chrom, []).append(r)
            motifs_seen.add(r.motifs)
        for chrom_regions in self._by_chrom.values():
            chrom_regions.sort(key=lambda r: r.start)
        self.chromosomes: list[str] = sorted(self._by_chrom.keys(), key=chrom_sort_key)
        # Flattened once, in canonical chromosome order, so full-scan consumers
        # (Run QC) don't need their own separate concatenation/sort logic.
        self.all: list[CatalogRegion] = [r for chrom in self.chromosomes for r in self._by_chrom[chrom]]
        self.total: int = len(self.all)
        self.unique_motifs: int = len(motifs_seen)

    def filter(self, chrom: str = "all", motif: str = "", min_len: int = 0) -> list[CatalogRegion]:
        rows = self._by_chrom.get(chrom, []) if chrom != "all" else self.all
        if motif:
            needle = motif.upper()
            rows = [r for r in rows if needle in r.motifs.upper()]
        if min_len:
            rows = [r for r in rows if r.length >= min_len]
        return rows

    @classmethod
    def from_bed(cls, bed_path: str | Path) -> "CatalogIndex":
        return cls(list(iter_catalog_regions(bed_path)))
