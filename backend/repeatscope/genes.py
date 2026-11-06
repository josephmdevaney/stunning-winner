"""Nearest/overlapping-gene annotation for genome-wide repeat catalog regions.

Loads a gene-body reference BED (chrom, start, stop, symbol, strand -- see
`project/uploads/genes.hg38.bed`, derived from UCSC's RefSeq `refGene` table,
transcripts collapsed to one row per real gene locus; see that file's own
provenance note in scripts/real_data if present) and answers, for any genomic
interval: which real gene(s) it falls inside, or -- if none -- the single
nearest gene, its distance, and which side (upstream/downstream, relative to
that gene's own strand, not raw genomic left/right).
"""

from __future__ import annotations

import bisect
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class GeneRegion:
    chrom: str
    start: int
    stop: int
    symbol: str
    strand: str


def load_genes(bed_path: str | Path) -> list[GeneRegion]:
    genes: list[GeneRegion] = []
    with open(bed_path, encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            chrom, start, stop, symbol, strand = line.rstrip("\n").split("\t")
            genes.append(GeneRegion(chrom, int(start), int(stop), symbol, strand))
    return genes


@dataclass(frozen=True)
class GeneAnnotation:
    """Either `genes` is non-empty (the region falls inside one or more real gene
    bodies -- rare but real, e.g. antisense/overlapping transcripts) and
    `nearest_gene` is None, or the region is intergenic: `genes` is empty and
    `nearest_gene`/`nearest_distance_bp`/`nearest_direction` describe the single
    closest gene. Direction is relative to that gene's own strand ("upstream" of
    a minus-strand gene is the *higher* genomic coordinate, not the lower one),
    since that's what "upstream/downstream" means biologically -- a raw
    left/right label would be actively misleading for half of all genes."""

    genes: tuple[str, ...]
    nearest_gene: str | None = None
    nearest_distance_bp: int | None = None
    nearest_direction: str | None = None  # "upstream" | "downstream"


def _direction_for(gene: GeneRegion, region_is_left_of_gene: bool) -> str:
    """Which side of `gene` a region sits on, in upstream/downstream terms.

    Plus-strand genes transcribe start->stop, so their TSS/upstream side is the
    lower-coordinate (start) side. Minus-strand genes transcribe stop->start, so
    their TSS/upstream side is the *higher*-coordinate (stop) side instead.
    """
    left_is_upstream = gene.strand == "+"
    if region_is_left_of_gene:
        return "upstream" if left_is_upstream else "downstream"
    return "downstream" if left_is_upstream else "upstream"


class GeneIndex:
    """Chromosome-bucketed, start-sorted index for overlap + nearest-gene lookup
    against a real gene set that can itself contain overlapping/nested entries
    (unlike a repeat catalog's own non-overlapping regions -- antisense pairs and
    readthrough transcripts are real biology, not a data error), so lookups can't
    assume a small fixed backward-scan window is always enough."""

    def __init__(self, genes: list[GeneRegion]):
        self._by_chrom: dict[str, list[GeneRegion]] = {}
        for g in genes:
            self._by_chrom.setdefault(g.chrom, []).append(g)
        for chrom_genes in self._by_chrom.values():
            chrom_genes.sort(key=lambda g: g.start)
        self._starts: dict[str, list[int]] = {
            chrom: [g.start for g in chrom_genes] for chrom, chrom_genes in self._by_chrom.items()
        }
        # The longest real gene body on each chromosome bounds how far back an
        # overlap scan must look: no gene starting earlier than
        # (query_start - max_span) can possibly still reach forward into the
        # query interval, since no gene on this chromosome is longer than that.
        # This bound is only valid for *overlap* detection -- NOT for the
        # "nearest gene before the query" search below, where a gene's start
        # can legitimately sit further back than that (a long gene ending well
        # short of the query, not reaching into it, can still start earlier
        # than query_start - max_span while remaining the closest non
        # -overlapping candidate).
        self._max_span: dict[str, int] = {
            chrom: max((g.stop - g.start for g in chrom_genes), default=0)
            for chrom, chrom_genes in self._by_chrom.items()
        }
        # Prefix maximum of `stop`, in start-sorted order: running_max_end[i] is
        # the gene with the largest end among genes[0..i] (storing the gene
        # itself, not just its end, so a lookup never needs a second search to
        # recover which gene achieved it). Once overlap has been ruled out for a
        # query, every gene with start <= query.start is guaranteed to also have
        # end <= query.start (otherwise it would itself be an overlap) -- so the
        # nearest-before candidate is exactly this prefix max at the index of
        # the last gene with start <= query.start, found in O(1) after an
        # O(log n) bisect, rather than rescanning a window per query.
        self._running_max_end_gene: dict[str, list[GeneRegion]] = {}
        for chrom, chrom_genes in self._by_chrom.items():
            running: list[GeneRegion] = []
            best = chrom_genes[0]
            for g in chrom_genes:
                if g.stop > best.stop:
                    best = g
                running.append(best)
            self._running_max_end_gene[chrom] = running

    def annotate(self, chrom: str, start: int, stop: int) -> GeneAnnotation:
        chrom_genes = self._by_chrom.get(chrom)
        if not chrom_genes:
            return GeneAnnotation(genes=())
        starts = self._starts[chrom]
        max_span = self._max_span[chrom]

        # Overlap requires gene.start < stop AND gene.end > start; given
        # gene.end - gene.start <= max_span, gene.end > start implies
        # gene.start > start - max_span, so only genes in this window need
        # checking. `hi` uses bisect_left (not bisect_right): a gene whose start
        # sits exactly at `stop` is adjacent, not overlapping (BED half-open
        # convention) -- and using bisect_left here also means `chrom_genes[hi]`
        # below is exactly the nearest touching-or-later gene, with the correct
        # zero distance for an exact touch.
        lo = bisect.bisect_left(starts, start - max_span)
        hi = bisect.bisect_left(starts, stop)
        window = chrom_genes[lo:hi]

        overlapping = [g.symbol for g in window if g.stop > start]
        if overlapping:
            # Preserve genomic order, dedupe (a gene can only appear once in
            # `window`, but be defensive rather than assume the source BED
            # never has an exact duplicate row).
            seen: list[str] = []
            for symbol in overlapping:
                if symbol not in seen:
                    seen.append(symbol)
            return GeneAnnotation(genes=tuple(seen))

        # No overlap. Nearest-before: the last gene with start <= start (there is
        # one iff hi > 0, since `hi` is the first index with start >= stop >
        # start) is guaranteed non-overlapping here, and the prefix-max-end
        # array gives its true nearest-before candidate's end in O(1) --
        # correctly picking the gene with the largest end even if an earlier
        # -starting, longer gene reaches closer than a later-starting, shorter
        # one. Nearest-after: `chrom_genes[hi]`, the first gene starting at or
        # after `stop`, directly (no equivalent subtlety on this side, since
        # among non-overlapping "after" genes the one starting soonest is
        # necessarily the closest).
        nearest_before = self._running_max_end_gene[chrom][hi - 1] if hi > 0 else None
        nearest_after = chrom_genes[hi] if hi < len(chrom_genes) else None

        # `nearest_before`'s gene sits to the *left* of the region (gene.stop <=
        # region.start) -- so the region itself is to the gene's right, i.e.
        # `region_is_left_of_gene=False`. `nearest_after`'s gene sits to the
        # region's right, so the region is to that gene's left: `True`.
        candidates = []
        if nearest_before is not None:
            candidates.append((start - nearest_before.stop, nearest_before, False))
        if nearest_after is not None:
            candidates.append((nearest_after.start - stop, nearest_after, True))
        if not candidates:
            return GeneAnnotation(genes=())

        distance, gene, region_is_left_of_gene = min(candidates, key=lambda c: c[0])
        return GeneAnnotation(
            genes=(),
            nearest_gene=gene.symbol,
            nearest_distance_bp=distance,
            nearest_direction=_direction_for(gene, region_is_left_of_gene),
        )
