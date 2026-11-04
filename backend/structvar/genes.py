"""Gene models for structural-variant impact, from the table built by
``scripts/real_data/build_sv_genes.py``: the MANE Select transcript (UCSC's MANE track) where a gene has
one, otherwise the Ensembl 101 transcript with the longest coding sequence.

An SV is compared with each gene it touches and described by what it does to that transcript's
exons: the whole gene, one or more numbered exons (coding or not), or only introns. This
is overlap with the gene model, not a predicted consequence: whether a deletion of exons 3-5 keeps the
reading frame is not worked out here."""

from __future__ import annotations

import bisect
import gzip
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Tuple

from .normalize import canon_chrom

# severity of a gene hit; the largest across an SV's genes is the SV's own impact
IMPACT_RANK = {"whole_gene": 5, "coding_exon": 4, "noncoding_exon": 3, "intronic": 2, "none": 0}
IMPACT_LABEL = {
    "whole_gene": "Whole gene", "coding_exon": "Coding exon", "noncoding_exon": "Non-coding exon / UTR",
    "intronic": "Intronic", "none": "No gene",
}


class Gene(NamedTuple):
    gene_id: str
    symbol: str
    biotype: str
    chrom: str
    start: int
    end: int
    strand: str
    transcript: str  # the display accession: a RefSeq NM_ for MANE Select, else the Ensembl ENST
    ensembl_transcript: str
    model_source: str  # "MANE Select" or the Ensembl 101 fallback
    exons: Tuple[Tuple[int, int], ...]  # genomic order
    cds_start: Optional[int]
    cds_end: Optional[int]
    pli: Optional[float]
    oe_lof: Optional[float]
    inheritance: Tuple[str, ...]
    hpo_n: int
    hpo_terms: Tuple[str, ...]
    diseases: Tuple[str, ...]


def _f(v: str) -> Optional[float]:
    try:
        return float(v)
    except ValueError:
        return None


class GeneModels:
    def __init__(self, genes: List[Gene], stamp: str = ""):
        self.stamp = stamp
        self.n_genes = len(genes)
        self.by_symbol: Dict[str, Gene] = {}
        by_chrom: Dict[str, List[Gene]] = defaultdict(list)
        for g in genes:
            by_chrom[g.chrom].append(g)
            if g.symbol not in self.by_symbol or (g.biotype == "protein_coding") > (self.by_symbol[g.symbol].biotype == "protein_coding"):
                self.by_symbol[g.symbol] = g
        self._c: Dict[str, Tuple[List[int], List[int], List[Gene]]] = {}
        for chrom, items in by_chrom.items():
            items.sort(key=lambda g: g.start)
            m, prefix = -1, []
            for g in items:
                m = max(m, g.end)
                prefix.append(m)
            self._c[chrom] = ([g.start for g in items], prefix, items)

    @classmethod
    def load(cls, path: str | Path) -> "GeneModels":
        p = Path(path)
        genes: List[Gene] = []
        with gzip.open(p, "rt", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#") or line.startswith("gene_id"):
                    continue
                v = line.rstrip("\n").split("\t")
                exons = tuple(tuple(int(x) for x in e.split("-")) for e in v[10].split(",") if e)  # type: ignore[misc]
                genes.append(Gene(
                    v[0], v[1], v[2], v[3], int(v[4]), int(v[5]), v[6], v[7], v[8], v[9], exons,
                    int(v[11]) if v[11] else None, int(v[12]) if v[12] else None, _f(v[13]), _f(v[14]),
                    tuple(x for x in v[15].split(",") if x), int(v[16] or 0),
                    tuple(x for x in v[17].split("; ") if x), tuple(x for x in v[18].split("; ") if x),
                ))
        st = p.stat()
        return cls(genes, f"{p.name}:{st.st_size}:{st.st_mtime_ns}")

    def overlapping(self, chrom: str, start: int, end: int) -> List[Gene]:
        entry = self._c.get(canon_chrom(chrom))
        if entry is None:
            return []
        starts, prefix, items = entry
        j = bisect.bisect_right(starts, end) - 1
        out = []
        while j >= 0 and prefix[j] >= start:
            if items[j].end >= start:
                out.append(items[j])
            j -= 1
        out.reverse()
        return out


def _exon_numbers(g: Gene, hits: List[int]) -> List[int]:
    n = len(g.exons)
    return sorted(n - i if g.strand == "-" else i + 1 for i in hits)


def _coding(g: Gene, exon: Tuple[int, int]) -> bool:
    return g.cds_start is not None and g.cds_end is not None and exon[1] >= g.cds_start and exon[0] <= g.cds_end


def gene_hit(svtype: str, start: int, end: int, g: Gene) -> Dict[str, object]:
    """How one SV touches one gene. For INS/BND ``start == end`` is the breakpoint."""
    point = svtype in ("INS", "BND")
    if not point and svtype != "INV" and start <= g.start and end >= g.end:
        relation, hits = "whole_gene", [i for i in range(len(g.exons))]
    elif point or svtype == "INV":
        pts = [start] if point else [start, end]
        hits = sorted({i for i, e in enumerate(g.exons) for p in pts if e[0] <= p <= e[1]})
        inside = any(g.start <= p <= g.end for p in pts)
        if not inside and svtype == "INV" and start <= g.start and end >= g.end:
            relation, hits = "none", []  # the gene lies wholly inside the inverted segment
        elif hits:
            relation = "coding_exon" if any(_coding(g, g.exons[i]) for i in hits) else "noncoding_exon"
        else:
            relation = "intronic" if inside else "none"
    else:
        hits = [i for i, e in enumerate(g.exons) if e[1] >= start and e[0] <= end]
        if hits:
            relation = "coding_exon" if any(_coding(g, g.exons[i]) for i in hits) else "noncoding_exon"
        else:
            relation = "intronic"
    return {
        "symbol": g.symbol, "gene_id": g.gene_id, "biotype": g.biotype, "relation": relation,
        "exons": _exon_numbers(g, hits) if relation != "whole_gene" else [], "n_exons": len(g.exons), "strand": g.strand,
        "transcript": g.transcript, "model_source": g.model_source, "rank": IMPACT_RANK[relation],
    }


def impacts(models: GeneModels, svtype: str, chrom: str, start: int, end: int,
            mate: Optional[Tuple[str, int]] = None) -> List[Dict[str, object]]:
    """Every gene the SV touches, most severe first (protein-coding before the rest on ties)."""
    spans = [(chrom, start, end)]
    if svtype == "BND":
        spans = [(chrom, start, start)] + ([(mate[0], mate[1], mate[1])] if mate else [])
    seen, out = set(), []
    for c, s, e in spans:
        probe_s, probe_e = (s, e) if svtype not in ("INV",) else (s, e)
        for g in models.overlapping(c, probe_s, probe_e):
            if g.gene_id in seen:
                continue
            h = gene_hit(svtype, s, e, g)
            if h["relation"] == "none":
                continue
            seen.add(g.gene_id)
            out.append(h)
    out.sort(key=lambda h: (-int(h["rank"]), h["biotype"] != "protein_coding", str(h["symbol"])))  # type: ignore[arg-type]
    return out
