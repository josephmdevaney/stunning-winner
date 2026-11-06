#!/usr/bin/env python3
"""Join RExPRT's real genome-wide tandem-repeat pathogenicity predictions onto
RepeatScope's genome-wide hg38 repeat catalog.

RExPRT (Danzi et al. 2024, https://pmc.ncbi.nlm.nih.gov/articles/PMC10832122/,
code+data at https://github.com/ZuchnerLab/RExPRT) is a published ensemble
SVM+XGBoost model that scores tandem repeats for pathogenicity potential using
real genomic-context features: GERP conservation, pLI/LOEUF gene constraint,
exon/UTR/promoter location, transcription-factor-binding-site overlap, TAD
boundaries, eSTR (expression-associated) status, nervous-system expression,
GC content, and Alu proximity. The lab publishes real precomputed scores for
~836,546 reference TRs (motif length 3-8bp) on hg19 coordinates in their repo's
"Reference TR scores/reference_TRsAnnotated_RExPRTscores.txt" (Git LFS, 333MB).

This script lifts those hg19 coordinates to hg38 (UCSC hg19ToHg38 chain, via
pyliftover) and matches them by coordinate overlap against RepeatScope's own
genome-wide catalog (project/uploads/repeat_catalog.hg38.bed -- TRExplorer
v2.0, 5,599,658 regions, real overlapping regions in places, see
project/uploads/README.md), emitting one row per successfully matched catalog
region.

Coverage is intentionally partial and reported to stderr, not papered over:
only a fraction of the catalog's regions have a 3-8bp motif (RExPRT's trained
range) -- the rest are longer-motif repeats RExPRT was never applied to and are
left unmatched. Regions within that range that still don't match are ones
where liftover failed or RExPRT's reference set didn't include that exact
locus. Callers must treat "no match" as "not evaluated by RExPRT", never as
"benign" -- see repeatscope/rexprt.py.

Requires:
    pip install pyliftover
    A local checkout with the LFS blob pulled:
        git clone https://github.com/ZuchnerLab/RExPRT.git
        cd RExPRT && git lfs pull --include="Reference TR scores/reference_TRsAnnotated_RExPRTscores.txt"

Usage:
    python3 build_rexprt_catalog_scores.py \
        "/path/to/RExPRT/Reference TR scores/reference_TRsAnnotated_RExPRTscores.txt" \
        > project/uploads/rexprt_scores.hg38.tsv
"""

from __future__ import annotations

import bisect
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
CATALOG_BED = REPO_ROOT / "project" / "uploads" / "repeat_catalog.hg38.bed"

# Source column name -> our output column name. Values are copied through
# verbatim (only whitespace-stripped); nothing here is recomputed or estimated.
# RExPRT's raw S2SNet topological graph indices (eSh*/eTr*/eH/eW/eS6/eS/J/eX*)
# are real but not clinically interpretable on their own, so they're
# intentionally left out of this join -- the full source file remains the
# reference if anyone needs them.
COLUMN_MAP = [
    ("motif", "motif"),
    ("gene", "gene"),
    ("gene_type", "gene_type"),
    ("region", "region"),
    ("location", "location"),
    ("gerp", "gerp"),
    ("loeuf", "loeuf"),
    ("pLi", "pli"),
    ("promoter", "promoter"),
    ("UTR_5", "utr_5"),
    ("UTR_3", "utr_3"),
    ("TAD", "tad"),
    ("eSTR", "estr"),
    ("opReg", "open_reg_region"),
    ("RAD21", "rad21"),
    ("SMC3", "smc3"),
    ("CTCF", "ctcf"),
    ("alu", "alu"),
    ("alu_distance", "alu_distance"),
    ("gc_content", "gc_content"),
    ("tissue_simple", "tissue"),
    ("SVM", "svm"),
    ("XGB", "xgb"),
    ("ensembleScore", "ensemble_score"),
    ("ensembleMax", "ensemble_max"),
    ("ensembleBinary", "ensemble_binary"),
]


def load_catalog_index():
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    with open(CATALOG_BED, encoding="utf-8") as f:
        for line in f:
            chrom, start, stop, _info = line.rstrip("\n").split("\t")
            by_chrom.setdefault(chrom, []).append((int(start), int(stop)))
    for regions in by_chrom.values():
        regions.sort()
    starts = {chrom: [r[0] for r in regions] for chrom, regions in by_chrom.items()}
    max_spans = {
        chrom: max((r[1] - r[0] for r in regions), default=0) for chrom, regions in by_chrom.items()
    }
    return by_chrom, starts, max_spans


def find_catalog_overlap(by_chrom, starts, max_spans, chrom: str, start: int, stop: int):
    """Best-overlap catalog region for [start, stop) on chrom, if any. Bounds
    the scan window by the chromosome's longest region (max_span) rather than a
    fixed element count: TRExplorer v2.0's catalog is NOT guaranteed
    non-overlapping (real runs of 1000+ consecutive overlapping regions were
    measured directly against the file), where a small fixed window would
    silently miss a real overlap. Same technique as GeneIndex.annotate in
    repeatscope/genes.py and the equivalent fix in TrgtIndex.find_overlap
    (repeatscope/trgt.py)."""
    regions = by_chrom.get(chrom)
    if not regions:
        return None
    chrom_starts = starts[chrom]
    max_span = max_spans[chrom]
    lo = bisect.bisect_left(chrom_starts, start - max_span)
    hi = bisect.bisect_right(chrom_starts, stop)
    best = None
    for j in range(hi - 1, lo - 1, -1):
        r_start, r_stop = regions[j]
        if r_start <= stop and r_stop >= start:
            overlap = min(r_stop, stop) - max(r_start, start)
            if best is None or overlap > best[2]:
                best = (r_start, r_stop, overlap)
    return best[:2] if best else None


def main():
    if len(sys.argv) != 2:
        print(__doc__, file=sys.stderr)
        sys.exit(1)
    src_path = sys.argv[1]

    from pyliftover import LiftOver
    lo = LiftOver("hg19", "hg38")

    by_chrom, starts, max_spans = load_catalog_index()

    n_total = n_lifted = n_matched = n_collisions = 0
    # catalog (chrom, start, stop) -> chosen (overlap_bp, out_row) -- keep the
    # best-overlapping RExPRT source row if more than one lifts onto the same
    # catalog region (rare, but the catalog and RExPRT's reference set weren't
    # built from the identical repeat-finder run).
    best_by_region: dict[tuple[str, int, int], tuple[int, list[str]]] = {}

    with open(src_path, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        col_idx = {name: i for i, name in enumerate(header)}
        missing = [src for src, _ in COLUMN_MAP if src not in col_idx]
        if missing:
            print(f"error: source file missing expected columns: {missing}", file=sys.stderr)
            sys.exit(1)

        for line in f:
            n_total += 1
            fields = line.rstrip("\n").split("\t")
            hg19_chrom = fields[col_idx["chr"]]
            hg19_start = int(fields[col_idx["start"]])
            hg19_end = int(fields[col_idx["end"]])

            lifted_start = lo.convert_coordinate(hg19_chrom, hg19_start)
            lifted_end = lo.convert_coordinate(hg19_chrom, hg19_end)
            if not lifted_start or not lifted_end:
                continue
            hg38_chrom, hg38_start, _, _ = lifted_start[0]
            hg38_chrom_end, hg38_end, _, _ = lifted_end[0]
            if hg38_chrom_end != hg38_chrom or hg38_end < hg38_start:
                continue
            n_lifted += 1

            match = find_catalog_overlap(by_chrom, starts, max_spans, hg38_chrom, hg38_start, hg38_end)
            if match is None:
                continue
            n_matched += 1
            cat_start, cat_stop = match
            overlap = min(cat_stop, hg38_end) - max(cat_start, hg38_start)

            out_row = [fields[col_idx[src]].strip() for src, _ in COLUMN_MAP]
            out_row += [hg19_chrom, str(hg19_start), str(hg19_end)]

            key = (hg38_chrom, cat_start, cat_stop)
            prev = best_by_region.get(key)
            if prev is not None:
                n_collisions += 1
                if prev[0] >= overlap:
                    continue
            best_by_region[key] = (overlap, out_row)

    out_cols = ["chrom", "start", "stop"] + [dst for _, dst in COLUMN_MAP] + ["hg19_chrom", "hg19_start", "hg19_end"]
    print("\t".join(out_cols))
    for (chrom, start, stop), (_overlap, out_row) in sorted(best_by_region.items()):
        print("\t".join([chrom, str(start), str(stop), *out_row]))

    print(
        f"source_rows={n_total} lifted={n_lifted} matched_catalog_regions={len(best_by_region)} "
        f"collisions_resolved={n_collisions}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
