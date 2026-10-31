#!/usr/bin/env python3
"""Build the sc-eTR expression annotation for the STRchive disease panel.

The TenK10K paper (Tanudisastro et al., bioRxiv 10.1101/2024.11.02.621562, v4) tested
how common repeat-length variation relates to gene expression in 28 immune cell types.
This script joins its released tables to the STRchive disease loci, so the locus screen
can show, for each disease repeat, whether it was polymorphic in that cohort and which
repeats at or near it were the lead expression association for a gene.

Inputs are two files from the paper's Zenodo record (doi 10.5281/zenodo.15009519),
unzipped:

  qval_egenes/*_qval.tsv                    one table per cell type. One row per tested gene,
                                            with list-valued columns holding that gene's lead
                                            repeat (or tied leads) and a gene-level q-value.
                                            `CD4_TCM_permuted_qval.tsv` is the permutation
                                            null and is skipped.
  filtered_polymorphic_catalog.json         the 2,638,529 repeat loci the authors found
                                            polymorphic (ExpansionHunter catalog format).

How a lead repeat is placed on the genome: the qval tables give `chr`, `pos` and `motif`
only. `pos` is the 0-based start of the repeat, the same as the catalog's `ReferenceRegion`
start and as this project's BED files (checked against all 60,571 lead repeats of the
released tables: every one resolves). Each lead is matched to the catalog repeat that starts
there and has its motif, which gives its end. Multi-repeat catalog loci are split into
their component repeats first, because each component is its own tested variant.

For each STRchive locus the output records:
  polymorphic   a catalog repeat overlaps the locus,
  leads         lead repeats that overlap the locus (relation "overlaps") or lie within
                --adjacent-bp of it ("adjacent", a different repeat next door), each with
                the genes it is the lead for and the per-cell-type effect.

Only lead repeats are in the released tables, so a locus with no lead is "not a lead",
not "no effect". The effect is the fixed-effects meta-analysis coefficient per summed
diploid repeat length, on normalized expression; the q-value is gene-level.

The source is CC BY-NC-ND 4.0, so the output stays out of the repo (the API reads it from
~/sc_etr_annotations by default).

Usage:
    python3 build_sc_etr_disease_annotation.py \
        --qval-dir /path/to/qval_egenes \
        --catalog-json /path/to/filtered_polymorphic_catalog.json \
        --strchive project/uploads/STRchive-loci.json \
        --out ~/sc_etr_annotations/disease_panel.hg38.json
"""

from __future__ import annotations

import argparse
import ast
import csv
import datetime
import json
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Iterator

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from repeatscope.motifs import is_same_repeat_unit  # noqa: E402

csv.field_size_limit(sys.maxsize)

SCHEMA_VERSION = 1
MAX_Q = 0.05
ADJACENT_BP = 1000

PAPER = {
    "title": "Tandem repeat variation shapes immune cell type-specific gene expression",
    "authors": "Tanudisastro et al. (TenK10K)",
    "doi": "10.1101/2024.11.02.621562",
    "version": "v4, posted 2026-03-03",
    "license": "CC BY-NC-ND 4.0",
    "data_doi": "10.5281/zenodo.15009519",
    "cohort": "TenK10K Phase 1: 1,925 people (Tasmanian Ophthalmic Biobank and BioHEART), 28 immune cell types",
}

_QVAL_LIST_COLUMNS = ("chr", "pos", "motif", "coeff_meta_fixed", "se_meta_fixed", "pval_meta_fixed")
_MOTIF_RE = re.compile(r"\(([ACGTN]+)\)")
_CATALOG_KEYS = {
    "ReferenceRegion", "LocusStructure", "LocusId", "VariantId", "ManeGeneName", "ManeGeneId",
}


def _sig(x: float) -> float:
    """Four significant digits: the released values carry far more than anyone reads."""
    return float(f"{x:.4g}")


# ---------------------------------------------------------------- lead repeats
def read_leads(qval_dir: str | Path, max_q: float = MAX_Q) -> tuple[dict, int]:
    """Gene-level significant lead repeats from the per-cell-type tables.

    Returns ({(chrom, pos, motif): {gene_id: {cell_type: effect}}}, number of significant
    gene-cell-type rows). `pos` is kept as released (0-based start)."""
    files = sorted(Path(qval_dir).glob("*_qval.tsv"))
    if not files:
        raise SystemExit(f"no *_qval.tsv files in {qval_dir}")
    leads: dict[tuple, dict] = {}
    n_sig_rows = 0
    for path in files:
        cell_type = path.name[: -len("_qval.tsv")]
        if cell_type.endswith("_permuted"):
            continue
        with open(path, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                q = float(row["qval"])
                if q >= max_q:
                    continue
                n_sig_rows += 1
                cols = {c: ast.literal_eval(row[c]) for c in _QVAL_LIST_COLUMNS}
                for i in range(len(cols["chr"])):
                    key = (cols["chr"][i], int(cols["pos"][i]), cols["motif"][i])
                    leads.setdefault(key, {}).setdefault(row["gene_name"], {})[cell_type] = {
                        "coeff": _sig(cols["coeff_meta_fixed"][i]),
                        "se": _sig(cols["se_meta_fixed"][i]),
                        "pval": _sig(cols["pval_meta_fixed"][i]),
                        "qval": _sig(q),
                    }
    return leads, n_sig_rows


# ---------------------------------------------------------------- catalog
def iter_catalog_repeats(path: str | Path) -> Iterator[dict]:
    """One dict per repeat of the paper's catalog JSON.

    The file is pretty-printed (indent 2), so a line state machine streams it without
    holding 2.6M objects. A multi-repeat locus has list-valued ReferenceRegion/VariantId
    and yields one repeat per component; its LocusStructure lists the motifs in order."""
    cur: dict | None = None
    arr_key: str | None = None
    arr: list = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if arr_key is not None:
                if s in ("]", "],"):
                    if arr_key in _CATALOG_KEYS:
                        cur[arr_key] = arr
                    arr_key = None
                elif arr_key in _CATALOG_KEYS:
                    arr.append(s.rstrip(",").strip('"'))
                continue
            if s == "{":
                cur = {}
                continue
            if s in ("}", "},"):
                if cur is not None:
                    yield from _repeats_of(cur)
                cur = None
                continue
            if cur is None or s in ("[", "]"):
                continue
            key, _, value = s.partition(": ")
            key = key.strip('"')
            if value == "[":
                arr_key, arr = key, []
            elif key in _CATALOG_KEYS:
                cur[key] = value.rstrip(",").strip('"')


def _repeats_of(locus: dict) -> Iterator[dict]:
    regions = locus["ReferenceRegion"]
    regions = regions if isinstance(regions, list) else [regions]
    motifs = _MOTIF_RE.findall(locus.get("LocusStructure", ""))
    if len(motifs) == len(regions):
        per_repeat = motifs
    elif len(regions) == 1:
        per_repeat = ["|".join(motifs)]
    else:
        per_repeat = [""] * len(regions)
    gene_id = locus.get("ManeGeneId", "null")
    gene = locus.get("ManeGeneName", "null")
    for region, motif in zip(regions, per_repeat):
        chrom, _, span = region.partition(":")
        start, _, end = span.partition("-")
        yield {
            "chrom": chrom,
            "start": int(start),
            "end": int(end),
            "motif": motif,
            "gene_id": "" if gene_id == "null" else gene_id.split(".")[0],
            "gene": "" if gene == "null" else gene,
        }


# ---------------------------------------------------------------- STRchive panel
def load_panel(strchive_json: str | Path) -> list[dict]:
    with open(strchive_json, encoding="utf-8") as f:
        raw = json.load(f)
    panel = []
    for r in raw:
        motifs = set()
        for key in (
            "reference_motif_reference_orientation",
            "pathogenic_motif_reference_orientation",
            "benign_motif_reference_orientation",
        ):
            motifs.update(r.get(key) or [])
        panel.append(
            {
                "id": r["id"],
                "gene": r["gene"],
                "chrom": r["chrom"],
                "start": r["start_hg38"],
                "stop": r["stop_hg38"],
                "motifs": sorted(motifs),
            }
        )
    return panel


# ---------------------------------------------------------------- join
def build_annotations(
    leads: dict,
    repeats: Iterable[dict],
    panel: list[dict],
    n_sig_rows: int = 0,
    adjacent_bp: int = ADJACENT_BP,
    built: str | None = None,
) -> dict:
    lead_starts: dict[tuple, list] = defaultdict(list)
    for key in leads:
        lead_starts[(key[0], key[1])].append(key)
    panel_by_chrom: dict[str, list] = defaultdict(list)
    for p in panel:
        panel_by_chrom[p["chrom"]].append(p)

    resolved: dict[tuple, tuple[int, int]] = {}
    symbols: dict[str, str] = {}
    polymorphic: set[str] = set()
    n_repeats = 0
    for r in repeats:
        n_repeats += 1
        if r["gene_id"] and r["gene"]:
            symbols[r["gene_id"]] = r["gene"]
        for key in lead_starts.get((r["chrom"], r["start"]), ()):
            if key not in resolved and key[2] in r["motif"].split("|"):
                resolved[key] = (r["start"], r["end"])
        for p in panel_by_chrom.get(r["chrom"], ()):
            if r["start"] < p["stop"] and r["end"] > p["start"]:
                polymorphic.add(p["id"])

    resolved_by_chrom: dict[str, list] = defaultdict(list)
    for key, (start, end) in resolved.items():
        resolved_by_chrom[key[0]].append((key, start, end))

    loci: dict[str, dict] = {}
    for p in panel:
        found = []
        for key, start, end in resolved_by_chrom.get(p["chrom"], ()):
            if start < p["stop"] and end > p["start"]:
                relation, distance = "overlaps", 0
            else:
                distance = max(start - p["stop"], p["start"] - end)
                if distance > adjacent_bp:
                    continue
                relation = "adjacent"
            targets = []
            for gene_id, by_cell in leads[key].items():
                effects = [{"cell_type": ct, **eff} for ct, eff in sorted(by_cell.items())]
                targets.append({"gene_id": gene_id, "gene": symbols.get(gene_id), "effects": effects})
            targets.sort(key=lambda t: (min(e["qval"] for e in t["effects"]), t["gene_id"]))
            found.append(
                {
                    "chrom": key[0],
                    "start": start,
                    "stop": end,
                    "motif": key[2],
                    "relation": relation,
                    "distance_bp": distance,
                    "same_motif": any(is_same_repeat_unit(key[2], m) for m in p["motifs"]),
                    "targets": targets,
                }
            )
        found.sort(
            key=lambda l: (
                l["relation"] != "overlaps",
                not l["same_motif"],
                min(e["qval"] for t in l["targets"] for e in t["effects"]),
                l["start"],
            )
        )
        loci[p["id"]] = {"polymorphic": p["id"] in polymorphic, "leads": found}

    at_repeat = sum(1 for e in loci.values() if any(l["relation"] == "overlaps" for l in e["leads"]))
    nearby_only = sum(
        1
        for e in loci.values()
        if e["leads"] and not any(l["relation"] == "overlaps" for l in e["leads"])
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "source": {
            "paper": PAPER,
            "genome_build": "GRCh38",
            "lead_only": True,
            "significance": f"gene-level q < {MAX_Q} (as released)",
            "effect": "fixed-effects meta-analysis coefficient per summed diploid repeat length, on normalized expression",
            "adjacent_bp": adjacent_bp,
            "built": built or datetime.date.today().isoformat(),
            "counts": {
                "catalog_repeats": n_repeats,
                "significant_gene_cell_type_rows": n_sig_rows,
                "lead_repeats": len(leads),
                "lead_repeats_placed": len(resolved),
                "egenes": len({g for by_gene in leads.values() for g in by_gene}),
                "panel_loci": len(panel),
                "panel_loci_polymorphic": len(polymorphic),
                "panel_loci_with_lead_at_repeat": at_repeat,
                "panel_loci_with_lead_nearby_only": nearby_only,
            },
        },
        "loci": loci,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--qval-dir", required=True, help="unzipped qval_egenes directory")
    ap.add_argument("--catalog-json", required=True, help="unzipped filtered_polymorphic_catalog.json")
    ap.add_argument("--strchive", required=True, help="STRchive-loci.json")
    ap.add_argument("--out", required=True, help="annotation JSON to write (outside the repo)")
    ap.add_argument("--adjacent-bp", type=int, default=ADJACENT_BP)
    args = ap.parse_args()

    leads, n_sig_rows = read_leads(args.qval_dir)
    print(f"{len(leads):,} lead repeats from {n_sig_rows:,} significant gene-cell-type rows", file=sys.stderr)
    panel = load_panel(args.strchive)
    result = build_annotations(
        leads, iter_catalog_repeats(args.catalog_json), panel, n_sig_rows, adjacent_bp=args.adjacent_bp
    )
    counts = result["source"]["counts"]
    unplaced = counts["lead_repeats"] - counts["lead_repeats_placed"]
    if unplaced:
        print(f"warning: {unplaced:,} lead repeats matched no catalog repeat and were left out", file=sys.stderr)
    out = Path(args.out).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {out}", file=sys.stderr)
    for k, v in counts.items():
        print(f"  {k}: {v:,}", file=sys.stderr)


if __name__ == "__main__":
    main()
