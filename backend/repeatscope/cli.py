"""Command-line entry point: annotate a TRGT VCF against STRchive / a genome-wide repeat catalog.

Examples:
    repeatscope annotate --vcf sample.trgt.vcf.gz --strchive-json STRchive-loci.json --out panel.json
    repeatscope catalog --vcf sample.trgt.vcf.gz --bed repeat_catalog.hg38.bed --out catalog.json
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import sys
from pathlib import Path

from .annotate import annotate_strchive_panel, iter_annotated_catalog, simulate_pileup_reads
from .trgt import sample_names

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_STRCHIVE_JSON = REPO_ROOT / "project" / "uploads" / "STRchive-loci.json"
DEFAULT_CATALOG_BED = REPO_ROOT / "project" / "uploads" / "repeat_catalog.hg38.bed"


def _cmd_annotate(args: argparse.Namespace) -> int:
    sample = args.sample
    if sample is None:
        names = sample_names(args.vcf)
        if len(names) != 1:
            print(
                f"error: VCF has {len(names)} samples ({names}); pass --sample to pick one",
                file=sys.stderr,
            )
            return 1
        sample = names[0]

    panel = annotate_strchive_panel(args.vcf, args.strchive_json, sample=sample)
    loci = []
    for locus in panel:
        d = dataclasses.asdict(locus)
        d["region"] = locus.region
        if args.include_pileup:
            d["pileup_reads"] = simulate_pileup_reads(locus)
        loci.append(d)

    genotyped = sum(1 for l in panel if l.genotyped)
    pathogenic = sum(1 for l in panel if l.status == "Full expansion")
    out = {
        "sample": sample,
        "source_vcf": str(args.vcf),
        "loci_analyzed": len(panel),
        "loci_genotyped": genotyped,
        "pathogenic_findings": pathogenic,
        "loci": loci,
    }

    text = json.dumps(out, indent=2 if args.pretty else None)
    if args.out:
        Path(args.out).write_text(text)
        print(f"wrote {len(panel)} annotated loci ({genotyped} genotyped, {pathogenic} pathogenic) to {args.out}")
    else:
        print(text)
    return 0


def _cmd_catalog(args: argparse.Namespace) -> int:
    sample = args.sample
    if sample is None:
        names = sample_names(args.vcf)
        if len(names) != 1:
            print(
                f"error: VCF has {len(names)} samples ({names}); pass --sample to pick one",
                file=sys.stderr,
            )
            return 1
        sample = names[0]

    regions = list(iter_annotated_catalog(args.bed, args.vcf, sample=sample))
    genotyped = sum(1 for r in regions if r.genotyped)
    out = {
        "sample": sample,
        "source_vcf": str(args.vcf),
        "total_regions": len(regions),
        "genotyped_regions": genotyped,
        "regions": [dataclasses.asdict(r) for r in regions],
    }
    text = json.dumps(out, indent=2 if args.pretty else None)
    if args.out:
        Path(args.out).write_text(text)
        print(f"wrote {len(regions)} regions ({genotyped} genotyped) to {args.out}")
    else:
        print(text)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="repeatscope")
    sub = parser.add_subparsers(dest="command", required=True)

    p_annotate = sub.add_parser("annotate", help="annotate STRchive disease loci against a TRGT VCF")
    p_annotate.add_argument("--vcf", required=True, help="TRGT VCF (.vcf or .vcf.gz)")
    p_annotate.add_argument("--sample", default=None, help="sample name (defaults to the VCF's only sample)")
    p_annotate.add_argument("--strchive-json", default=str(DEFAULT_STRCHIVE_JSON), help="STRchive loci JSON")
    p_annotate.add_argument("--out", default=None, help="output JSON path (defaults to stdout)")
    p_annotate.add_argument("--pretty", action="store_true", help="pretty-print JSON output")
    p_annotate.add_argument(
        "--include-pileup", action="store_true",
        help="include per-locus estimated pileup reads (derived from SD/ALLR, not raw per-read data)",
    )
    p_annotate.set_defaults(func=_cmd_annotate)

    p_catalog = sub.add_parser("catalog", help="annotate the genome-wide repeat catalog against a TRGT VCF")
    p_catalog.add_argument("--vcf", required=True, help="TRGT VCF (.vcf or .vcf.gz)")
    p_catalog.add_argument("--sample", default=None, help="sample name (defaults to the VCF's only sample)")
    p_catalog.add_argument("--bed", default=str(DEFAULT_CATALOG_BED), help="genome-wide repeat catalog BED")
    p_catalog.add_argument("--out", default=None, help="output JSON path (defaults to stdout)")
    p_catalog.add_argument("--pretty", action="store_true", help="pretty-print JSON output")
    p_catalog.set_defaults(func=_cmd_catalog)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
