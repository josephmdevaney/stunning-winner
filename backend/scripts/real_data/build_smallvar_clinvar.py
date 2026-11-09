#!/usr/bin/env python3
"""Build the reduced ClinVar table the small-variant suite matches against.

Input is NCBI's ClinVar ``variant_summary.txt.gz`` (https://ftp.ncbi.nlm.nih.gov/pub/clinvar/
tab_delimited/variant_summary.txt.gz). The full file is about 450 MB compressed and several
million rows, too many to hold as a lookup, so this keeps the rows the review screens act on:

  * GRCh38 rows with VCF-style alleles (``PositionVCF``, ``ReferenceAlleleVCF``,
    ``AlternateAlleleVCF``) that are plain sequence and small (indel size < 50 bp), and
  * a germline classification that is pathogenic, likely pathogenic, or "Conflicting
    classifications of pathogenicity".

Alleles are trimmed to the parsimonious form (normalize.variant_key) so they compare equal
to a sample's trimmed alleles. ClinVar is public data from NCBI; the output is derived, large,
and release-specific, so it lives outside the repo (default ~/smallvar_data).

Usage:
    python3 build_smallvar_clinvar.py \
        --variant-summary variant_summary_2026-10.txt.gz \
        --out ~/smallvar_data/clinvar_plp.grch38.tsv.gz
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from smallvar.clinvar import COLUMNS, classification_group, review_stars  # noqa: E402
from smallvar.normalize import variant_key  # noqa: E402

csv.field_size_limit(sys.maxsize)
_ACGT = re.compile(r"^[ACGT]+$")
MAX_INDEL = 50
MAX_ALLELE = 200


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variant-summary", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--release", help="release label, e.g. 2026-10 (default: from the input file name)")
    args = ap.parse_args()

    src = Path(args.variant_summary)
    release = args.release or (re.search(r"(\d{4}-\d{2})", src.name) or [None, "unknown"])[1]
    out = Path(args.out).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)

    n_in = n_grch38 = n_kept = 0
    groups = {"plp": 0, "conflicting": 0}
    seen = set()
    opener = gzip.open if src.suffix == ".gz" else open
    with opener(src, "rt", encoding="utf-8") as fin, gzip.open(out, "wt", encoding="utf-8") as fout:
        fout.write(f"# release: {release}\n# source: ClinVar variant_summary ({src.name}); GRCh38 small variants, P/LP and conflicting\n")
        fout.write("\t".join(COLUMNS) + "\n")
        reader = csv.reader(fin, delimiter="\t")
        header = [h.lstrip("#") for h in next(reader)]
        ix = {name: i for i, name in enumerate(header)}
        for row in reader:
            n_in += 1
            if row[ix["Assembly"]] != "GRCh38":
                continue
            n_grch38 += 1
            ref, alt, pos = row[ix["ReferenceAlleleVCF"]], row[ix["AlternateAlleleVCF"]], row[ix["PositionVCF"]]
            if not (_ACGT.match(ref) and _ACGT.match(alt) and pos.isdigit()):
                continue
            if len(ref) > MAX_ALLELE or len(alt) > MAX_ALLELE or abs(len(ref) - len(alt)) >= MAX_INDEL:
                continue
            classification = row[ix["ClinicalSignificance"]]
            group = classification_group(classification)
            if group == "other":
                continue
            key = variant_key(row[ix["Chromosome"]], int(pos), ref, alt)
            if key in seen:  # variant_summary can repeat an allele across records
                continue
            seen.add(key)
            conditions = "; ".join(c for c in row[ix["PhenotypeList"]].split("|") if c)[:300]
            fout.write(
                "\t".join(
                    [key[0], str(key[1]), key[2], key[3], row[ix["GeneSymbol"]], classification, group,
                     str(review_stars(row[ix["ReviewStatus"]])), row[ix["ReviewStatus"]], row[ix["VariationID"]],
                     row[ix["Name"]], conditions, row[ix["LastEvaluated"]], row[ix["NumberSubmitters"]] or "0"]
                )
                + "\n"
            )
            n_kept += 1
            groups[group] += 1
    print(f"ClinVar {release}: {n_in:,} rows, {n_grch38:,} GRCh38, kept {n_kept:,} "
          f"({groups['plp']:,} P/LP, {groups['conflicting']:,} conflicting) -> {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
