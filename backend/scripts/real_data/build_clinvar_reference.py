"""Build a project/uploads/clinvar_<gene>.hg38.tsv from ClinVar's own current
variant classifications for one gene's region.

(Copied from the carrier tool, puretarget-carrier-screen. Here the tables are used only for the
stale-havanno-variant check in paraphase.py's `_classify_havanno_gene`; the pseudogene-pair caveat
below concerns that tool's separate small-variant classification, which this app does not have.)

ClinVar (NCBI, https://www.ncbi.nlm.nih.gov/clinvar/) is the real, current,
publicly citable clinical variant classification database. This script pulls
just the given gene's region (the exact span PTCP's own Paraphase config
realigns and calls small variants against -- see each gene's
`realign_region` in Paraphase's own GRCh38 config,
https://raw.githubusercontent.com/PacificBiosciences/paraphase/main/data/38/config.yaml)
from ClinVar's own published VCF, using an indexed remote region query so
this never has to download the full multi-hundred-MB file.

IMPORTANT caveat this project's own classification logic must respect (see
paraphase.py's `_classify_small_variants`) for any gene with a same-size
pseudogene pair (CYP21A2, GBA -- see `_HAS_PSEUDOGENE_PAIR`): a handful of
well-known "gene-pseudogene differentiating" positions are BOTH the
pseudogene's own completely normal sequence AND (when gene-converted onto
the real functional gene copy) the most common real pathogenic mutations --
the exact same allele is benign on one physical copy and pathogenic on the
other. ClinVar's own CLNSIG for those specific positions is "Pathogenic"
(that's correct when the position really is evaluated on the gene copy)
with no way to tell from the VCF alone which copy a given call sits on.
This reference file is therefore only sound to consult for calls that this
project's own per-haplotype VCF parsing has already determined are NOT part
of that recurring differentiating-site block (isolated, cohort-rare
calls) -- consulting it for a cluster-pattern call would silently reproduce
the exact false-positive bug this project already found and fixed for
havanno.

Also joins in ClinVar's own preferred HGVS transcript notation (its
tab-delimited `variant_summary.txt.gz`, "Name" column, e.g.
"NM_000157.4(GBA1):c.1448T>G (p.Leu483Arg)") for each entry, reformatted to
this project's own "<gene> <c.> <p.>" display style -- the same style
`_format_curated_variant_name` in paraphase.py produces from PacBio's own
curated variant names. `_classify_small_variants` prefers a PacBio curated
name when one exists for a call's exact position (that project's own
purpose-built curation for this specific carrier panel); this HGVS field is
what it falls back to for the (far more common) case where ClinVar has
classified a call but PacBio's own small, hand-curated list hasn't -- so a
displayed variant is (almost) never just a bare genomic coordinate.

Usage:
    python -m scripts.real_data.build_clinvar_reference \\
        --gene cyp21 --region 6:32038085-32045442 \\
        --out project/uploads/clinvar_cyp21.hg38.tsv

Requires `bcftools` on PATH (used for the remote indexed region query).
Re-run this periodically to pick up new ClinVar submissions/reclassifications.

Offline: pass `--variant-summary <variant_summary_YYYY-MM.txt.gz>` (a ClinVar release's tab-delimited
variant summary, e.g. from a folder of monthly releases) and everything is built from that one file,
with no bcftools and no download. Its human-readable labels ("Likely benign", "criteria provided,
single submitter") become the VCF's ("Likely_benign", "criteria_provided,_single_submitter"), so the
tables come out in the same format either way.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import subprocess
import urllib.request
from collections import Counter
from pathlib import Path

CLINVAR_VCF_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
CLINVAR_VARIANT_SUMMARY_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"

# ClinVar's "Name" column is consistently "<transcript>(<gene>):<c. HGVS>"
# optionally followed by " (<p. HGVS>)" -- e.g. "NM_000157.4(GBA1):c.1448T>G
# (p.Leu483Arg)" or, for a variant with no protein consequence to report
# (intronic, UTR), just "NM_000157.4(GBA1):c.203-10T>C" with no trailing
# parenthetical. Validated against real GBA1/HBB entries (see this module's
# tests) -- the trailing "$" anchor stops the c. group from swallowing the
# optional " (p...)" suffix when one is present.
_CLINVAR_NAME_RE = re.compile(r"\(([^)]+)\):(c\.\S+?)(?:\s+\((p\.[^)]+)\))?$")


def _format_clinvar_name(name: str) -> str:
    match = _CLINVAR_NAME_RE.search(name)
    if not match:
        return ""
    gene, cdna, protein = match.groups()
    return f"{gene} {cdna}" + (f" {protein}" if protein else "")


def _parse_info(info: str) -> dict[str, str]:
    fields: dict[str, str] = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            fields[k] = v
    return fields


def fetch_region_vcf_lines(region: str) -> list[str]:
    result = subprocess.run(
        ["bcftools", "view", "-r", region, CLINVAR_VCF_URL],
        capture_output=True, text=True, check=True,
    )
    return [line for line in result.stdout.splitlines() if line and not line.startswith("#")]


def fetch_hgvs_by_position(chrom: str, positions: set[tuple[str, int, str, str]]) -> dict[tuple[str, int, str, str], str]:
    """Stream ClinVar's full `variant_summary.txt.gz` (there's no
    region-indexed version of this file, unlike the VCF) and pull just the
    GRCh38, matching-chromosome rows whose exact (chrom, pos, ref, alt) --
    ClinVar's own VCF-normalized columns, `PositionVCF`/`*AlleleVCF` -- is
    in `positions` (the calls this gene's region query already found).
    Streamed rather than fully buffered since the file is several hundred
    MB uncompressed."""
    hgvs_by_position: dict[tuple[str, int, str, str], str] = {}
    with urllib.request.urlopen(CLINVAR_VARIANT_SUMMARY_URL) as response:
        with gzip.open(response, mode="rt", encoding="utf-8", errors="replace") as f:
            header = next(f).rstrip("\n").split("\t")
            name_i, chrom_i, assembly_i = header.index("Name"), header.index("Chromosome"), header.index("Assembly")
            pos_i, ref_i, alt_i = header.index("PositionVCF"), header.index("ReferenceAlleleVCF"), header.index("AlternateAlleleVCF")
            for line in f:
                fields = line.rstrip("\n").split("\t")
                if fields[chrom_i] != chrom or fields[assembly_i] != "GRCh38":
                    continue
                if fields[pos_i] == "-":
                    continue
                key = (f"chr{fields[chrom_i]}", int(fields[pos_i]), fields[ref_i], fields[alt_i])
                if key in positions and key not in hgvs_by_position:
                    formatted = _format_clinvar_name(fields[name_i])
                    if formatted:
                        hgvs_by_position[key] = formatted
    return hgvs_by_position


def build(region: str) -> list[dict]:
    rows: list[dict] = []
    for line in fetch_region_vcf_lines(region):
        fields = line.split("\t")
        chrom, pos, _id, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
        info = _parse_info(fields[7])
        clnsig = info.get("CLNSIG", "")
        if not clnsig:
            continue
        rows.append({
            "chrom": f"chr{chrom}",
            "pos": int(pos),
            "ref": ref,
            "alt": alt,
            "clnsig": clnsig,
            "clnrevstat": info.get("CLNREVSTAT", ""),
            "gene": info.get("GENEINFO", "").split(":")[0],
            "condition": info.get("CLNDN", ""),
            "rsid": info.get("RS", ""),
        })
    rows.sort(key=lambda r: r["pos"])

    chrom = region.split(":", 1)[0].removeprefix("chr")
    positions = {(r["chrom"], r["pos"], r["ref"], r["alt"]) for r in rows}
    hgvs_by_position = fetch_hgvs_by_position(chrom, positions)
    for r in rows:
        r["hgvs"] = hgvs_by_position.get((r["chrom"], r["pos"], r["ref"], r["alt"]), "")
    return rows


def _vcf_style(label: str) -> str:
    """variant_summary's label in the VCF's spelling: multiple values joined by "|", spaces as "_"."""
    return label.replace("; ", "|").replace(" ", "_")


def build_from_variant_summary(region: str, variant_summary_path: Path) -> list[dict]:
    """Rows for `region` (GRCh38) from a local ClinVar `variant_summary` file, one streamed pass.

    A variant is in the region if its VCF-normalized span [PositionVCF, PositionVCF + len(ref) - 1]
    overlaps it, as `bcftools view -r` selects. Rows ClinVar gives no classification ("-") are left
    out, as in the VCF build."""
    chrom, span = region.split(":", 1)
    chrom = chrom.removeprefix("chr")
    start, stop = (int(x) for x in span.replace(",", "").split("-"))
    rows: dict[tuple, dict] = {}
    with gzip.open(variant_summary_path, mode="rt", encoding="utf-8", errors="replace") as f:
        header = next(f).rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if fields[col["Assembly"]] != "GRCh38" or fields[col["Chromosome"]] != chrom:
                continue
            position = fields[col["PositionVCF"]]
            if position == "-":
                continue
            pos, ref, alt = int(position), fields[col["ReferenceAlleleVCF"]], fields[col["AlternateAlleleVCF"]]
            if pos > stop or pos + len(ref) - 1 < start:
                continue
            clnsig = fields[col["ClinicalSignificance"]]
            if not clnsig or clnsig == "-":
                continue
            key = (f"chr{chrom}", pos, ref, alt)
            if key in rows:
                continue
            rsid = fields[col["RS# (dbSNP)"]]
            rows[key] = {
                "chrom": key[0], "pos": pos, "ref": ref, "alt": alt,
                "clnsig": _vcf_style(clnsig),
                "clnrevstat": _vcf_style(fields[col["ReviewStatus"]]),
                "gene": fields[col["GeneSymbol"]].split(";")[0],
                "condition": _vcf_style(fields[col["PhenotypeList"]]),
                "rsid": "" if rsid in ("-1", "-", "") else rsid,
                "hgvs": _format_clinvar_name(fields[col["Name"]]),
            }
    return sorted(rows.values(), key=lambda r: r["pos"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gene", required=True,
        help="Paraphase gene key (e.g. cyp21, gba, hba, hbb, smn1, rpgr) -- used only to name the default --out file.",
    )
    parser.add_argument(
        "--region", required=True,
        help="GRCh38 region to query, e.g. 6:32038085-32045442 -- use the gene's own realign_region from "
        "Paraphase's config so the reference never falls short of what the per-haplotype VCF calls against.",
    )
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument(
        "--variant-summary", type=Path, default=None,
        help="Build from this local ClinVar variant_summary_*.txt.gz instead of querying NCBI.",
    )
    args = parser.parse_args()

    out = args.out or (
        Path(__file__).resolve().parents[3] / "project" / "uploads" / f"clinvar_{args.gene}.hg38.tsv"
    )
    rows = build_from_variant_summary(args.region, args.variant_summary) if args.variant_summary else build(args.region)
    out.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["chrom", "pos", "ref", "alt", "clnsig", "clnrevstat", "gene", "condition", "rsid", "hgvs"]
    with open(out, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    counts = Counter(r["clnsig"] for r in rows)
    n_with_hgvs = sum(1 for r in rows if r["hgvs"])
    print(f"wrote {len(rows)} ClinVar entries for {args.region} to {out} ({n_with_hgvs} with HGVS c./p. notation)")
    for sig, n in counts.most_common(10):
        print(f"  {sig}: {n}")


if __name__ == "__main__":
    main()
