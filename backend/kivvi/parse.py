"""Parsing and interpretation of Kivvi output (see Kivvi's docs/json.md, docs/vcf.md).

Pure functions over already-loaded JSON / VCF text -- no filesystem access, so
they are unit-testable with small synthetic inputs.
"""

from __future__ import annotations

import re
from typing import Optional

FLANKS = ("LeftFlank", "RightFlank")

# Kivvi's docs: FSHD1 is caused by contraction of D4Z4 to 1-10 repeat units on a
# qA haplotype with an intact PolyA signal. Kivvi defines no other cutoffs, so
# none are applied here (in particular no hypomethylation or KIV2 cutoff).
FSHD1_MAX_UNITS = 10
PATHOGENIC_DISTAL = "qAIntactPolyA"
# chr4 upstream group "Group2.1" is the 4A166 haplotype, "not pathogenic when
# contracted" (Kivvi docs/json.md).
NON_PATHOGENIC_UPSTREAM_GROUP = "Group2.1"
# Kivvi's D4Z4 preset "expects four alleles" (two chr4 + two chr10).
EXPECTED_D4Z4_ALLELES = 4
EXPECTED_KIV2_ALLELES = 2


def unit_ids(allele_name: str) -> list[str]:
    """Ordered repeat-unit IDs of an allele name ('LeftFlank-1-14-...-RightFlank')."""
    return [t for t in allele_name.split("-") if t and t not in FLANKS]


def parse_allele_size(raw: str) -> tuple[int, bool]:
    """'28' -> (28, False); '>14' -> (14, True) (a lower bound for partial alleles)."""
    text = str(raw).strip()
    lower_bound = text.startswith(">")
    return int(text.lstrip(">")), lower_bound


def split_chromosome(raw: str) -> tuple[str, Optional[str]]:
    """'chr4:Group1.1' -> ('chr4', 'Group1.1'); 'unknown' -> ('unknown', None)."""
    chrom, _, group = (raw or "unknown").partition(":")
    return chrom, (group or None)


def classify_d4z4_allele(
    size: int, lower_bound: bool, chromosome: str, group: Optional[str], distal: str
) -> dict:
    """FSHD1 relevance of one D4Z4 allele, with the reason spelled out."""
    if lower_bound:
        if size > FSHD1_MAX_UNITS:
            return _call("normal", "Not contracted", f"At least {size} units (partially assembled) — above the 1–{FSHD1_MAX_UNITS} contraction range.")
        return _call("indeterminate", "Indeterminate", f"Only a lower bound of {size} units was assembled, so contraction cannot be excluded.")

    if size > FSHD1_MAX_UNITS:
        return _call("normal", "Not contracted", f"{size} units — above the 1–{FSHD1_MAX_UNITS} contraction range.")

    # From here the allele is contracted (<= FSHD1_MAX_UNITS); whether that matters
    # depends on the haplotype it sits on.
    if distal != PATHOGENIC_DISTAL:
        return _call("normal", "Contracted, non-pathogenic", f"{size} units, but the distal haplotype is {distal}; only {PATHOGENIC_DISTAL} alleles can express DUX4.")
    if chromosome == "chr4" and group == NON_PATHOGENIC_UPSTREAM_GROUP:
        return _call("normal", "Contracted, non-pathogenic", f"{size} units on the chr4 {group} (4A166) haplotype, which is not pathogenic when contracted.")
    if chromosome == "chr4":
        return _call("pathogenic", "FSHD1-range contraction", f"{size} units on chr4 with a qA haplotype and intact PolyA signal.")
    return _call("review", "Review", f"{size} units with an intact-PolyA qA haplotype, but the allele is assigned to {chromosome}; confirm its chromosome of origin before interpreting.")


def _call(status: str, label: str, reason: str) -> dict:
    return {"status": status, "label": label, "reason": reason}


def _methylation_by_allele(data: dict) -> dict[str, list[float]]:
    """allele name -> per-unit methylation, from methylation.all_allele_ends
    (which also covers partial alleles)."""
    out: dict[str, list[float]] = {}
    per_unit = (data.get("methylation") or {}).get("all_allele_ends", {}).get("median_methylation_per_unit", {})
    for name, csv in per_unit.items():
        out[name] = [float(v) for v in str(csv).split(",") if v != ""]
    return out


def shape_d4z4(data: dict) -> dict:
    """D4Z4 alleles, interpretation and methylation from a Kivvi D4Z4 JSON."""
    per_unit_meth = _methylation_by_allele(data)
    alleles = []
    for info in data.get("allele_info", []):
        size, lower_bound = parse_allele_size(info["allele_size"])
        chromosome, group = split_chromosome(info.get("chromosome", "unknown"))
        distal = info.get("distal_haplotype", "unknown")
        name = info["allele_name"]
        alleles.append(
            {
                "name": name,
                "units": unit_ids(name),
                "size": size,
                "size_is_lower_bound": lower_bound,
                "size_label": f"≥{size}" if lower_bound else str(size),
                "chromosome": chromosome,
                "upstream_group": group,
                "distal_haplotype": distal,
                "allele_type": info.get("allele_type", ""),
                "methylation": info.get("methylation"),
                "ending_in_qal": bool(info.get("ending_in_qal", False)),
                "per_unit_methylation": per_unit_meth.get(name, []),
                "interpretation": classify_d4z4_allele(size, lower_bound, chromosome, group, distal),
            }
        )

    additional = data.get("additional") or {}
    return {
        "alleles": alleles,
        "summary": _d4z4_summary(alleles),
        "overall_methylation": additional.get("methylation_all_sites"),
        "in_cis_arrays": additional.get("d4z4_arrays_in_cis") or [],
        "qc": additional.get("qc_metrics") or {},
    }


def _d4z4_summary(alleles: list[dict]) -> dict:
    by_status = [a["interpretation"]["status"] for a in alleles]
    complete = sum(1 for a in alleles if not a["size_is_lower_bound"])
    if "pathogenic" in by_status:
        status, label = "pathogenic", "FSHD1-range contraction detected"
    elif "review" in by_status:
        status, label = "review", "Review: contracted qA allele with unconfirmed chromosome"
    elif "indeterminate" in by_status:
        status, label = "indeterminate", "No contraction detected, but not every allele is fully resolved"
    elif alleles:
        status, label = "normal", "No FSHD1-range contraction"
    else:
        status, label = "indeterminate", "No D4Z4 allele was assembled"
    notes = []
    if len(alleles) < EXPECTED_D4Z4_ALLELES:
        notes.append(f"Only {len(alleles)} of {EXPECTED_D4Z4_ALLELES} expected D4Z4 alleles were resolved.")
    return {
        "status": status,
        "label": label,
        "complete_alleles": complete,
        "partial_alleles": len(alleles) - complete,
        "notes": notes,
    }


def shape_kiv2(data: dict) -> dict:
    """KIV2 copy number per assembled allele from a Kivvi KIV2 JSON.

    Kivvi emits no per-allele metadata for KIV2 (allele_info is empty), so alleles
    are read straight off complete_alleles; supporting_reads is keyed by the same
    allele names. No clinical cutoff is applied -- Kivvi only notes that a short
    KIV2 allele is associated with higher cardiovascular risk.
    """
    supporting = data.get("supporting_reads") or {}
    alleles = []
    for index, name in enumerate(data.get("complete_alleles", []), start=1):
        units = unit_ids(name)
        alleles.append(
            {
                "index": index,
                "name": name,
                "units": units,
                "copy_number": len(units),
                "supporting_reads": len(supporting.get(name, [])),
            }
        )
    partial = [{"name": n, "units": unit_ids(n), "lower_bound": len(unit_ids(n))} for n in data.get("partial_alleles", [])]

    notes = []
    if len(alleles) < EXPECTED_KIV2_ALLELES:
        notes.append(
            f"Only {len(alleles)} of {EXPECTED_KIV2_ALLELES} expected KIV2 alleles were assembled; "
            "a second allele of the same length may have collapsed into one, or assembly failed."
        )
    shortest = min((a["copy_number"] for a in alleles), default=None)
    return {
        "alleles": alleles,
        "partial_alleles": partial,
        "shortest_copy_number": shortest,
        "allele_cn": data.get("allele_cn", ""),
        "notes": notes,
    }


_RU_RE = re.compile(r"^(?P<unit>[^:]+):(?P<where>[^:]+):(?P<depth>\d+):(?P<reads>\d+)$")


def parse_vcf(text: str) -> dict:
    """Parse a Kivvi VCF: header facts plus one record per variant.

    Kivvi repurposes the sample columns so each is one assembled allele (GT 0/1),
    and adds INFO RU=<unit>:<allele.index>:<depth>:<reads>[,...]. Records whose RU
    position is 'Unknown' sit on repeat units that were not assembled into a
    complete allele, so no allele carries them.
    """
    version = None
    command = None
    allele_columns: list[str] = []
    records = []
    for line in text.splitlines():
        if line.startswith("##kivviVersion="):
            version = line.split("=", 1)[1]
        elif line.startswith("##kivviCommand="):
            command = line.split("=", 1)[1]
        elif line.startswith("#CHROM"):
            allele_columns = line.split("\t")[9:]
        elif line and not line.startswith("#"):
            f = line.split("\t")
            if len(f) < 10:
                continue
            units = []
            for ru in _info_value(f[7], "RU").split(","):
                m = _RU_RE.match(ru)
                if m:
                    units.append(
                        {
                            "unit": m["unit"],
                            "where": m["where"],
                            "depth": int(m["depth"]),
                            "reads": int(m["reads"]),
                        }
                    )
            genotypes = f[9:]
            records.append(
                {
                    "chrom": f[0],
                    "pos": int(f[1]),
                    "ref": f[3],
                    "alt": f[4],
                    "filter": f[6],
                    "genotypes": genotypes,
                    "on_assembled_allele": any(g not in ("0", ".") for g in genotypes),
                    "units": units,
                    "n_units": len(units),
                    # Weakest read support across the units carrying the variant.
                    "min_vaf": min((u["reads"] / u["depth"] for u in units if u["depth"]), default=None),
                }
            )
    return {"version": version, "command": command, "allele_columns": allele_columns, "records": records}


def _info_value(info: str, key: str) -> str:
    for part in info.split(";"):
        if part.startswith(key + "="):
            return part.split("=", 1)[1]
    return ""
