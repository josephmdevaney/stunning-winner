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
# Chen et al. (bioRxiv 2026, the Kivvi preprint) call 11-20 units "the range commonly indicated for FSHD2" and note
# that Kivvi's accuracy there was not measured. It only labels an allele's size for the FSHD2 evidence view; no call
# depends on it.
FSHD2_RANGE = (11, 20)
# Kivvi's D4Z4 preset "expects four alleles" (two chr4 + two chr10).
EXPECTED_D4Z4_ALLELES = 4
EXPECTED_KIV2_ALLELES = 2


def unit_ids(allele_name: str) -> list[str]:
    """Ordered repeat-unit IDs of an allele name ('LeftFlank-1-14-...-RightFlank')."""
    return [t for t in allele_name.split("-") if t and t not in FLANKS]


_ALLELE_SIZE = re.compile(r"(>=?)?\s*(\d+)")


def parse_allele_size(raw: str) -> tuple[int, bool]:
    """'28' -> (28, False); '>14' and '>=14' -> (14, True).

    Kivvi writes a lower bound for partial and merged alleles as '>=N' ("at least N") or '>N' ("more than
    N"); the docs show only '>N'. Both are kept as a lower bound of N, which for '>N' is one unit short, so
    an allele's size is never overstated. Anything else is refused rather than guessed at.
    """
    match = _ALLELE_SIZE.fullmatch(str(raw).strip())
    if match is None:
        raise ValueError(f"Unrecognised Kivvi allele_size {raw!r}; expected N, >N or >=N")
    return int(match.group(2)), match.group(1) is not None


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
        # Not yet known to be above the range, so the haplotype is read first: one that cannot cause FSHD1 settles
        # the allele whatever its final size (the same two rules as for an exact size, below).
        if distal != PATHOGENIC_DISTAL:
            why = f"the distal haplotype is {distal} and only {PATHOGENIC_DISTAL} alleles can express DUX4"
        elif chromosome == "chr4" and group == NON_PATHOGENIC_UPSTREAM_GROUP:
            why = f"it is on the chr4 {group} (4A166) haplotype, which is not pathogenic when contracted"
        else:
            return _call("indeterminate", "Indeterminate", f"Only a lower bound of {size} units was assembled, so contraction cannot be excluded.")
        return _call("normal", "Non-pathogenic haplotype", f"At least {size} units (partially assembled), but {why}, so it is not pathogenic whatever its final size.")

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


def is_permissive(chromosome: str, group: Optional[str], distal: str) -> bool:
    """A chr4 qAIntactPolyA allele outside the 4A166 haplotype: the kind that can express DUX4. The same haplotype
    test as the FSHD1 rule above, without the size condition."""
    return chromosome == "chr4" and distal == PATHOGENIC_DISTAL and group != NON_PATHOGENIC_UPSTREAM_GROUP


def size_range(size: int, lower_bound: bool) -> str:
    """Where an allele's size sits: 'contracted' (1-10 units, FSHD1 range), 'fshd2_range' (11-20), 'above' (more
    than 20) or 'unresolved' (a partially assembled allele whose lower bound does not settle it)."""
    low, high = FSHD2_RANGE
    if lower_bound:
        return "above" if size > high else "unresolved"
    if size <= FSHD1_MAX_UNITS:
        return "contracted"
    return "fshd2_range" if size <= high else "above"


def shape_fshd2(alleles: list[dict], sample_methylation: Optional[float]) -> dict:
    """The facts that bear on FSHD2, laid side by side. Nothing is concluded from them.

    Kivvi defines no FSHD2 rule and no methylation cut-off, and its preprint validated it on FSHD1 and control
    samples only, so this gathers what Kivvi reports (each allele's haplotype, size and methylation, and the sample's
    overall methylation) and leaves the reading to the reader. Variants in SMCHD1, DNMT3B or LRIF1, the other part of
    an FSHD2 finding, are not something Kivvi reports and are not assessed.
    """
    rows = []
    for index, a in enumerate(alleles, start=1):
        rows.append(
            {
                "index": index,
                "size": a["size"],
                "size_label": a["size_label"],
                "size_is_lower_bound": a["size_is_lower_bound"],
                "size_range": size_range(a["size"], a["size_is_lower_bound"]),
                "chromosome": a["chromosome"],
                "upstream_group": a["upstream_group"],
                "distal_haplotype": a["distal_haplotype"],
                "permissive": is_permissive(a["chromosome"], a["upstream_group"], a["distal_haplotype"]),
                "methylation": a["methylation"],
            }
        )
    permissive = [r for r in rows if r["permissive"]]
    levels = [r["methylation"] for r in rows if r["methylation"] is not None]
    return {
        "alleles": rows,
        "sample_methylation": sample_methylation,
        "summary": {
            "permissive_alleles": len(permissive),
            "permissive_contracted": sum(r["size_range"] == "contracted" for r in permissive),
            "permissive_fshd2_range": sum(r["size_range"] == "fshd2_range" for r in permissive),
            "alleles_with_methylation": len(levels),
            "lowest_methylation": min(levels, default=None),
            "highest_methylation": max(levels, default=None),
        },
    }


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
        "fshd2": shape_fshd2(alleles, additional.get("methylation_all_sites")),
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
            # With no complete allele Kivvi writes no FORMAT or sample columns at all (8 fields); those records
            # are still variants, on units no allele carries.
            if len(f) < 8:
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
