"""A sample's Paraphase JSON as a browsable table of regions -- no network and no
wait, unlike the annotation run (runner.py), which spends minutes to hours on
Ensembl/ClinVar lookups.

Everything here is read straight off fields Paraphase itself wrote; nothing is
re-derived or interpreted (no pathogenicity, no carrier calls -- that belongs to
the annotation run and to the Carrier panel). A "flag" only restates something the
JSON says (a fusion haplotype was called, a region had no usable data, ...).
Works on both generations of Paraphase output by normalizing to the v4 layout first
(see calls.py).
"""

from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path
from typing import Any, Optional

from . import calls, runner, samples

_MAX_VARIANTS_SHOWN = 25
_HAP_NAME = re.compile(r"^[^_]+_(?P<label>.+?)hap\d+$")

# Documented in Paraphase's docs/vcf.md: haplotypes are labelled gene1 / gene2
# (gene2 = the paralog or pseudogene) where Paraphase can tell them apart, plus
# "fusion" and "homology" haplotypes.
_KIND_BY_LABEL = {
    "gene1": "gene",
    "gene2": "paralog / pseudogene",
    "fusion": "fusion",
    "homology": "homology",
}


@lru_cache(maxsize=4)
def _load(path: str, mtime_ns: int) -> tuple[Optional[str], dict[str, dict]]:
    """(schema of the file as written, v4-normalized regions), cached per file
    version (the key includes the modification time, so an updated file is
    re-read)."""
    raw = calls.load_calls(Path(path))
    return calls.schema(raw), calls.to_v4(raw)


def _loaded(sample: str) -> tuple[Optional[str], dict[str, dict]]:
    path = samples.json_path(sample)
    return _load(str(path), path.stat().st_mtime_ns)


def _calls_for(sample: str) -> dict[str, dict]:
    return _loaded(sample)[1]


def _median_depth(region: dict) -> Optional[float]:
    return (region.get("region_depth") or {}).get("median")


def _specific(region: dict) -> dict:
    return region.get("region_specific_info") or {}


def _fusion_names(region: dict) -> list[str]:
    fusions = _specific(region).get("fusions_called")
    return sorted(fusions) if isinstance(fusions, dict) else []


def _sv_text(region: dict) -> Optional[str]:
    sv = _specific(region).get("sv_called")
    if not sv:
        return None
    if isinstance(sv, dict):
        return ", ".join(f"{k}: {v}" for k, v in sv.items())
    return str(sv)


def _cn(value: Any) -> str:
    return "?" if value is None else str(value)


def _headline(name: str, region: dict) -> Optional[str]:
    """One line of the region's own headline numbers, for the regions that have
    them. Strings only restate JSON fields."""
    spec = _specific(region)
    if "smn1_cn" in spec:
        text = f"SMN1 {_cn(spec.get('smn1_cn'))} · SMN2 {_cn(spec.get('smn2_cn'))}"
        if spec.get("smn_del78_cn"):
            text += f" · SMNΔ7-8 {spec['smn_del78_cn']}"
        return text
    if spec.get("genotype"):
        return f"genotype {spec['genotype']}"
    if "phasing_success" in spec:
        return "phasing " + ("succeeded" if spec["phasing_success"] else "failed")
    gene_cn, total_cn = spec.get("gene_cn"), region.get("total_cn")
    if gene_cn is not None and total_cn is not None:
        return f"{gene_cn} gene cop{'y' if gene_cn == 1 else 'ies'} of {total_cn} total"
    return None


def _flags(name: str, region: dict) -> list[dict]:
    flags = []
    depth = _median_depth(region)
    if not depth and not region.get("final_haplotypes"):
        flags.append({"kind": "no_data", "text": "No usable data: no depth and no phased haplotypes."})
    fusions = _fusion_names(region)
    if fusions:
        types = sorted({str((_specific(region)["fusions_called"][f] or {}).get("type", "fusion")) for f in fusions})
        flags.append({"kind": "fusion", "text": f"Fusion haplotype called ({', '.join(types)}): {', '.join(fusions)}."})
    sv = _sv_text(region)
    if sv:
        flags.append({"kind": "sv", "text": f"Structural variant called: {sv}."})
    if _specific(region).get("phasing_success") is False:
        flags.append({"kind": "phasing_failed", "text": "Paraphase reported that phasing did not succeed for this region."})
    if "smn1_cn" in _specific(region) and _specific(region).get("smn1_cn") is None:
        flags.append({"kind": "cn_unresolved", "text": "SMN1 and SMN2 copy numbers were not reported for this sample."})
    return flags


def _row(name: str, region: dict) -> dict:
    depth = _median_depth(region)
    genome = region.get("genome_depth")
    return {
        "region": name,
        "genes": [g for g in (region.get("genes_in_region") or "").split(",") if g],
        "medical": name.lower() in {r.lower() for r in runner.MEDICAL_REGIONS},
        "total_cn": region.get("total_cn"),
        "gene_cn": _specific(region).get("gene_cn"),
        "n_haplotypes": len(region.get("final_haplotypes") or {}),
        "two_copy_haplotypes": len(region.get("two_copy_haplotypes") or []),
        "depth_median": depth,
        "depth_p80": (region.get("region_depth") or {}).get("percentile80"),
        "genome_depth": genome,
        "depth_vs_genome": round(depth / genome, 2) if depth is not None and genome else None,
        "headline": _headline(name, region),
        "flags": _flags(name, region),
        "phase_region": region.get("phase_region"),
    }


def table(sample: str) -> dict:
    """Every region in the sample's JSON, plus what produced it."""
    source_schema, data = _loaded(sample)
    provenance = calls.vcf_provenance(samples.vcfs_dir(sample))
    rows = [_row(name, r) for name, r in data.items() if isinstance(r, dict)]
    sexes = {r.get("sample_sex") for r in data.values() if isinstance(r, dict) and r.get("sample_sex")}
    return {
        "sample": sample,
        "source_schema": source_schema,
        "paraphase_version": provenance["paraphase_version"],
        "paraphase_command": provenance["paraphase_command"],
        "sample_sex": sorted(sexes)[0] if len(sexes) == 1 else None,
        "n_regions": len(rows),
        "regions": rows,
    }


def _haplotype(name: str, region: dict) -> dict:
    details = (region.get("haplotype_details") or {}).get(name) or {}
    variants = details.get("variants") or []
    m = _HAP_NAME.match(name)
    label = m.group("label") if m else None
    unique = (region.get("unique_supporting_reads") or {}).get(name)
    return {
        "name": name,
        "group": label,
        "kind": _KIND_BY_LABEL.get(label) if label else None,
        "boundary": details.get("boundary"),
        "boundary_gene2": details.get("boundary_gene2"),
        "truncated": details.get("is_truncated") or [],
        "haplogroup": details.get("haplogroup"),
        "two_copy": name in (region.get("two_copy_haplotypes") or []),
        "unique_reads": len(unique) if isinstance(unique, list) else None,
        "n_variants": len(variants),
        "variants": variants[:_MAX_VARIANTS_SHOWN],
    }


def detail(sample: str, region_name: str) -> Optional[dict]:
    """One region in full: its summary row, haplotypes, alleles, fusion calls and
    every scalar region-specific field. None if the region isn't in the JSON."""
    region = _calls_for(sample).get(region_name)
    if not isinstance(region, dict):
        return None
    names = list((region.get("final_haplotypes") or {}).values())
    spec = _specific(region)
    fusions = spec.get("fusions_called")
    return {
        "sample": sample,
        "summary": _row(region_name, region),
        "haplotypes": [_haplotype(n, region) for n in names],
        "alleles": spec.get("alleles_final") if isinstance(spec.get("alleles_final"), list) else None,
        "fusions": (
            [
                {"name": n, "type": (f or {}).get("type"), "breakpoint": (f or {}).get("breakpoint")}
                for n, f in sorted(fusions.items())
            ] if isinstance(fusions, dict) else []
        ),
        # Scalar region-specific fields (copy numbers, genotype, read counts, ...),
        # shown as Paraphase wrote them; nested structures are left to the JSON.
        "region_fields": {
            k: v for k, v in sorted(spec.items())
            if isinstance(v, (int, float, str, bool)) or v is None
        },
    }
