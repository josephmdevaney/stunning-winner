"""Sample registry and data shaping for the Kivvi suite.

Kivvi writes `<prefix>.kivvi.<target>.{json,vcf,svg,bam}` per target (`kiv2`,
`d4z4`); the VCF may be bgzipped (`.vcf.gz`). Samples are discovered from a results
root laid out as

    <root>/<sample>/<target>/<sample>.kivvi.<target>.json

Dropping a new Kivvi run into that layout is all it takes to add a sample; nothing
is hardcoded. The root lives outside the repo -- real per-participant results are
never committed (like the other real samples in this project) -- and defaults to
~/kivvi_tool/runs, overridable with KIVVI_RESULTS_ROOT.
"""

from __future__ import annotations

import gzip
import json
import os
from pathlib import Path
from typing import Optional

from . import parse

TARGETS = ("d4z4", "kiv2")
TARGET_LABELS = {"d4z4": "D4Z4", "kiv2": "LPA KIV2"}

RESULTS_ROOT = Path(os.environ.get("KIVVI_RESULTS_ROOT", Path.home() / "kivvi_tool" / "runs"))

# (path, mtime_ns) -> parsed file, so re-running Kivvi invalidates on its own.
_json_cache: dict[tuple, dict] = {}
_vcf_cache: dict[tuple, dict] = {}


def _target_dir(sample: str, target: str) -> Path:
    return RESULTS_ROOT / sample / target


def _file(sample: str, target: str, ext: str) -> Path:
    return _target_dir(sample, target) / f"{sample}.kivvi.{target}.{ext}"


def sample_keys() -> list[str]:
    if not RESULTS_ROOT.is_dir():
        return []
    return sorted(
        d.name
        for d in RESULTS_ROOT.iterdir()
        if d.is_dir() and any(_file(d.name, t, "json").is_file() for t in TARGETS)
    )


def default_sample() -> str:
    keys = sample_keys()
    return keys[0] if keys else ""


def resolve(sample: Optional[str]) -> str:
    """Validate a sample key (None -> default). Raises KeyError if unknown."""
    key = sample or default_sample()
    if key not in sample_keys():
        raise KeyError(key)
    return key


def has_target(sample: str, target: str) -> bool:
    return _file(sample, target, "json").is_file()


def svg_path(sample: str, target: str) -> Optional[Path]:
    path = _file(sample, target, "svg")
    return path if path.is_file() else None


def _stamp(path: Path) -> tuple:
    return (str(path), path.stat().st_mtime_ns)


def _load_json(sample: str, target: str) -> Optional[dict]:
    path = _file(sample, target, "json")
    if not path.is_file():
        return None
    key = _stamp(path)
    if key not in _json_cache:
        _json_cache[key] = json.loads(path.read_text())
    return _json_cache[key]


def _vcf_path(sample: str, target: str) -> Optional[Path]:
    """Kivvi's VCF, plain (`.vcf`) or bgzipped (`.vcf.gz`, as the GeneDx pipeline writes it)."""
    for ext in ("vcf", "vcf.gz"):
        path = _file(sample, target, ext)
        if path.is_file():
            return path
    return None


def _read_vcf(path: Path) -> str:
    # bgzip is a series of concatenated gzip members, which gzip reads as one stream.
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as f:
            return f.read()
    return path.read_text()


def _load_vcf(sample: str, target: str) -> Optional[dict]:
    path = _vcf_path(sample, target)
    if path is None:
        return None
    key = _stamp(path)
    if key not in _vcf_cache:
        _vcf_cache[key] = parse.parse_vcf(_read_vcf(path))
    return _vcf_cache[key]


def _genome_depth(sample: str) -> Optional[float]:
    for target in TARGETS:
        data = _load_json(sample, target)
        depth = (data or {}).get("depth_summary", {}).get("genome_depth")
        if depth is not None:
            return depth
    return None


def meta(sample: str) -> dict:
    version = None
    for target in TARGETS:
        vcf = _load_vcf(sample, target)
        if vcf and vcf["version"]:
            version = vcf["version"]
            break
    # Kivvi can't report genome depth from a bamlet; surface that rather than hide it.
    bamlet = _genome_depth(sample) is None
    return {
        "sample": sample,
        "reference": "GRCh38",
        "kivvi_version": version,
        "targets": {t: has_target(sample, t) for t in TARGETS},
        "input_note": "Bamlet input: genome-wide depth unavailable" if bamlet else None,
    }


def _allele_labels(data: dict, shaped: dict, columns: list[str]) -> list[str]:
    """Human labels for the VCF's allele columns.

    The VCF's sample columns follow ``complete_alleles`` order (not allele_info
    order), so each column is matched to its allele by name to show its size.
    """
    by_name = {a["name"]: a for a in shaped["alleles"]}
    complete = data.get("complete_alleles", [])
    if len(complete) != len(columns):
        return list(columns)
    labels = []
    for i, name in enumerate(complete, start=1):
        allele = by_name.get(name)
        if allele is None:
            labels.append(f"A{i}")
        elif "chromosome" in allele:  # D4Z4
            labels.append(f"A{i} · {allele['size_label']} · {allele['chromosome']}")
        else:  # KIV2
            labels.append(f"A{i} · {allele['copy_number']}")
    return labels


def _vcf_variants(sample: str, target: str, data: dict, shaped: dict) -> dict:
    vcf = _load_vcf(sample, target)
    if vcf is None:
        return {"allele_columns": [], "allele_labels": [], "records": [], "n_total": 0, "n_on_assembled": 0}
    records = vcf["records"]
    return {
        "allele_columns": vcf["allele_columns"],
        "allele_labels": _allele_labels(data, shaped, vcf["allele_columns"]),
        "records": records,
        "n_total": len(records),
        "n_on_assembled": sum(1 for r in records if r["on_assembled_allele"]),
    }


def target_detail(sample: str, target: str) -> Optional[dict]:
    data = _load_json(sample, target)
    if data is None:
        return None
    shaped = parse.shape_d4z4(data) if target == "d4z4" else parse.shape_kiv2(data)
    return {
        "sample": sample,
        "target": target,
        "label": TARGET_LABELS[target],
        **shaped,
        "variants": _vcf_variants(sample, target, data, shaped),
        "svg_available": svg_path(sample, target) is not None,
    }


def overview(sample: str) -> dict:
    d4z4 = target_detail(sample, "d4z4")
    kiv2 = target_detail(sample, "kiv2")
    return {
        **meta(sample),
        "d4z4": None
        if d4z4 is None
        else {
            "summary": d4z4["summary"],
            "overall_methylation": d4z4["overall_methylation"],
            "alleles": [
                {
                    "size_label": a["size_label"],
                    "chromosome": a["chromosome"],
                    "upstream_group": a["upstream_group"],
                    "distal_haplotype": a["distal_haplotype"],
                    "methylation": a["methylation"],
                    "interpretation": a["interpretation"],
                }
                for a in d4z4["alleles"]
            ],
        },
        "kiv2": None
        if kiv2 is None
        else {
            "allele_cn": kiv2["allele_cn"],
            "shortest_copy_number": kiv2["shortest_copy_number"],
            "alleles": [{"index": a["index"], "copy_number": a["copy_number"]} for a in kiv2["alleles"]],
            "notes": kiv2["notes"],
        },
    }


def sample_qc(sample: str) -> dict:
    targets = []
    for target in TARGETS:
        data = _load_json(sample, target)
        if data is None:
            continue
        depth = data.get("depth_summary") or {}
        qc = (data.get("additional") or {}).get("qc_metrics") or {}
        detail = target_detail(sample, target)
        expected = parse.EXPECTED_D4Z4_ALLELES if target == "d4z4" else parse.EXPECTED_KIV2_ALLELES
        resolved = len(detail["alleles"])
        targets.append(
            {
                "target": target,
                "label": TARGET_LABELS[target],
                "repeat_depth": depth.get("repeat_depth"),
                "genome_depth": depth.get("genome_depth"),
                "median_read_length": qc.get("median_read_length"),
                "per_allele_depth": qc.get("per_allele_depth"),
                "reads_linking_units": len(data.get("read_info") or {}),
                "alleles_resolved": resolved,
                "alleles_expected": expected,
                "partial_alleles": len(data.get("partial_alleles") or []),
                "variants_total": detail["variants"]["n_total"],
                "variants_on_assembled": detail["variants"]["n_on_assembled"],
            }
        )
    return {**meta(sample), "target_qc": targets}
