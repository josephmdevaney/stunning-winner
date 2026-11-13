"""Reading and normalizing a Paraphase JSON, across Paraphase's two output schemas.

Paraphase v3.x writes one flat dict per region: the general fields (haplotypes,
depth, ...) and the fields only some regions have (`smn1_cn`, `fusions_called`,
`annotated_alleles`, ...) sit side by side. Paraphase v4.0.0 (the Rust rewrite)
keeps the general fields where they were and moves every region-specific field
under a `region_specific_info` sub-dict (see Paraphase's docs/json.md).

PacBio's own `interpret_paraphase.py` follows that change: v3.5.0's script reads
the flat fields, v4.x's reads `region_specific_info`. `to_v4()` lets the newest
script run on either generation of output, so one vendored script covers every
sample regardless of which Paraphase version produced it.

The rule is the one docs/json.md states: the general fields stay at the top
level, and everything else is region-specific. A key not listed below is
therefore treated as region-specific -- which is also what keeps this correct
when a newer Paraphase adds another per-region field.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Optional

# Paraphase docs/json.md "General fields" and "Information on phased haplotypes",
# plus `highest_total_cn` (documented there as an internal-use field).
GENERAL_FIELDS = frozenset({
    "total_cn",
    "two_copy_haplotypes",
    "phase_region",
    "genome_depth",
    "region_depth",
    "sample_sex",
    "genes_in_region",
    "region_specific_info",
    "sites_for_phasing",
    "final_haplotypes",
    "haplotype_details",
    "read_details",
    "unique_supporting_reads",
    "nonunique_supporting_reads",
    "assembled_haplotypes",
    "het_sites_not_used_in_phasing",
    "heterozygous_sites",
    "homozygous_sites",
    "highest_total_cn",
})

SCHEMA_V3 = "v3"
SCHEMA_V4 = "v4"


def load_calls(path: Path) -> dict[str, dict]:
    """The per-region dict from a Paraphase JSON, region name -> region dict."""
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise ValueError(f"{path}: expected a JSON object of regions, got {type(data).__name__}")
    return data


def schema(calls: dict[str, dict]) -> Optional[str]:
    """`v4` when any region carries `region_specific_info`, `v3` when regions
    exist but none do, None for an empty file."""
    regions = [r for r in calls.values() if isinstance(r, dict)]
    if not regions:
        return None
    return SCHEMA_V4 if any("region_specific_info" in r for r in regions) else SCHEMA_V3


def to_v4(calls: dict[str, dict]) -> dict[str, dict]:
    """A v4-shaped copy of `calls`. Regions already in the v4 shape are returned
    unchanged; a v3 region has its non-general fields moved under
    `region_specific_info`. The input is not modified."""
    out: dict[str, dict] = {}
    for name, region in calls.items():
        if not isinstance(region, dict) or "region_specific_info" in region:
            out[name] = region
            continue
        general = {k: v for k, v in region.items() if k in GENERAL_FIELDS}
        specific = {k: v for k, v in region.items() if k not in GENERAL_FIELDS}
        out[name] = {**general, "region_specific_info": specific}
    return out


_VERSION_LINE = re.compile(r"^##paraphase_version=(.+)$")
_COMMAND_LINE = re.compile(r"^##paraphase_command=(.+)$")


def vcf_provenance(vcfs_dir: Path) -> dict:
    """Which Paraphase produced this output, read from the header of the first
    per-region VCF in `<sample>_paraphase_vcfs/` (the JSON itself records no
    version). Empty values when the directory or the header lines are missing."""
    info: dict = {"paraphase_version": None, "paraphase_command": None}
    try:
        vcfs = sorted(vcfs_dir.glob("*.vcf"))
    except OSError:
        return info
    for vcf in vcfs[:1]:
        with open(vcf, encoding="utf-8", errors="replace") as f:
            for line in f:
                if not line.startswith("##"):
                    break
                m = _VERSION_LINE.match(line.strip())
                if m:
                    info["paraphase_version"] = m.group(1)
                m = _COMMAND_LINE.match(line.strip())
                if m:
                    info["paraphase_command"] = m.group(1)
    return info
