"""Single-cell expression-TR (sc-eTR) annotations for the STRchive disease panel.

The TenK10K paper (Tanudisastro et al., bioRxiv 10.1101/2024.11.02.621562) tested how
common repeat-length variation relates to gene expression in 28 immune cell types, in
1,925 people. This module reads the small per-locus file that
`scripts/real_data/build_sc_etr_disease_annotation.py` derives from the paper's released
tables: for each STRchive disease locus, whether the repeat was polymorphic in that cohort
and which repeats at or near it were the lead expression association for a gene.

This is a research annotation about common-range length variation and expression. It says
nothing about pathogenic expansions, and an absent entry means "not a lead association",
never "no effect" -- the released tables carry only the lead repeat per gene per cell type.

The paper and its data are CC BY-NC-ND 4.0, so the derived file is not part of the repo; it
is read from a path outside it (see `SC_ETR_ANNOTATIONS_JSON` in the API) and the feature
is off when the file is absent.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

SCHEMA_VERSION = 1


def load_sc_etr_annotations(path: str | Path) -> dict:
    """Parse the annotation file. Raises ValueError (including JSON errors) if it is not a
    file this code understands, so a stale or truncated build is reported, not guessed at."""
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(data, dict) or data.get("schema_version") != SCHEMA_VERSION:
        found = data.get("schema_version") if isinstance(data, dict) else type(data).__name__
        raise ValueError(f"sc-eTR annotation file has schema_version {found!r}, expected {SCHEMA_VERSION}")
    if not isinstance(data.get("loci"), dict) or not isinstance(data.get("source"), dict):
        raise ValueError("sc-eTR annotation file is missing its 'loci' or 'source' section")
    return data


def locus_annotation(data: dict, locus_id: str) -> Optional[dict]:
    """The annotation for one STRchive locus id, or None if the file has no entry for it."""
    entry = data["loci"].get(locus_id)
    if entry is None:
        return None
    return {
        "available": True,
        "polymorphic": bool(entry["polymorphic"]),
        "leads": entry["leads"],
        "source": data["source"],
    }
