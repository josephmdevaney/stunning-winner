"""Sample registry for the Epigenetics suite.

A sample is a genome's pb-CpG-tools pileup: three bgzipped, tabix-indexed beds (combined, hap1, hap2).
Which samples exist, and where their files live, comes from a local config file -- never from the repo --
because the paths are machine-specific and the genomes may be sensitive:

    $EPIGENETICS_CONFIG, else backend/config/epigenetics.local.json  (see epigenetics.local.json.example)

Each entry is either explicit file paths (`combined`, `hap1`, `hap2`) or `dir` + `prefix`, where the files
are `<dir>/<prefix>.<track>.bed.gz`. `hap1`/`hap2` are optional; without them calls fall back to the pooled
value alone.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

BACKEND_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = BACKEND_ROOT / "config" / "epigenetics.local.json"
TRACKS = ("combined", "hap1", "hap2")


@dataclass(frozen=True)
class Sample:
    key: str
    label: str
    paths: dict[str, Path]      # track -> bed path (only tracks that exist)
    note: Optional[str] = None


def config_path() -> Path:
    return Path(os.environ.get("EPIGENETICS_CONFIG", str(DEFAULT_CONFIG)))


def _expand(entry: dict) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    for track in TRACKS:
        if entry.get(track):
            paths[track] = Path(entry[track]).expanduser()
        elif entry.get("dir") and entry.get("prefix"):
            paths[track] = Path(entry["dir"]).expanduser() / f"{entry['prefix']}.{track}.bed.gz"
    return paths


def load() -> dict[str, Sample]:
    """The registered samples whose combined bed (and its .tbi index) exist. A missing optional haplotype
    track is dropped; a missing combined track drops the sample."""
    path = config_path()
    try:
        raw = json.loads(path.read_text())
    except (OSError, ValueError):
        return {}
    out: dict[str, Sample] = {}
    for key, entry in (raw.get("samples") or {}).items():
        paths = _expand(entry)
        combined = paths.get("combined")
        if combined is None or not combined.exists() or not Path(str(combined) + ".tbi").exists():
            continue
        usable = {t: p for t, p in paths.items() if p.exists() and Path(str(p) + ".tbi").exists()}
        out[key] = Sample(key=key, label=entry.get("label") or key, paths=usable, note=entry.get("note"))
    return out
