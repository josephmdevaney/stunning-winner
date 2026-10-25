"""Where the fetched EpiSignaLR files live, and whether they are all there.

The files are fetched by `backend/scripts/setup_episignalr.py` into `vendor/` (gitignored: the EpiSignaLR
repository has no licence, so its files are never committed). Under it:

    VERSION.json              the commit fetched, when, and a SHA-256 for every file
    raw/<Label>/              the model's own files, as published
    models/<Label>.json       the scoring model converted from the pickle (see convert.py)
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

VENDOR_DIR = Path(__file__).resolve().parent / "vendor"
VERSION_FILE = VENDOR_DIR / "VERSION.json"
RAW_DIR = VENDOR_DIR / "raw"
MODELS_DIR = VENDOR_DIR / "models"

REPO = "millerlaboratory/EpiSignaLR"
SETUP_COMMAND = "python backend/scripts/setup_episignalr.py"

# key -> (file label in the repository, display name, gene the signature is anchored on)
CATALOG = {
    "setd1b": ("SETD1B", "SETD1B-related neurodevelopmental disorder", "SETD1B"),
    "kleefstra": ("Kleefstra", "Kleefstra syndrome", "EHMT1"),
    "kabuki": ("Kabuki", "Kabuki syndrome", "KMT2D"),
}

# The model's own files within raw/<Label>/
def model_files(label: str) -> dict[str, Path]:
    d = RAW_DIR / label
    return {
        "metadata": d / "enet_svm_metadata.json",
        "training": d / "enet_svm_training_data.csv",
        "loo": d / "loo_predictions.csv",
        "dmrs": d / f"{label}_enet_svm_full_panel_DMPs.bed",
        "converted": MODELS_DIR / f"{label}.json",
    }


def info() -> Optional[dict]:
    """The recorded fetch (commit, time, hashes), or None when nothing has been fetched."""
    try:
        return json.loads(VERSION_FILE.read_text())
    except (OSError, ValueError):
        return None


def installed() -> bool:
    if info() is None:
        return False
    return all(p.exists() for key in CATALOG for p in model_files(CATALOG[key][0]).values())
