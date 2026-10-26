"""Where the fetched EpiSignaLR files live, and whether they are all there.

The files are fetched by `backend/scripts/setup_episignalr.py` into `vendor/` (gitignored: the EpiSignaLR
repository has no licence, so its files are never committed). Under it:

    VERSION.json              the commit fetched, when, and a SHA-256 for every file
    raw/<Label>/              the model's own files, as published
    raw/array/                the array-derived reference data and CpG coordinates, as published
    models/<Label>.json       the scoring model converted from the pickle (see convert.py)
    models/array_references.json   the array-derived reference tables converted from the pickle (convert_array.py)
    .venv/                    made by setup when needed: numpy, joblib and scikit-learn for the two converters
                              and the array-derived scorer
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


# The array-derived classifiers' files: reference tables (pickle), the CpGs a pileup is intersected with (bed) and
# the converted tables.
ARRAY_FILES = {
    "pickle": RAW_DIR / "array" / "no_strand_all_points_dict.pickle",
    "coordinates": RAW_DIR / "array" / "collapsed_coordinates_sorted.bed",
    "converted": MODELS_DIR / "array_references.json",
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


def array_installed() -> bool:
    """The array-derived classifiers need their files and a Python with scikit-learn to run them."""
    return info() is not None and all(p.exists() for p in ARRAY_FILES.values()) and scoring_python() is not None


def scoring_python() -> Optional[str]:
    """The interpreter setup found or made with numpy and scikit-learn, if it is still there."""
    recorded = (info() or {}).get("python")
    return recorded if recorded and Path(recorded).exists() else None
