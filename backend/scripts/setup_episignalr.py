#!/usr/bin/env python3
"""Fetch the EpiSignaLR episignature models and convert them for the Epigenetics suite.

    python backend/scripts/setup_episignalr.py            # the pinned commit below
    python backend/scripts/setup_episignalr.py --force    # fetch and convert again

EpiSignaLR (https://github.com/millerlaboratory/EpiSignaLR, Miller laboratory) scores a long-read genome
against three condition-specific methylation signatures (SETD1B-related disorder, Kleefstra syndrome, Kabuki
syndrome) and against 34 array-derived classifiers adapted from NSBEpi. Its repository publishes no licence, so nothing from it is committed here: each machine
fetches its own copy into `epigenetics/episignature/vendor/` (gitignored), like the Paraphase annotation
script (see setup_paraphase_annotate.py).

What happens:
  1. Each file is downloaded from one pinned commit and its SHA-256 must match the value pinned below,
     so what is converted and scored is exactly what was reviewed.
  2. The models and the array-derived reference tables are pickles, and a pickle can run code when it is
     loaded. They are never unpickled in the API process: `convert.py` and `convert_array.py` read them in a
     separate `python -I` process that allows only numpy's array reconstruction and stands in for every
     scikit-learn and pandas class, and write plain JSON.
  3. The converted models are checked against the weights the authors publish in each model's metadata, and
     the array reference tables against the CpG coordinates EpiSignaLR intersects a pileup with.
  4. The array-derived classifiers are retrained for every genome (see array_svm.py), so scoring them needs
     scikit-learn. The converters and that scorer run in a Python with numpy, joblib and scikit-learn: the one
     running this script if it has them, otherwise a small virtualenv created inside the vendor folder. Its
     path is recorded in VERSION.json.

To move to a newer EpiSignaLR, change COMMIT and the hashes after reading the diff.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path

BACKEND = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BACKEND))

from epigenetics.episignature import fetched  # noqa: E402

COMMIT = "870e45be06d558d519ebf8e5acd5bf8092d4ff0c"      # main, 2026-10-21 ("Removed mention of R...")

# label -> {file name in data/<label>_enet_svm_model/: sha256}
SHA256 = {
    "Kabuki": {
        "Kabuki_enet_svm_full_panel_DMPs.bed": "24b28a4ff80c9595141eb013905875461ed47ce5525587adb3745c3bc622b472",
        "enet_svm_metadata.json": "0c47ceca6198f040f6c26d0ae6de2ddb2883545eb28fb31948c15a2bb3bd2b82",
        "enet_svm_model.joblib": "b74c0bbf9b5c33b17ccf71ba7bbff002b538054acc9087ba1fcd3b140eee34f0",
        "enet_svm_training_data.csv": "01bfb2fdb47274d071afc709b13817eb0d8d31ee11d98b3527d4c55ae2a9ddbd",
        "loo_predictions.csv": "a84b6980c58ece00c9baba493421fa6cf9a37f44149cfafe077d1bc52ff745cb",
    },
    "Kleefstra": {
        "Kleefstra_enet_svm_full_panel_DMPs.bed": "456b1fb463975632b0df97c7ee65b3e4a5c7ee75be2975691eda94edec1dc82a",
        "enet_svm_metadata.json": "b003650bf1c43e2e2dfeea18472d21158eae936c2f695b36f0f36d69f6939b80",
        "enet_svm_model.joblib": "3411d078c8c8578cd273a550525edda3c6346997e0f65ed589f9c80d221aaadf",
        "enet_svm_training_data.csv": "2fd1fef3bf3bf1b3edbfe21f47ffe59f1d0215eca1c519bc662c567834bd1475",
        "loo_predictions.csv": "b1c96522a92a996527c0293b1fb6624763fdb3d641700ca242404ef77f1163ec",
    },
    "SETD1B": {
        "SETD1B_enet_svm_full_panel_DMPs.bed": "c76ce4276e6735954cfd0d8d41c6054fb9ec14d5b30290bd94e68b2a2bb06b22",
        "enet_svm_metadata.json": "097ea6500cfb2671ddeaeb24ac246e551adccc9ddf9d2fd2e767363766d830c3",
        "enet_svm_model.joblib": "376627389598638b2200ab1f4f56754279c824bc0c412a4f6b3b4f9d8973cdb1",
        "enet_svm_training_data.csv": "0a6d8623de7257030d8e44f96f7cc0deb09733372f3706078814d370c02d5474",
        "loo_predictions.csv": "703feb13eaad202d274c31b6eb2aeef0113c97bd6cecabe3f451a2630e0427fc",
    },
}
# The array-derived classifiers' files, under data/ in the repository.
ARRAY_SHA256 = {
    "collapsed_coordinates_sorted.bed": "1a93af15230f2e03caa9293b6b4fb6eb864f25f1e0a0dde2866a406d348194df",
    "no_strand_all_points_dict.pickle": "1bacef610b69b00c122bb0e8a0d2f1d1058ba5bb69ff0b887074c556856e91e3",
}

EPISIGNATURE = BACKEND / "epigenetics" / "episignature"
CONVERT = EPISIGNATURE / "convert.py"
CONVERT_ARRAY = EPISIGNATURE / "convert_array.py"
ARRAY_SVM = EPISIGNATURE / "array_svm.py"
VENV = fetched.VENDOR_DIR / ".venv"
PACKAGES = ("numpy", "joblib", "scikit-learn")


def download(repo_path: str, sha256: str) -> bytes:
    url = f"https://raw.githubusercontent.com/{fetched.REPO}/{COMMIT}/{repo_path}"
    print(f"Downloading {repo_path} ...")
    with urllib.request.urlopen(url, timeout=60) as resp:
        data = resp.read()
    got = hashlib.sha256(data).hexdigest()
    if got != sha256:
        raise SystemExit(f"{repo_path}: SHA-256 is {got}, expected {sha256}. "
                         "Nothing was written for this file; review the change before updating the pin.")
    return data


def can_run(python: str) -> bool:
    probe = "import numpy, joblib, sklearn"
    return subprocess.run([python, "-I", "-c", probe], capture_output=True).returncode == 0


def tools_python() -> str:
    """A Python with numpy, joblib and scikit-learn: this one if it has them, else the virtualenv made for the purpose."""
    if can_run(sys.executable):
        return sys.executable
    venv_python = VENV / "bin" / "python"
    if not venv_python.exists():
        print(f"Creating {VENV} ...")
        subprocess.run([sys.executable, "-m", "venv", str(VENV)], check=True)
    if not can_run(str(venv_python)):
        print(f"Installing {', '.join(PACKAGES)} into {VENV} ...")
        subprocess.run([str(venv_python), "-m", "pip", "install", "-q", *PACKAGES], check=True)
    return str(venv_python)


def convert(python: str, label: str) -> None:
    files = fetched.model_files(label)
    files["converted"].parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as cwd:      # not the data folder: -I keeps the cwd off the import path anyway
        subprocess.run([python, "-I", str(CONVERT), str(fetched.RAW_DIR / label / "enet_svm_model.joblib"), str(files["converted"])],
                       check=True, cwd=cwd)


def convert_arrays(python: str) -> None:
    fetched.ARRAY_FILES["converted"].parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as cwd:
        subprocess.run([python, "-I", str(CONVERT_ARRAY), str(fetched.ARRAY_FILES["pickle"]), str(fetched.ARRAY_FILES["coordinates"]),
                        str(fetched.ARRAY_FILES["converted"])], check=True, cwd=cwd)


def smoke_test_scorer(python: str) -> str:
    """Score a genome with no values at all, which exercises scikit-learn end to end; returns its version."""
    done = subprocess.run([python, "-I", str(ARRAY_SVM), str(fetched.ARRAY_FILES["converted"])], input="{}", capture_output=True, text=True)
    if done.returncode != 0:
        raise SystemExit(f"the array-derived scorer does not run with {python}:\n{done.stderr}")
    result = json.loads(done.stdout)
    if len(result["scores"]) != 34:
        raise SystemExit(f"the array-derived scorer returned {len(result['scores'])} conditions, expected 34")
    return result["sklearn"]


def complete() -> bool:
    info = fetched.info()
    return bool(info) and info.get("commit") == COMMIT and fetched.installed() and fetched.array_installed()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--force", action="store_true", help="fetch and convert again even if the pinned commit is already in place")
    args = ap.parse_args()

    if complete() and not args.force:
        print(f"EpiSignaLR {COMMIT[:7]} already in place.")
        return 0

    for label, files in SHA256.items():
        for name, digest in files.items():
            dest = fetched.RAW_DIR / label / name
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(download(f"data/{label}_enet_svm_model/{name}", digest))
    for name, digest in ARRAY_SHA256.items():
        dest = fetched.RAW_DIR / "array" / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(download(f"data/{name}", digest))

    python = tools_python()
    for label in SHA256:
        convert(python, label)
    convert_arrays(python)

    from epigenetics.episignature import model  # loads each model and cross-checks it against its metadata
    model.load.cache_clear()
    for key in fetched.CATALOG:
        m = model.load(key)
        print(f"{m.label}: {len(m.features)} DMRs, {len(m.selected)} selected; weights match the published metadata.")
    sklearn_version = smoke_test_scorer(python)
    print(f"Array-derived classifiers: 34 conditions, scored with scikit-learn {sklearn_version} ({python}).")

    fetched.VERSION_FILE.write_text(json.dumps({
        "repo": fetched.REPO,
        "commit": COMMIT,
        "fetched_at": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "source": f"https://github.com/{fetched.REPO}/tree/{COMMIT}",
        "licence": "none published in the repository; files are fetched for local use and never committed",
        "python": python,
        "sha256": {**{f"{label}/{name}": digest for label, files in SHA256.items() for name, digest in files.items()},
                   **{f"array/{name}": digest for name, digest in ARRAY_SHA256.items()}},
    }, indent=2) + "\n")
    print(f"EpiSignaLR setup complete ({COMMIT[:7]}).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
