#!/usr/bin/env python3
"""Download and locally vendor PacBio's `interpret_paraphase.py` (the Paraphase
project's own annotation script, see
https://github.com/PacificBiosciences/paraphase/tree/main/annotation) plus the
two small data files it needs, so paraphase_annotate/runner.py can invoke it.

    python backend/scripts/setup_paraphase_annotate.py             # the pinned tag below
    python backend/scripts/setup_paraphase_annotate.py --latest    # newest GitHub release
    python backend/scripts/setup_paraphase_annotate.py --ref v4.1.0 --force

What was fetched, from which tag and when, is recorded in
`paraphase_annotate/vendor/VERSION.json`; the runner reads it, so an upgrade is
visible in the app and a different tag is re-fetched (not skipped) on the next
run of this script.

Vendored files are NOT committed to source control (see .gitignore) -- like
haplogrep3.jar (see setup_haplogrep.py), this is a fetched build dependency,
not project source. That distinction matters here specifically: Paraphase ships
under PacBio's own proprietary EULA (LICENSE.md in that repo), which grants
"distribute in binary form only" -- republishing its source in this (public)
repo's git history would not be a license-compliant redistribution, so each
machine fetches its own copy instead.

Version notes:
  * v4.x reads each region's region-specific fields from a `region_specific_info`
    sub-dict (Paraphase's v4.0.0 Rust rewrite changed the JSON layout). The
    runner adapts older v3.x JSON to that layout (paraphase_annotate/calls.py),
    so one v4.x script serves output from either generation of Paraphase.
  * v4.1.0's script uses an f-string grammar (PEP 701) that needs Python 3.12+.
    The runner launches it with a separate interpreter (it finds one, or use
    PARAPHASE_ANNOTATE_PYTHON); the API server itself can stay on an older Python.
  * v3.5.0 is still fetchable with --ref v3.5.0 (plain Python 3.9, flat JSON only).
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

VENDOR_DIR = Path(__file__).resolve().parents[1] / "paraphase_annotate" / "vendor"
DATA_DIR = VENDOR_DIR / "data"
VERSION_FILE = VENDOR_DIR / "VERSION.json"

DEFAULT_TAG = "v4.1.0"
REPO = "PacificBiosciences/paraphase"


def _major(tag: str) -> int:
    return int(tag.lstrip("v").split(".")[0])


def files_for(tag: str) -> dict[str, Path]:
    """raw.githubusercontent URL -> destination, for one release tag."""
    base = f"https://raw.githubusercontent.com/{REPO}/{tag}"
    # The v4.0.0 Rust rewrite flattened the repo layout: the GRCh38 config moved
    # from `paraphase/data/38/` to `data/38/`. Saved as config_grch38.yaml either
    # way -- this project only ever runs GRCh38.
    config = f"{base}/data/38/config.yaml" if _major(tag) >= 4 else f"{base}/paraphase/data/38/config.yaml"
    return {
        f"{base}/annotation/interpret_paraphase.py": VENDOR_DIR / "interpret_paraphase.py",
        f"{base}/annotation/data/paraphase_genes_transcript_start_end.txt": DATA_DIR / "paraphase_genes_transcript_start_end.txt",
        config: DATA_DIR / "config_grch38.yaml",
    }


def latest_tag() -> str:
    with urllib.request.urlopen(f"https://api.github.com/repos/{REPO}/releases/latest", timeout=20) as resp:
        return json.load(resp)["tag_name"]


def _download(url: str, dest: Path) -> str:
    print(f"Downloading {url} ...")
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=60) as resp:
        data = resp.read()
    dest.write_bytes(data)
    return hashlib.sha256(data).hexdigest()


def installed_tag() -> str | None:
    try:
        return json.loads(VERSION_FILE.read_text())["paraphase_tag"]
    except (OSError, ValueError, KeyError):
        return None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--ref", default=None, help=f"Paraphase release tag to vendor (default {DEFAULT_TAG})")
    ap.add_argument("--latest", action="store_true", help="vendor the newest GitHub release instead of the pinned tag")
    ap.add_argument("--force", action="store_true", help="re-download even if this tag is already vendored")
    args = ap.parse_args()

    tag = latest_tag() if args.latest else (args.ref or DEFAULT_TAG)
    current = installed_tag()
    if current == tag and not args.force and all(p.exists() for p in files_for(tag).values()):
        print(f"Paraphase annotation {tag} already vendored.")
        return 0
    if current and current != tag:
        print(f"Replacing vendored {current} with {tag}.")

    hashes = {}
    for url, dest in files_for(tag).items():
        hashes[dest.name] = _download(url, dest)
    (VENDOR_DIR / "__init__.py").touch()
    VERSION_FILE.write_text(json.dumps({
        "paraphase_tag": tag,
        "fetched_at": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "source": f"https://github.com/{REPO}/tree/{tag}/annotation",
        "sha256": hashes,
    }, indent=2) + "\n")
    print(f"Paraphase annotation setup complete ({tag}).")
    if _major(tag) >= 4:
        print("Note: this script needs Python 3.12+ with `requests` and `pyyaml` (see PARAPHASE_ANNOTATE_PYTHON).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
