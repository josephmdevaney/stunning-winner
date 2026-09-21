#!/usr/bin/env python3
"""Download and locally vendor PacBio's `interpret_paraphase.py` (the Paraphase
project's own proof-of-concept annotation script, see
https://github.com/PacificBiosciences/paraphase/tree/main/annotation) plus the
two small data files it needs, so paraphase_annotate/runner.py can invoke it. Run
with:

    python backend/scripts/setup_paraphase_annotate.py

Idempotent: skips anything already present. Vendored files are NOT committed to
source control (see .gitignore) -- like haplogrep3.jar (see setup_haplogrep.py),
this is a fetched build dependency, not project source. That distinction matters
here specifically: Paraphase ships under PacBio's own proprietary EULA
(LICENSE.md in that repo), which grants "distribute in binary form only" --
republishing its source in this (public) repo's git history would not be a
license-compliant redistribution, so each machine fetches its own copy instead.

Pinned to the v3.5.0 tag, not main: main's interpret_paraphase.py reads
region-specific JSON fields from a `region_specific_info` sub-key that
Paraphase only introduced in its v4.0.0 Rust rewrite. Real hackathon JSON this
project has on hand predates that schema change (flat fields on the region
dict), and v3.5.0's script matches it -- confirmed against a real sample's
Paraphase output. v3.5.0 also happens to be plain Python-3.9-compatible; main's
script needs Python 3.12+ (PEP 701 f-string grammar) for an unrelated reason.
"""

from __future__ import annotations

import urllib.request
from pathlib import Path

VENDOR_DIR = Path(__file__).resolve().parents[1] / "paraphase_annotate" / "vendor"
DATA_DIR = VENDOR_DIR / "data"

RAW_BASE = "https://raw.githubusercontent.com/PacificBiosciences/paraphase/v3.5.0"
FILES = {
    f"{RAW_BASE}/annotation/interpret_paraphase.py": VENDOR_DIR / "interpret_paraphase.py",
    f"{RAW_BASE}/annotation/data/paraphase_genes_transcript_start_end.txt": DATA_DIR / "paraphase_genes_transcript_start_end.txt",
    # v3.5.0 predates the v4.0.0 Rust rewrite's repo-layout flattening, so its
    # own default config still lives under a `paraphase/` package dir (not
    # `data/38/config.yaml` at repo root, where main's copy is). Renamed here
    # to config_grch38.yaml -- this project only ever runs GRCh38, and the
    # plain name avoids implying a build-number directory structure that
    # doesn't exist here.
    f"{RAW_BASE}/paraphase/data/38/config.yaml": DATA_DIR / "config_grch38.yaml",
}


def _download(url: str, dest: Path) -> None:
    print(f"Downloading {url} ...")
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url) as resp, open(dest, "wb") as f:
        f.write(resp.read())


if __name__ == "__main__":
    for url, dest in FILES.items():
        if dest.exists():
            print(f"Already present: {dest}")
            continue
        _download(url, dest)
    (VENDOR_DIR / "__init__.py").touch()
    print("Paraphase annotation setup complete.")
