#!/usr/bin/env python3
"""Download and locally vendor Haplogrep 3 (real mtDNA haplogroup classifier) plus
the one PhyloTree reference tree MitoLR classifies against, so real haplogroup
classification (mitolr/haplogroup.py) runs fully offline -- no live network call
per classification request, just this one-time setup fetch. Run with:

    python backend/scripts/setup_haplogrep.py

Idempotent: skips anything already present. Vendored files are NOT committed to
source control (see .gitignore) -- like a virtualenv or node_modules, this is a
fetched build dependency, not project source.
"""

from __future__ import annotations

import glob
import os
import shutil
import subprocess
import tarfile
import urllib.request
from pathlib import Path

VENDOR_DIR = Path(__file__).resolve().parents[1] / "mitolr" / "vendor"
JAVA_PATH_FILE = VENDOR_DIR / "java_path.txt"

# Common install locations for a real JDK that a plain 'java' on PATH might not
# resolve to (e.g. Homebrew's openjdk is keg-only by design -- it deliberately
# isn't symlinked onto PATH so it doesn't fight macOS's own Java shim).
_JAVA_CANDIDATES_GLOBS = [
    "/opt/homebrew/opt/openjdk*/bin/java",  # Homebrew, Apple Silicon
    "/usr/local/opt/openjdk*/bin/java",  # Homebrew, Intel
    "/usr/lib/jvm/*/bin/java",  # common Linux package locations
]

HAPLOGREP_VERSION = "3.3.2"
HAPLOGREP_URL = (
    f"https://github.com/genepi/haplogrep3/releases/download/"
    f"v{HAPLOGREP_VERSION}/haplogrep3-{HAPLOGREP_VERSION}-linux.tar.gz"
)
# The Linux/Windows release assets both bundle the same platform-independent
# haplogrep3.jar (plus a native launcher script + JRE this project doesn't use --
# mitolr/haplogroup.py runs the jar with the system/Homebrew JVM instead).

TREE_ID = "phylotree-fu-rcrs"
TREE_VERSION = "1.3"  # current default recommended by Haplogrep 3, GRCh38/rCRS-based


def _download(url: str, dest: Path) -> None:
    print(f"Downloading {url} ...")
    with urllib.request.urlopen(url) as resp, open(dest, "wb") as f:
        shutil.copyfileobj(resp, f)


def setup_jar() -> None:
    jar_path = VENDOR_DIR / "haplogrep3.jar"
    if jar_path.exists():
        print(f"Already present: {jar_path}")
        return
    VENDOR_DIR.mkdir(parents=True, exist_ok=True)
    tmp_tar = VENDOR_DIR / "_haplogrep3.tar.gz"
    _download(HAPLOGREP_URL, tmp_tar)
    with tarfile.open(tmp_tar) as tar:
        # `data/examples/*` is only needed transiently, to give setup_tree() a
        # tiny real VCF to run a bootstrap classification against below.
        for member in tar.getmembers():
            if member.name in ("haplogrep3.jar", "LICENSE") or member.name.startswith("data/examples/"):
                tar.extract(member, path=VENDOR_DIR)
    tmp_tar.unlink()
    # Trim the bundled config's `phylotrees:` list down to the one tree this
    # build actually uses, so haplogrep3 doesn't fetch the other 7 supported
    # trees it lists by default (phylotree-rcrs@17.x, phylotree-rsrs@17.x,
    # etc.) that MitoLR never classifies against. The `repositories:` key
    # MUST stay -- without it, haplogrep3 has nowhere to resolve/download an
    # unrecognized tree ID from and fails outright, even for a tree ID that's
    # listed in `phylotrees:`.
    (VENDOR_DIR / "haplogrep3.yaml").write_text(
        "repositories:\n"
        "  - https://raw.githubusercontent.com/genepi/haplogrep-trees/main/trees.json\n"
        "\n"
        "phylotrees:\n"
        f"- {TREE_ID}@{TREE_VERSION}\n"
    )
    print(f"Installed: {jar_path}")


def setup_tree(java_bin: str) -> None:
    """Haplogrep 3 tracks "installed" trees via its own internal registry file
    (not just the presence of a tree.xml on disk), so the tree must be
    installed by actually running the jar -- not by hand-extracting its
    release zip -- or it won't recognize the tree as available."""
    tree_dir = VENDOR_DIR / "trees" / TREE_ID / TREE_VERSION
    if (tree_dir / "tree.xml").exists():
        print(f"Already present: {tree_dir}")
        return
    example_vcf = VENDOR_DIR / "data" / "examples" / "example-wgs.vcf"
    tmp_out = VENDOR_DIR / "_bootstrap_out.txt"
    subprocess.run(
        [
            java_bin, "-jar", str(VENDOR_DIR / "haplogrep3.jar"), "classify",
            "--tree", f"{TREE_ID}@{TREE_VERSION}",
            "--in", str(example_vcf),
            "--out", str(tmp_out),
        ],
        cwd=VENDOR_DIR,
        check=True,
    )
    tmp_out.unlink(missing_ok=True)
    shutil.rmtree(VENDOR_DIR / "data")
    # The jar's own install re-bundles the tree's contents as package.zip
    # inside the tree dir too -- redundant once extracted, and it roughly
    # doubles this tree's footprint.
    (tree_dir / "package.zip").unlink(missing_ok=True)
    print(f"Installed: {tree_dir}")


def _java_works(java_bin: str) -> bool:
    try:
        subprocess.run([java_bin, "-version"], capture_output=True, check=True)
        return True
    except (subprocess.CalledProcessError, OSError):
        return False


def resolve_java(explicit: str | None) -> str:
    """A working Java 11+ binary, or "" if none could be found. Tries, in order:
    an explicitly-passed --java path; a plain 'java' on PATH (works out of the
    box on most Linux deployments); JAVA_HOME; and common install locations a
    plain 'java' might not resolve to (e.g. Homebrew's openjdk is deliberately
    keg-only, not symlinked onto PATH, so it doesn't fight macOS's own Java
    shim). The winner is recorded to java_path.txt so haplogroup.py finds the
    SAME binary at runtime without needing HAPLOGREP_JAVA set in every shell
    that starts the server -- see haplogroup.py's _java_bin()."""
    candidates = []
    if explicit:
        candidates.append(explicit)
    candidates.append("java")
    java_home = os.environ.get("JAVA_HOME")
    if java_home:
        candidates.append(str(Path(java_home) / "bin" / "java"))
    for pattern in _JAVA_CANDIDATES_GLOBS:
        candidates.extend(sorted(glob.glob(pattern), reverse=True))  # newest version first

    for candidate in candidates:
        if _java_works(candidate):
            return candidate
    return ""


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--java", default=None, help="Path to a Java 11+ executable (default: auto-detect)")
    args = parser.parse_args()

    java_bin = resolve_java(args.java)
    if not java_bin:
        print(
            "ERROR: no working Java 11+ found (checked --java, 'java' on PATH, JAVA_HOME, "
            "and common Homebrew/Linux install locations).\n"
            "Install one, e.g.:\n"
            "  macOS:  brew install openjdk\n"
            "  Linux:  apt install openjdk-17-jre-headless  (or your distro's equivalent)\n"
            "Then re-run this script (or pass --java /path/to/java)."
        )
        raise SystemExit(1)
    print(f"Using Java: {java_bin}")

    setup_jar()
    setup_tree(java_bin)

    if java_bin != "java":
        # Record it so haplogroup.py finds this same binary at runtime without
        # needing HAPLOGREP_JAVA set in whatever environment starts the server.
        JAVA_PATH_FILE.write_text(java_bin)
        print(f"Recorded {java_bin} to {JAVA_PATH_FILE} for mitolr/haplogroup.py to use at runtime.")

    print("Haplogrep setup complete.")
