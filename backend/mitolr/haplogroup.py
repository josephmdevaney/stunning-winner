"""Real mtDNA haplogroup classification via Haplogrep 3 (genepi/haplogrep3),
run locally against vendored PhyloTree data -- an actual phylogenetic
classification of the sample's real variant calls, not a placeholder.

Requires `python backend/scripts/setup_haplogrep.py` to have been run once
(fetches the Haplogrep 3 jar + the one PhyloTree tree this module uses; not
committed to source control, see .gitignore) and a Java 11+ runtime on PATH
(override with the HAPLOGREP_JAVA env var). Runs fully offline after that
one-time setup -- no live network call per classification request.

If the classifier isn't installed, or classification fails for any reason,
this reports that honestly via the `region` field (e.g. "not installed", or
the specific failure) rather than crashing the request or fabricating a
haplogroup call.
"""

from __future__ import annotations

import os
import subprocess
import tempfile
from functools import lru_cache
from pathlib import Path
from typing import Optional, TypedDict

from . import vcf as vcf_mod
from .reference import GENOME_LENGTH

VENDOR_DIR = Path(__file__).resolve().parent / "vendor"
JAR_PATH = VENDOR_DIR / "haplogrep3.jar"
JAVA_PATH_FILE = VENDOR_DIR / "java_path.txt"
TREE_ID = "phylotree-fu-rcrs@1.3"
JAVA_BIN_ENV = "HAPLOGREP_JAVA"


class ClassifyResult(TypedDict):
    name: str
    confidence: Optional[float]
    region: str


def available() -> bool:
    """Whether the vendored classifier is installed (see setup_haplogrep.py)."""
    return JAR_PATH.exists()


def _java_bin() -> str:
    """A Java 11+ binary to run the classifier with: HAPLOGREP_JAVA if set
    (explicit override always wins); otherwise the binary setup_haplogrep.py
    recorded as working on this machine (see JAVA_PATH_FILE) -- this is what
    makes classification work out of the box after setup, without needing
    HAPLOGREP_JAVA exported in whatever environment starts the server; falls
    back to a plain 'java' on PATH, which is correct on environments where
    that's already a working JDK (e.g. most Linux deployments)."""
    env = os.environ.get(JAVA_BIN_ENV)
    if env:
        return env
    if JAVA_PATH_FILE.exists():
        recorded = JAVA_PATH_FILE.read_text().strip()
        if recorded:
            return recorded
    return "java"


def _build_vcf(sample: str, records: list[dict]) -> str:
    """A minimal VCF Haplogrep 3 can parse: contig 'MT' (matching its own
    example data's naming, not MitoLR's usual 'chrM') with one row per real
    call, GT/DP/AD/VAF preserved so Haplogrep reads heteroplasmy fraction from
    AD rather than treating every call as a binary present/absent genotype."""
    lines = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID=MT,length={GENOME_LENGTH}>",
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}",
    ]
    for r in records:
        ad = ",".join(str(x) for x in r["ad"]) if r["ad"] else "."
        vaf = r["vaf"] if r["vaf"] is not None else "."
        dp = r["dp"] if r["dp"] is not None else "."
        sample_field = f"{r['gt']}:{dp}:{ad}:{vaf}"
        lines.append(
            "\t".join(["MT", str(r["pos"]), ".", r["ref"], r["alt"], ".", ".", ".", "GT:DP:AD:VAF", sample_field])
        )
    return "\n".join(lines) + "\n"


def _parse_output(text: str, sample: str) -> Optional[tuple[str, float]]:
    """(haplogroup, quality) from Haplogrep's tab-separated --out file, or None
    if this sample's row isn't in it."""
    lines = text.strip().splitlines()
    if len(lines) < 2:
        return None
    header = [h.strip('"') for h in lines[0].split("\t")]
    idx = {name: i for i, name in enumerate(header)}
    for line in lines[1:]:
        cols = [c.strip('"') for c in line.split("\t")]
        if cols[idx["SampleID"]] == sample:
            return cols[idx["Haplogroup"]], float(cols[idx["Quality"]])
    return None


@lru_cache(maxsize=None)
def classify(sample: str, vcf_path: str) -> ClassifyResult:
    """Real Haplogrep 3 classification for this sample's real VCF calls."""
    if not available():
        return {
            "name": "Unclassified",
            "confidence": None,
            "region": "Haplogrep 3 not installed -- run `python backend/scripts/setup_haplogrep.py`",
        }

    records = vcf_mod.parse_records(vcf_path)
    if not records:
        return {"name": "Unclassified", "confidence": None, "region": "No variant calls to classify from."}

    vcf_text = _build_vcf(sample, records)
    with tempfile.TemporaryDirectory() as tmp:
        in_path = Path(tmp) / "sample.vcf"
        out_path = Path(tmp) / "result.txt"
        in_path.write_text(vcf_text, encoding="utf-8")
        try:
            subprocess.run(
                [
                    _java_bin(), "-jar", str(JAR_PATH), "classify",
                    "--tree", TREE_ID,
                    "--in", str(in_path),
                    "--out", str(out_path),
                ],
                cwd=VENDOR_DIR,
                capture_output=True,
                text=True,
                timeout=60,
                check=True,
            )
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError) as exc:
            return {
                "name": "Unclassified",
                "confidence": None,
                "region": f"Haplogrep 3 classification failed: {exc}",
            }
        result = _parse_output(out_path.read_text(encoding="utf-8"), sample) if out_path.exists() else None

    if result is None:
        return {
            "name": "Unclassified",
            "confidence": None,
            "region": "Haplogrep 3 ran but returned no result for this sample.",
        }

    name, quality = result
    return {
        "name": name,
        "confidence": round(quality * 100, 1),
        "region": f"Haplogrep 3 · PhyloTree {TREE_ID}",
    }
