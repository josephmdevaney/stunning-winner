"""Scoring a genome against EpiSignaLR's 34 array-derived classifiers.

The classifiers are NSBEpi's (Geysens et al.), adapted by EpiSignaLR: for each of 34 conditions a linear SVM
separates that condition's published median array methylation profile from the other profiles and controls, over
that condition's 90 to 150 signature CpGs (3,642 distinct CpGs in all). A genome's methylation at those CpGs
is read from its pileup, as EpiSignaLR does -- the percent methylated of the pileup row at the CpG, no coverage
filter -- and scored by array_svm.py (a separate process: it needs scikit-learn, which the API does not).
"""

from __future__ import annotations

import json
import subprocess
import tempfile
from functools import lru_cache
from pathlib import Path

from .. import cpg
from ..samples import Sample
from . import fetched

CONDITIONS_FILE = Path(__file__).resolve().parent / "data" / "array_conditions.json"
ARRAY_SVM = Path(__file__).resolve().parent / "array_svm.py"

# EpiSignaLR's script and methods use 0.35; the preprint's results section says 0.3.
THRESHOLD = 0.35
CONTROL_PROFILE = "Control"
# The scikit-learn the EpiSignaLR environment pins; scores agree across versions (the libsvm in it has not changed).
REFERENCE_SKLEARN = "1.8.0"


class ScoringError(RuntimeError):
    pass


@lru_cache(maxsize=1)
def conditions() -> dict:
    """{id: {name, genes}} for the 34 conditions, from the preprint's Supplemental Table 4."""
    return json.loads(CONDITIONS_FILE.read_text())["conditions"]


@lru_cache(maxsize=1)
def references() -> dict:
    return json.loads(fetched.ARRAY_FILES["converted"].read_text())


def sites() -> list[tuple[str, int]]:
    """The distinct signature CpGs (chromosome without "chr", 0-based start), in file order."""
    seen: dict[tuple[str, int], None] = {}
    for cond in references()["conditions"].values():
        for chrom, start in cond["coords"]:
            seen.setdefault((chrom, start), None)
    return list(seen)


def genome_values(bed: Path) -> dict[str, float | None]:
    """The genome's methylation fraction (0-1) at every signature CpG, None where the pileup has no row."""
    found = cpg.fetch_sites(bed, sites())
    return {f"{c}:{s}": (found[(c, s)].mod_score / 100 if (c, s) in found else None) for c, s in sites()}


def run_scorer(values: dict[str, float | None]) -> dict:
    python = fetched.scoring_python()
    if python is None:
        raise ScoringError("the Python made by setup (with scikit-learn) is missing; re-run " + fetched.SETUP_COMMAND)
    with tempfile.TemporaryDirectory() as cwd:
        done = subprocess.run([python, "-I", str(ARRAY_SVM), str(fetched.ARRAY_FILES["converted"])], input=json.dumps(values),
                              capture_output=True, text=True, timeout=300, cwd=cwd)
    if done.returncode != 0:
        raise ScoringError(f"the array-derived scorer failed: {done.stderr.strip().splitlines()[-1] if done.stderr.strip() else done.returncode}")
    return json.loads(done.stdout)


def _shape(scores: dict[str, float], values: dict[str, float | None]) -> list[dict]:
    meta = conditions()
    out = []
    for cid, cond in references()["conditions"].items():
        profiles = cond["profiles"]
        j_cond, j_ctl = profiles.index(cond["positive_profile"]), profiles.index(CONTROL_PROFILE)
        keys = [f"{c}:{s}" for c, s in cond["coords"]]
        sample = [values[k] for k in keys]
        found = sum(v is not None for v in sample)
        out.append({
            "id": cid, "name": meta[cid]["name"], "genes": meta[cid]["genes"],
            "score": scores[cid], "matches": scores[cid] >= THRESHOLD,
            "n_cpgs": len(keys), "n_found": found,
            "cpgs": {
                "coords": [[c, s + 1] for c, s in cond["coords"]],      # 1-based positions for display
                "control": [row[j_ctl] for row in cond["values"]],
                "condition": [row[j_cond] for row in cond["values"]],
                "sample": sample,
            },
        })
    return sorted(out, key=lambda c: -c["score"])


def _provenance_warnings(provenance: dict) -> list[dict]:
    mode, sites_mode = provenance.get("pileup_mode"), provenance.get("modsites_mode")
    if mode is None and sites_mode is None:
        text = ("The pileup does not record how it was made. EpiSignaLR reads count-mode pileups at reference CpG sites; "
                "other modes give different methylation values.")
    elif (mode, sites_mode) != ("count", "reference"):
        text = (f"This pileup was made in {mode or 'unknown'} mode with {sites_mode or 'unknown'} sites; EpiSignaLR uses count mode with "
                "reference sites, so these scores are not directly comparable with the published ones.")
    else:
        return []
    return [{"id": "pileup_mode", "severity": "caution", "text": text}]


_cache: dict[tuple, dict] = {}


def overview(sample: Sample) -> dict:
    """All 34 conditions scored for a registered genome (cached until its pileup changes)."""
    combined = sample.paths["combined"]
    stamp = (sample.key, str(combined), combined.stat().st_mtime_ns)
    if stamp in _cache:
        return _cache[stamp]
    values = genome_values(combined)
    result = run_scorer(values)
    shaped = _shape(result["scores"], values)
    matches = [c["id"] for c in shaped if c["matches"]]
    out = {
        "sample": sample.key, "label": sample.label,
        "threshold": THRESHOLD,
        "call": "control" if not matches else "one" if len(matches) == 1 else "several",
        "matches": matches,
        "conditions": shaped,
        "cpgs": {"total": len(values), "found": sum(v is not None for v in values.values())},
        "scorer": {"scikit_learn": result["sklearn"], "reference_scikit_learn": REFERENCE_SKLEARN},
        "warnings": _provenance_warnings(cpg.header(combined)),
    }
    _cache[stamp] = out
    return out
