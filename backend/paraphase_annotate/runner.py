"""Runs PacBio's own `interpret_paraphase.py` (vendored, see
scripts/setup_paraphase_annotate.py) against a sample's real Paraphase JSON +
per-gene VCFs, and tracks that run's status.

The script makes live Ensembl VEP + NCBI ClinVar API calls per variant (rate
limited to ~3/sec by the script itself), so a run can take anywhere from
seconds to a few minutes depending on variant count -- too slow to run inline
in a request handler. It's launched as a background subprocess instead; the
API layer (router.py) polls status() rather than blocking on it.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import Optional

from . import samples

BACKEND_ROOT = Path(__file__).resolve().parents[1]
VENDOR_DIR = Path(__file__).resolve().parent / "vendor"
VENDOR_SCRIPT = VENDOR_DIR / "interpret_paraphase.py"
VENDOR_CONFIG = VENDOR_DIR / "data" / "config_grch38.yaml"
# Pinned to the Paraphase v3.5.0 tag (see setup_paraphase_annotate.py), not
# main -- main's interpret_paraphase.py reads region-specific JSON fields
# from a `region_specific_info` sub-key that Paraphase only introduced in its
# v4.0.0 Rust rewrite. Real hackathon JSON on hand so far predates that
# schema change (flat fields on the region dict), and v3.5.0's script matches
# it. v3.5.0 also happens to sidestep a separate Python-version issue: main's
# script uses an f-string with a nested double-quoted dict key, valid only
# under PEP 701 (Python 3.12+) -- this project's backend venv is 3.9.
#
# Deliberately outside the repo (like the hackathon sample data itself, see
# samples.py) -- a real ~1.5-hour CDS27601 run once finished cleanly (log
# ended with the real "Completed annotating all regions" line) but left both
# output CSVs at 0 bytes, matching this same repo directory's independently
# confirmed exposure to Google Drive folder sync (a stray `.tmp.driveupload/`
# was found at the repo root during an unrelated `git push`). A long-lived
# process writing steadily to a file for that long, inside a folder Drive is
# also trying to sync, is exactly the shape of thing that goes wrong.
RESULTS_ROOT = Path("/Users/joedevaney/paraphase_results")

# The 14 regions interpret_paraphase.py annotates by default when -r is
# omitted -- passed explicitly here anyway so the behavior doesn't silently
# change if a future vendored version's default list does.
MEDICAL_REGIONS = [
    "rccx", "smn1", "pms2", "strc", "cfc1", "ikbkg", "ncf1",
    "neb", "f8", "hba", "TNXB", "OTOA", "GBA", "CFHclust",
]

# One subprocess.Popen per sample currently running, so a second start_run()
# call while one is in flight attaches to the same run instead of launching a
# duplicate. Cleared once the process exits (see status()).
_RUNNING: dict[str, subprocess.Popen] = {}


def results_dir(sample: str) -> Path:
    return RESULTS_ROOT / sample


def log_path(sample: str) -> Path:
    return results_dir(sample) / "run.log"


def summary_csv(sample: str) -> Path:
    return results_dir(sample) / f"{sample}_summary.csv"


def detailed_csv(sample: str) -> Path:
    return results_dir(sample) / f"{sample}_detailed.csv"


def preflight(sample: str) -> Optional[str]:
    """None if a run can proceed; otherwise the reason it can't."""
    if sample not in samples.SAMPLES:
        return f"unknown sample {sample!r}"
    if not VENDOR_SCRIPT.exists():
        return (
            "interpret_paraphase.py not set up on this machine -- run "
            "`python backend/scripts/setup_paraphase_annotate.py` once, then retry."
        )
    if not samples.json_path(sample).exists():
        return f"Paraphase JSON not found: {samples.json_path(sample)}"
    vcfs_dir = samples.vcfs_dir(sample)
    if not vcfs_dir.exists():
        return (
            f"Paraphase per-gene VCFs not found at {vcfs_dir} -- interpret_paraphase.py "
            "requires a `<sample>_paraphase_vcfs/` directory (Paraphase's own output "
            "layout) alongside the JSON."
        )
    return None


_COMPLETION_MARKER = "Completed annotating all regions"


def status(sample: str) -> dict:
    """One of: blocked (preflight failed), running, done, failed, not_run."""
    proc = _RUNNING.get(sample)
    if proc is not None and proc.poll() is None:
        return {"state": "running"}
    if proc is not None:
        del _RUNNING[sample]

    log = log_path(sample)
    if log.exists():
        log_text = log.read_text(errors="replace")
        if _COMPLETION_MARKER in log_text:
            # Belt-and-suspenders: a real run has left the log's own
            # completion marker with a 0-byte (or missing) CSV before --
            # not a bug in interpret_paraphase.py, which writes each file's
            # header the moment it opens it, but external interference with
            # a long-lived write (see RESULTS_ROOT's own comment). Don't
            # report "done" over genuinely empty results.
            summary, detailed = summary_csv(sample), detailed_csv(sample)
            if not summary.exists() or summary.stat().st_size == 0 or not detailed.exists() or detailed.stat().st_size == 0:
                return {
                    "state": "failed",
                    "log_tail": (
                        "Run log reports completion, but its output file(s) are empty or "
                        f"missing ({summary.name}: {summary.stat().st_size if summary.exists() else 'missing'} bytes, "
                        f"{detailed.name}: {detailed.stat().st_size if detailed.exists() else 'missing'} bytes). "
                        "Likely external interference with the write (not a script bug) -- retry the run.\n\n"
                    ) + log_text[-2000:],
                }
            return {"state": "done"}
        # A run was attempted but never reached its own completion marker --
        # failed, regardless of whether it left partial CSVs behind.
        # interpret_paraphase.py's main() opens both output files and writes
        # their headers *before* processing the first region (see its own
        # source), so a CSV existing on disk doesn't mean the run succeeded.
        return {"state": "failed", "log_tail": log_text[-4000:]}

    reason = preflight(sample)
    if reason is not None:
        return {"state": "blocked", "reason": reason}
    return {"state": "not_run"}


def start_run(sample: str) -> dict:
    reason = preflight(sample)
    if reason is not None:
        return {"state": "blocked", "reason": reason}
    if sample in _RUNNING and _RUNNING[sample].poll() is None:
        return {"state": "running"}

    out_dir = results_dir(sample)
    out_dir.mkdir(parents=True, exist_ok=True)
    for stale in (summary_csv(sample), detailed_csv(sample)):
        stale.unlink(missing_ok=True)

    log_file = open(log_path(sample), "w")
    proc = subprocess.Popen(
        [
            sys.executable, str(VENDOR_SCRIPT),
            "--input", str(samples.json_path(sample)),
            "--out", str(out_dir),
            "--region", ",".join(MEDICAL_REGIONS),
            "--config", str(VENDOR_CONFIG),
        ],
        stdout=log_file, stderr=subprocess.STDOUT,
        cwd=BACKEND_ROOT,
    )
    _RUNNING[sample] = proc
    return {"state": "running"}
