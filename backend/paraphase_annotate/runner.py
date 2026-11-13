"""Runs PacBio's own `interpret_paraphase.py` (vendored, see
scripts/setup_paraphase_annotate.py) against a sample's real Paraphase JSON +
per-gene VCFs, and tracks that run's status.

The script makes live Ensembl VEP + NCBI ClinVar API calls per variant (rate
limited to ~3/sec by the script itself), so a run can take anywhere from
seconds to a few minutes depending on variant count -- too slow to run inline
in a request handler. It's launched as a background subprocess instead; the
API layer (router.py) polls status() rather than blocking on it.

Which Paraphase produced a sample's JSON and which release of the script is
vendored are independent, so each run prepares its own input (see
`prepare_input`): the newest (v4.x) script reads region-specific fields from
`region_specific_info`, which older (v3.x) Paraphase output doesn't have, and
`calls.to_v4` bridges that. What was used for a run is recorded next to its
results (`run_meta.json`) so a result can always be traced to a tool version.
"""

from __future__ import annotations

import datetime
import json
import os
import shutil
import subprocess
import sys
from functools import lru_cache
from pathlib import Path
from typing import Optional

from . import calls, launch, samples

BACKEND_ROOT = Path(__file__).resolve().parents[1]
VENDOR_DIR = Path(__file__).resolve().parent / "vendor"
VENDOR_SCRIPT = VENDOR_DIR / "interpret_paraphase.py"
# Runs the script with a timeout on its HTTP requests (see launch.py); the script sets none of its own.
LAUNCHER = Path(__file__).resolve().parent / "launch.py"
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

VENDOR_VERSION_FILE = VENDOR_DIR / "VERSION.json"

# interpret_paraphase.py v4.x uses PEP 701 f-strings (Python 3.12+); v3.x runs on
# anything this project supports.
_V4_MIN_PYTHON = (3, 12)
_INTERPRETER_ENV = "PARAPHASE_ANNOTATE_PYTHON"


def vendor_info() -> Optional[dict]:
    """What setup_paraphase_annotate.py fetched (tag, when, checksums), or None
    when no script is vendored. A script vendored before VERSION.json existed was
    the v3.5.0 this project originally pinned."""
    if not VENDOR_SCRIPT.exists():
        return None
    try:
        return json.loads(VENDOR_VERSION_FILE.read_text())
    except (OSError, ValueError):
        return {"paraphase_tag": "v3.5.0", "fetched_at": None, "source": None, "sha256": {}}


def script_major(info: Optional[dict] = None) -> int:
    info = info if info is not None else vendor_info()
    tag = (info or {}).get("paraphase_tag") or "v3.5.0"
    try:
        return int(tag.lstrip("v").split(".")[0])
    except ValueError:
        return 3


def _probe(python: str, min_version: tuple[int, int]) -> bool:
    code = f"import sys, requests, yaml; sys.exit(0 if sys.version_info >= {min_version!r} else 1)"
    try:
        return subprocess.run([python, "-c", code], capture_output=True, timeout=30).returncode == 0
    except (OSError, subprocess.SubprocessError):
        return False


@lru_cache(maxsize=4)
def find_interpreter(min_version: tuple[int, int] = _V4_MIN_PYTHON) -> Optional[str]:
    """A Python >= min_version that has `requests` and `pyyaml` (what the script
    imports). The API server's own interpreter is tried first, then
    $PARAPHASE_ANNOTATE_PYTHON, then python3.14/3.13/3.12 on PATH and in the
    usual Homebrew/conda/framework locations."""
    candidates: list[Optional[str]] = [os.environ.get(_INTERPRETER_ENV), sys.executable]
    for name in ("python3.14", "python3.13", "python3.12", "python3"):
        candidates.append(shutil.which(name))
    for d in ("/opt/homebrew/bin", "/usr/local/bin", str(Path.home() / "miniforge3" / "bin"), str(Path.home() / "miniconda3" / "bin")):
        for name in ("python3.14", "python3.13", "python3.12"):
            candidates.append(str(Path(d) / name))
    seen: set[str] = set()
    for c in candidates:
        if not c or c in seen or not Path(c).exists():
            continue
        seen.add(c)
        if _probe(c, min_version):
            return c
    return None


def interpreter_for_script() -> Optional[str]:
    return find_interpreter(_V4_MIN_PYTHON if script_major() >= 4 else (3, 9))


# One subprocess.Popen per sample currently running, so a second start_run()
# call while one is in flight attaches to the same run instead of launching a
# duplicate. Cleared once the process exits (see status()).
_RUNNING: dict[str, subprocess.Popen] = {}


def results_dir(sample: str) -> Path:
    return RESULTS_ROOT / sample


def log_path(sample: str) -> Path:
    return results_dir(sample) / "run.log"


def input_dir(sample: str) -> Path:
    return results_dir(sample) / "input"


def run_meta_path(sample: str) -> Path:
    return results_dir(sample) / "run_meta.json"


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
    if interpreter_for_script() is None:
        need = "3.12+" if script_major() >= 4 else "3.9+"
        return (
            f"The vendored interpret_paraphase.py ({(vendor_info() or {}).get('paraphase_tag')}) needs Python {need} "
            "with `requests` and `pyyaml` installed, and none was found. Install them for a suitable "
            f"Python or point ${_INTERPRETER_ENV} at one, then retry."
        )
    try:
        source = calls.load_calls(samples.json_path(sample))
    except (OSError, ValueError) as e:
        return f"Paraphase JSON could not be read: {e}"
    if calls.schema(source) == calls.SCHEMA_V4 and script_major() < 4:
        return (
            "This Paraphase JSON is from Paraphase v4+, but the vendored interpret_paraphase.py is v3.x, which "
            "can't read it. Run `python backend/scripts/setup_paraphase_annotate.py` to vendor v4.1.0."
        )
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
    archive_previous(sample)

    input_json = prepare_input(sample)
    python = interpreter_for_script()
    assert python is not None  # preflight() already ensured it
    log_file = open(log_path(sample), "w")
    proc = subprocess.Popen(
        [
            python, str(LAUNCHER), str(VENDOR_SCRIPT),
            "--input", str(input_json),
            "--out", str(out_dir),
            "--region", ",".join(MEDICAL_REGIONS),
            "--config", str(VENDOR_CONFIG),
        ],
        stdout=log_file, stderr=subprocess.STDOUT,
        cwd=BACKEND_ROOT,
    )
    _RUNNING[sample] = proc
    return {"state": "running"}


def archive_previous(sample: str) -> Optional[Path]:
    """Move a previous run's results aside (to `archive/<UTC timestamp>/` in the
    same results directory) instead of deleting them. A real run can take over an
    hour of rate-limited lookups, so a re-run that then fails must not cost the
    results it was meant to improve on. Returns the archive directory, or None
    when there was nothing to keep."""
    previous = [p for p in (summary_csv(sample), detailed_csv(sample), log_path(sample), run_meta_path(sample)) if p.exists()]
    if not previous:
        return None
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    dest = results_dir(sample) / "archive" / stamp
    dest.mkdir(parents=True, exist_ok=True)
    for p in previous:
        shutil.move(str(p), dest / p.name)
    return dest


def prepare_input(sample: str) -> Path:
    """Write the JSON the vendored script will read, and record how this run was
    set up. The script finds a sample's VCFs by looking for
    `<sampleID>_paraphase_vcfs/` next to its input JSON, so the (normalized)
    JSON lives in a work directory beside a symlink to the real VCF directory;
    nothing under the sample's own directory is written to or modified."""
    source_path = samples.json_path(sample)
    source = calls.load_calls(source_path)
    source_schema = calls.schema(source)
    script_schema = calls.SCHEMA_V4 if script_major() >= 4 else calls.SCHEMA_V3
    data = calls.to_v4(source) if script_schema == calls.SCHEMA_V4 else source

    work = input_dir(sample)
    work.mkdir(parents=True, exist_ok=True)
    target = work / f"{sample}.paraphase.json"
    target.write_text(json.dumps(data))
    link = work / f"{sample}_paraphase_vcfs"
    if link.is_symlink() or link.exists():
        link.unlink()
    link.symlink_to(samples.vcfs_dir(sample))

    info = vendor_info() or {}
    present = {r.lower() for r in source}
    meta = {
        "sample": sample,
        "started_at": _now(),
        "script_tag": info.get("paraphase_tag"),
        "script_fetched_at": info.get("fetched_at"),
        "python": interpreter_for_script(),
        "input_schema": source_schema,
        "script_schema": script_schema,
        "adapted_json": source_schema != script_schema,
        "data_paraphase_version": calls.vcf_provenance(samples.vcfs_dir(sample))["paraphase_version"],
        "request_timeouts": {"connect": launch.CONNECT_TIMEOUT, "read": launch.READ_TIMEOUT},
        "regions_requested": MEDICAL_REGIONS,
        # The script skips (with a warning) any region absent from the JSON.
        "regions_missing_from_json": [r for r in MEDICAL_REGIONS if r.lower() not in present],
    }
    run_meta_path(sample).write_text(json.dumps(meta, indent=2) + "\n")
    return target


def _now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def run_meta(sample: str) -> Optional[dict]:
    """How this sample's most recent run was set up (script release, Python,
    whether its JSON was adapted, which Paraphase made the data), plus when it
    finished once it has. None if it has never been run."""
    try:
        meta = json.loads(run_meta_path(sample).read_text())
    except (OSError, ValueError):
        return None
    summary = summary_csv(sample)
    if status(sample).get("state") == "done" and summary.exists():
        meta["finished_at"] = datetime.datetime.fromtimestamp(
            summary.stat().st_mtime, datetime.timezone.utc
        ).isoformat(timespec="seconds")
    return meta
