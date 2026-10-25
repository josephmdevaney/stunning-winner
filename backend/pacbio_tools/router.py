"""FastAPI routes for the PacBio tool-versions report, mounted on the shared app
under /api/tools.

`GET /api/tools/versions` compares the version of each PacBio tool that produced
the data this app has loaded with that tool's newest GitHub release. The
observed versions are read from the files themselves (VCF headers, recorded
caller strings); `?refresh=true` forces a live GitHub check instead of serving
the cached/snapshot releases.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from fastapi import APIRouter

from . import observed, releases, report

router = APIRouter(prefix="/api/tools", tags=["pacbio-tools"])


def _candidates(base: Path, sample: str, folder_suffix: str, name_suffixes: tuple[str, ...]) -> list[Path]:
    """The flat and per-sample-subdirectory layouts PacBio's pipelines write (the
    same two the Carrier loader accepts)."""
    flat = [base / f"{sample}{s}" for s in name_suffixes]
    nested = [base / f"{sample}{folder_suffix}" / f"{sample}{s}" for s in name_suffixes]
    return flat + nested


def collect() -> dict[str, list[dict]]:
    """Observed versions per tool key, from whatever data the running app has."""
    # Imported here, not at module load: the app imports this router to mount it.
    from mitolr import samples as mito_samples
    from paraphase_annotate import calls as paraphase_calls, runner as paraphase_runner, samples as paraphase_samples
    from repeatscope.api import main as api_main
    from repeatscope.trgt import trgt_version

    trgt, sawfish, paraphase = [], [], []
    for mode in api_main._modes.values():
        label = mode["label"]
        for sample, vcf in mode["sample_paths"].items():
            trgt.append((label, trgt_version(vcf), sample))
        base = mode.get("paraphase_dir")
        if not base:
            continue
        base = Path(base)
        for sample in mode["sample_paths"]:
            sv = next((p for p in _candidates(base, sample, "_sawfish", (".sv.vcf", ".sv.vcf.gz")) if p.exists()), None)
            if sv is not None:
                sawfish.append((label, observed.sawfish_version(sv), sample))
            vcfs = base / f"{sample}_paraphase" / f"{sample}_paraphase_vcfs"
            paraphase.append((label, paraphase_calls.vcf_provenance(vcfs)["paraphase_version"], sample))

    for sample in paraphase_samples.sample_keys():
        version = paraphase_calls.vcf_provenance(paraphase_samples.vcfs_dir(sample))["paraphase_version"]
        paraphase.append(("Paraphase Annotation samples", version, sample))

    mito = [
        ("MitoLR samples", observed.mitorsaw_version(s.get("caller")), key)
        for key, s in mito_samples.SAMPLES.items() if s.get("source") == "vcf"
    ]

    out: dict[str, list[dict]] = {
        "trgt": observed.group(trgt),
        "sawfish": observed.group(sawfish),
        "paraphase": observed.group(paraphase),
        "mitorsaw": observed.group(mito),
    }
    script = paraphase_runner.vendor_info()
    if script:
        out["paraphase"].append({
            "label": "Vendored annotation script (interpret_paraphase.py)", "raw": script.get("paraphase_tag"),
            "count": 1, "examples": [],
        })
    carrier = api_main._modes.get("carrier")
    if carrier:
        out["ptcp"] = [{
            "label": "Carrier panel (documented for this cohort, not read from a file)", "raw": "3.2.1",
            "count": len(carrier["sample_paths"]), "examples": [],
        }]
    return out


@router.get("/versions")
def get_versions(refresh: bool = False):
    return report.build(collect(), releases.latest_releases(refresh=refresh))
