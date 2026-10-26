"""Sample registry for the Paraphase Annotation suite.

Each entry needs a real Paraphase JSON, named `<key>.paraphase.json` per
interpret_paraphase.py's own sample-ID convention (see vendor/interpret_paraphase.py's
`main()`: `sample_id = os.path.basename(json_file).split(".paraphase.json")[0]`).
The per-gene VCFs interpret_paraphase.py annotates are expected as a sibling
directory `<key>_paraphase_vcfs/` next to the JSON -- Paraphase's own output
layout, not something this project controls.
"""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path

# This machine's own samples: {"samples": {"<key>": "/path/to/<key>.paraphase.json"}}, in a gitignored file so
# machine-specific paths stay out of the repository (see paraphase.local.json.example). Read when this module is
# imported, so the server needs a restart to pick up a change.
LOCAL_CONFIG = Path(os.environ.get("PARAPHASE_SAMPLES_CONFIG", Path(__file__).resolve().parents[1] / "config" / "paraphase.local.json"))


def _local_samples(config: Path = LOCAL_CONFIG) -> dict[str, str]:
    """The samples listed in the local config whose JSON exists and is named `<key>.paraphase.json`, which is how
    interpret_paraphase.py derives the sample ID. Anything else is skipped with a warning rather than breaking the app."""
    try:
        raw = json.loads(config.read_text())
    except (OSError, ValueError):
        return {}
    found: dict[str, str] = {}
    for key, path in (raw.get("samples") or {}).items():
        path = Path(path).expanduser()
        if path.name != f"{key}.paraphase.json" or not path.is_file():
            logging.warning("paraphase sample %r in %s skipped: expected an existing file named %s.paraphase.json, got %s", key, config, key, path)
            continue
        found[key] = str(path)
    return found


# key -> path to that sample's <key>.paraphase.json. Real per-participant
# hackathon data, kept outside the repo (never committed) like this project's
# other real samples (see roh/analysis.py, mitolr/samples.py).
SAMPLES: dict[str, str] = {
    "CDS27601": "/Users/joedevaney/hackathon_cds2760x/CDS27601.paraphase.json",
    "CDS48001": "/Users/joedevaney/hackathon_cds48001/CDS48001.paraphase.json",
    "CDS57401": "/Users/joedevaney/hackathon_cds57401/CDS57401.paraphase.json",
    **_local_samples(),
}


def sample_keys() -> list[str]:
    return list(SAMPLES)


def json_path(sample: str) -> Path:
    return Path(SAMPLES[sample])


def vcfs_dir(sample: str) -> Path:
    """Paraphase's own fixed layout: `<sampleID>_paraphase_vcfs/` next to the JSON."""
    return json_path(sample).parent / f"{sample}_paraphase_vcfs"
