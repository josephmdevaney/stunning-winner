"""Sample registry for the Paraphase Annotation suite.

Each entry needs a real Paraphase JSON, named `<key>.paraphase.json` per
interpret_paraphase.py's own sample-ID convention (see vendor/interpret_paraphase.py's
`main()`: `sample_id = os.path.basename(json_file).split(".paraphase.json")[0]`).
The per-gene VCFs interpret_paraphase.py annotates are expected as a sibling
directory `<key>_paraphase_vcfs/` next to the JSON -- Paraphase's own output
layout, not something this project controls.
"""

from __future__ import annotations

from pathlib import Path

# key -> path to that sample's <key>.paraphase.json. Real per-participant
# hackathon data, kept outside the repo (never committed) like this project's
# other real samples (see repeatscope/roh.py, mitolr/samples.py).
SAMPLES: dict[str, str] = {
    "CDS27601": "/Users/joedevaney/hackathon_cds2760x/CDS27601.paraphase.json",
    "CDS48001": "/Users/joedevaney/hackathon_cds48001/CDS48001.paraphase.json",
    "CDS57401": "/Users/joedevaney/hackathon_cds57401/CDS57401.paraphase.json",
}


def sample_keys() -> list[str]:
    return list(SAMPLES)


def json_path(sample: str) -> Path:
    return Path(SAMPLES[sample])


def vcfs_dir(sample: str) -> Path:
    """Paraphase's own fixed layout: `<sampleID>_paraphase_vcfs/` next to the JSON."""
    return json_path(sample).parent / f"{sample}_paraphase_vcfs"
