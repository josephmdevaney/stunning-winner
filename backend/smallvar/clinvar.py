"""ClinVar lookup for the small-variant suite: an exact-allele match against a reduced
table built from NCBI's ``variant_summary.txt.gz`` by
``scripts/real_data/build_smallvar_clinvar.py``.

The reduced table keeps GRCh38 small variants that ClinVar classifies pathogenic or likely
pathogenic, or as conflicting classifications of pathogenicity. A variant that is not in
it is *not* shown to be absent from ClinVar: it may be benign, uncertain, or unreviewed.
The classification shown is ClinVar's own assertion with its review status; this suite
does not classify variants itself.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Dict, NamedTuple, Optional, Tuple

COLUMNS = (
    "chrom", "pos", "ref", "alt", "gene", "classification", "group", "stars", "review",
    "variation_id", "name", "conditions", "last_evaluated", "n_submitters",
)

_STARS = (
    ("practice guideline", 4),
    ("reviewed by expert panel", 3),
    ("criteria provided, multiple submitters, no conflicts", 2),
    ("criteria provided, single submitter", 1),
    ("criteria provided, conflicting", 1),
)


def review_stars(review_status: str) -> int:
    s = review_status.lower()
    for prefix, stars in _STARS:
        if s.startswith(prefix):
            return stars
    return 0


def classification_group(classification: str) -> str:
    """plp | conflicting | other, from ClinVar's ClinicalSignificance text."""
    for part in classification.split(";"):
        if part.strip().startswith(("Pathogenic", "Likely pathogenic")):
            return "plp"
    if classification.startswith("Conflicting"):
        return "conflicting"
    return "other"


class ClinvarRecord(NamedTuple):
    gene: str
    classification: str
    group: str
    stars: int
    review: str
    variation_id: str
    name: str
    conditions: str
    last_evaluated: str
    n_submitters: int


class ClinvarTable:
    def __init__(self, records: Dict[Tuple[str, int, str, str], ClinvarRecord], release: str, path: str, stamp: str):
        self._r = records
        self.release = release
        self.path = path
        self.stamp = stamp  # changes when the file does, so cached indexes know to rebuild

    def __len__(self) -> int:
        return len(self._r)

    def lookup(self, chrom: str, pos: int, ref: str, alt: str) -> Optional[ClinvarRecord]:
        """``pos``/``ref``/``alt`` must already be trimmed (see normalize.variant_key)."""
        return self._r.get((chrom, pos, ref, alt))

    @classmethod
    def load(cls, path: str | Path) -> "ClinvarTable":
        p = Path(path)
        records: Dict[Tuple[str, int, str, str], ClinvarRecord] = {}
        release = "unknown"
        opener = gzip.open if p.suffix == ".gz" else open
        with opener(p, "rt", encoding="utf-8") as f:  # type: ignore[arg-type]
            for line in f:
                if line.startswith("#"):
                    if line.startswith("# release:"):
                        release = line.split(":", 1)[1].strip()
                    continue
                v = line.rstrip("\n").split("\t")
                if len(v) < len(COLUMNS) or v[0] == "chrom":
                    continue
                records[(v[0], int(v[1]), v[2], v[3])] = ClinvarRecord(
                    v[4], v[5], v[6], int(v[7]), v[8], v[9], v[10], v[11], v[12], int(v[13] or 0)
                )
        st = p.stat()
        return cls(records, release, str(p), f"{p.name}:{st.st_size}:{st.st_mtime_ns}")
