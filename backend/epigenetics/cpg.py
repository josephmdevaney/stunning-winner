"""Reading PacBio pb-CpG-tools per-CpG methylation beds (`aligned_bam_to_cpg_scores`).

A genome's pileup is three bgzipped, tabix-indexed beds -- `combined`, `hap1`, `hap2` -- one row per
CpG: chrom, start, end, mod_score (percent methylated), type, coverage, est_mod_count,
est_unmod_count, discretized_mod_score. Region queries go through the `tabix` binary (a 30x genome's
combined bed is ~300 MB, so reading it whole is not an option) and cost well under a second.

The first lines of each file record how it was made (pb-CpG-tools version, pileup mode, modsites
mode, minimum coverage); `header()` returns them so a result can always say what it was computed from.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


class TabixError(RuntimeError):
    pass


@dataclass(frozen=True)
class CpG:
    pos: int            # 0-based start, as in the bed
    mod_score: float    # percent methylated, 0-100
    coverage: int
    est_mod: int
    est_unmod: int


def tabix_path() -> str:
    found = shutil.which("tabix")
    if not found:
        raise TabixError(
            "`tabix` (from htslib) was not found on PATH. Install it (for example `brew install htslib`) "
            "to read pb-CpG-tools beds."
        )
    return found


def _run(args: list[str]) -> str:
    try:
        proc = subprocess.run([tabix_path(), *args], capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired as e:
        raise TabixError(f"tabix timed out: {' '.join(args)}") from e
    if proc.returncode != 0:
        raise TabixError(proc.stderr.strip() or f"tabix exited {proc.returncode}")
    return proc.stdout


def header(bed: Path) -> dict:
    """pb-CpG-tools provenance from the `##` header lines (all optional)."""
    meta: dict = {"pb_cpg_tools_version": None, "pileup_mode": None, "modsites_mode": None, "min_coverage": None}
    for line in _run(["-H", str(bed), "chrM:1-1"]).splitlines():
        if not line.startswith("##"):
            continue
        if line.startswith("##pb-cpg-tools-version="):
            meta["pb_cpg_tools_version"] = line.split("=", 1)[1].strip()
        elif line.startswith("##pileup-mode="):
            meta["pileup_mode"] = line.split("=", 1)[1].strip()
        elif line.startswith("##modsites-mode="):
            meta["modsites_mode"] = line.split("=", 1)[1].strip()
        elif line.startswith("##cmdline="):
            m = re.search(r"--min-coverage\s+(\d+)", line)
            meta["min_coverage"] = int(m.group(1)) if m else None
    return meta


def fetch(bed: Path, chrom: str, start: int, end: int) -> list[CpG]:
    """CpGs overlapping [start, end) on `chrom` (0-based half-open, bed coordinates)."""
    if end <= start:
        return []
    out = _run([str(bed), f"{chrom}:{start + 1}-{end}"])
    rows = []
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 8 or line.startswith("#"):
            continue
        rows.append(CpG(int(f[1]), float(f[3]), int(f[5]), int(float(f[6])), int(float(f[7]))))
    return rows


def fetch_sites(bed: Path, sites: list[tuple[str, int]]) -> dict[tuple[str, int], CpG]:
    """The CpG at each (chromosome, 0-based start) in `sites`, read in one tabix call (`tabix -R`), for the few
    thousand positions of a signature. Chromosome names may be given without "chr"; they are matched to the
    file's own naming. Sites with no row in the pileup are absent from the result."""
    names = contigs(bed) or []
    prefix = "chr" if any(n.startswith("chr") for n in names) else ""
    fd, regions = tempfile.mkstemp(suffix=".bed")      # tabix reads a -R file as BED only from the extension
    try:
        with os.fdopen(fd, "w") as f:
            for chrom, start in sorted(set(sites), key=lambda s: (names.index(prefix + s[0]) if prefix + s[0] in names else len(names), s[1])):
                f.write(f"{prefix}{chrom}\t{start}\t{start + 1}\n")
        out = _run(["-R", regions, str(bed)])
    finally:
        os.unlink(regions)
    found: dict[tuple[str, int], CpG] = {}
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 8 or line.startswith("#"):
            continue
        chrom = f[0][3:] if prefix and f[0].startswith("chr") else f[0]
        found[(chrom, int(f[1]))] = CpG(int(f[1]), float(f[3]), int(f[5]), int(float(f[6])), int(float(f[7])))
    return found


def contigs(bed: Path) -> Optional[list[str]]:
    """Contig names in the tabix index (also a cheap check that the index is present and readable)."""
    try:
        return _run(["-l", str(bed)]).split()
    except TabixError:
        return None
