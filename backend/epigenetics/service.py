"""Computing an imprinting overview for a registered sample.

For every DMR in the panel: fetch the CpGs of its analysis window from each track, pool them, classify. The
panel's 48 DMRs cost about 150 tabix queries. The CpGs read are cached per (sample, window) and the finished
overview per (sample, window, thresholds), so changing a threshold re-classifies without re-reading the beds. Both
caches are keyed on the bed files' modification times, so a regenerated pileup is re-read.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Optional

from . import cpg
from .analysis import STATUS_LOW_COVERAGE, STATUS_NO_DATA, Thresholds, region_result
from .disorders import CAVEATS, evaluate_disorders, evaluate_mlid
from .samples import Sample

DATA = Path(__file__).resolve().parent / "data" / "imprinted_dmrs.json"


@lru_cache(maxsize=1)
def panel() -> dict:
    return json.loads(DATA.read_text())


def window_for(region: dict, window: str) -> tuple[str, int, int, str]:
    """(chrom, start, end, which) -- the tighter CpG-island core when asked for and available, else the DMR."""
    if window == "core" and region.get("core"):
        c = region["core"]
        return c["chrom"], c["start"], c["end"], "core"
    return region["chrom"], region["start"], region["end"], "dmr"


def _mtime_key(sample: Sample) -> tuple:
    return tuple((t, str(p), p.stat().st_mtime_ns) for t, p in sorted(sample.paths.items()))


_cache: dict[tuple, dict] = {}
_cpg_cache: dict[tuple, dict[str, dict]] = {}


def _window_cpgs(sample: Sample, window: str) -> dict[str, dict]:
    """{region id: {track: [CpG]}} for the analysis window of every panel DMR."""
    key = (sample.key, window, _mtime_key(sample))
    if key not in _cpg_cache:
        out = {}
        for region in panel()["regions"]:
            chrom, start, end, _ = window_for(region, window)
            out[region["id"]] = {t: (cpg.fetch(sample.paths[t], chrom, start, end) if t in sample.paths else []) for t in ("combined", "hap1", "hap2")}
        _cpg_cache[key] = out
    return _cpg_cache[key]


def overview(sample: Sample, window: str = "core", th: Thresholds = Thresholds()) -> dict:
    key = (sample.key, window, tuple(sorted(th.as_dict().items())), _mtime_key(sample))
    if key in _cache:
        return _cache[key]

    results: dict[str, dict] = {}
    window_cpgs = _window_cpgs(sample, window)
    for region in panel()["regions"]:
        chrom, start, end, which = window_for(region, window)
        tracks = window_cpgs[region["id"]]
        r = region_result(region, tracks["combined"], tracks["hap1"], tracks["hap2"], th, window=which)
        r["interval"] = {"chrom": chrom, "start": start, "end": end}
        r["dmr"] = {"chrom": region["chrom"], "start": region["start"], "end": region["end"]}
        r["n_cpgs_published"] = region["n_cpgs_published"]
        r["imprint_origin"] = region["imprint_origin"]
        r["coordinate_notes"] = region.get("coordinate_notes")
        results[region["id"]] = r

    called = [r for r in results.values() if r["status"] not in (STATUS_LOW_COVERAGE, STATUS_NO_DATA)]
    counts: dict[str, int] = {}
    for r in results.values():
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    combined_bed = sample.paths["combined"]
    out = {
        "sample": sample.key, "label": sample.label, "note": sample.note,
        "window": window, "thresholds": th.as_dict(),
        "provenance": {**cpg.header(combined_bed), "tracks": sorted(sample.paths), "genome_build": panel()["source"]["genome_build"]},
        "panel": {k: panel()["source"][k] for k in ("title", "journal", "pmcid", "license")},
        "summary": {"regions": len(results), "called": len(called), "by_status": counts},
        "regions": list(results.values()),
        "disorders": evaluate_disorders(results),
        "mlid": evaluate_mlid(results),
        "caveats": CAVEATS,
    }
    _cache[key] = out
    return out


def region_track(sample: Sample, region_id: str, flank: int = 3000, window: str = "core") -> Optional[dict]:
    """CpG-level methylation across a region and `flank` bp either side, for plotting."""
    region = next((r for r in panel()["regions"] if r["id"] == region_id), None)
    if region is None:
        return None
    chrom, start, end, which = window_for(region, window)
    lo, hi = max(0, min(start, region["start"]) - flank), max(end, region["end"]) + flank
    tracks = {}
    for t, p in sample.paths.items():
        tracks[t] = [
            {"pos": c.pos, "methylation": c.mod_score, "coverage": c.coverage}
            for c in cpg.fetch(p, chrom, lo, hi)
        ]
    return {
        "id": region["id"], "name": region["name"], "chrom": chrom,
        "dmr": {"start": region["start"], "end": region["end"]},
        "analysis_window": {"start": start, "end": end, "which": which},
        "view": {"start": lo, "end": hi}, "tracks": tracks,
    }
