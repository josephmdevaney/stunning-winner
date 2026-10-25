"""The tool-versions report: registry + newest releases + observed versions."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from . import observed, registry, releases, versions

# How serious a gap is, worst first.
_RANK = {"major": 3, "minor": 2, "patch": 1, "none": 0, "ahead": 0}


def _observation_status(raw: Optional[str], latest_version: Optional[str]) -> dict:
    gap = versions.gap(raw, latest_version) if raw else None
    return {"version": versions.clean(raw), "gap": gap}


def build(observations: dict[str, list[dict]], latest: dict) -> dict:
    """`observations` maps a tool key to `observed.group()` output (plus any
    hand-labelled entries); `latest` is `releases.latest_releases()`."""
    tools = []
    for tool in registry.TOOLS:
        rel = latest["releases"].get(tool.key)
        latest_version = (rel or {}).get("version")
        obs = []
        for o in observations.get(tool.key, []):
            s = _observation_status(o.get("raw"), latest_version)
            obs.append({**o, **s})

        gaps = [o["gap"] for o in obs if o["gap"] is not None]
        if not tool.used:
            status = "not_used"
        elif not obs:
            status = "not_recorded"
        elif not gaps:
            status = "unknown"
        else:
            status = "current" if all(g in ("none", "ahead") for g in gaps) else "behind"
        worst = max(gaps, key=lambda g: _RANK[g]) if gaps else None

        tools.append({
            "key": tool.key, "name": tool.name, "repo": tool.repo, "url": f"https://github.com/{tool.repo}",
            "role": tool.role, "used": tool.used, "suite_use": tool.suite_use,
            "latest": rel, "observed": obs, "status": status, "worst_gap": worst,
        })
    used = [t for t in tools if t["used"]]
    return {
        "source": latest["source"],
        "fetched_at": latest["fetched_at"],
        "errors": latest.get("errors", {}),
        "summary": {
            "used": len(used),
            "current": sum(t["status"] == "current" for t in used),
            "behind": sum(t["status"] == "behind" for t in used),
            "not_recorded": sum(t["status"] in ("not_recorded", "unknown") for t in used),
        },
        "tools": tools,
    }
