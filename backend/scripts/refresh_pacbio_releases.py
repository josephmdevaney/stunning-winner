#!/usr/bin/env python3
"""Rewrite pacbio_tools/releases_snapshot.json from GitHub: the offline fallback
the tool-versions page compares against when it can't reach GitHub.

    python backend/scripts/refresh_pacbio_releases.py
"""

from __future__ import annotations

import datetime
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pacbio_tools import releases  # noqa: E402


def main() -> int:
    live, errors = releases.fetch_live()
    if errors:
        for key, err in errors.items():
            print(f"  {key}: {err}", file=sys.stderr)
    if not live:
        print("No release could be fetched; snapshot left unchanged.", file=sys.stderr)
        return 1
    existing = json.loads(releases.SNAPSHOT_PATH.read_text())["releases"] if releases.SNAPSHOT_PATH.exists() else {}
    merged = {**existing, **live}  # a tool that failed this time keeps its previous entry
    releases.SNAPSHOT_PATH.write_text(json.dumps({
        "fetched_at": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "releases": merged,
    }, indent=2) + "\n")
    for key in sorted(merged):
        r = merged[key]
        print(f"{key:14s} {r['tag']:12s} {r.get('published') or ''}")
    return 0 if not errors else 2


if __name__ == "__main__":
    sys.exit(main())
