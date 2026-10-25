"""The newest release of each registered PacBio tool.

Three sources, freshest first: a live GitHub query, an on-disk cache of the last
live query (so a page load doesn't hit GitHub every time), and a snapshot
committed with the code (so an offline machine still has something to compare
against, with its date shown). `scripts/refresh_pacbio_releases.py` rewrites the
snapshot; the app only ever writes the cache.
"""

from __future__ import annotations

import datetime
import json
import os
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Optional

from . import registry, versions

SNAPSHOT_PATH = Path(__file__).with_name("releases_snapshot.json")
CACHE_PATH = Path(os.environ.get("PACBIO_TOOLS_CACHE", Path.home() / ".cache" / "stunning-winner" / "pacbio_releases.json"))
CACHE_TTL_SECONDS = 6 * 3600
_REQUEST_TIMEOUT = 6


def _now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def _read(path: Path) -> Optional[dict]:
    try:
        data = json.loads(path.read_text())
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) and isinstance(data.get("releases"), dict) else None


def _get_json(url: str):
    headers = {"Accept": "application/vnd.github+json", "User-Agent": "stunning-winner-tool-versions"}
    token = os.environ.get("GITHUB_TOKEN")
    if token:
        headers["Authorization"] = f"Bearer {token}"
    with urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=_REQUEST_TIMEOUT) as resp:
        return json.load(resp)


def fetch_one(repo: str) -> dict:
    """{"tag", "version", "published", "name", "url"} for a repo's latest release
    (falling back to its newest tag when it publishes no releases). Raises on a
    network or API failure."""
    try:
        rel = _get_json(f"https://api.github.com/repos/{repo}/releases/latest")
        tag = rel["tag_name"]
        return {
            "tag": tag, "version": versions.clean(tag), "published": (rel.get("published_at") or "")[:10] or None,
            "name": rel.get("name") or tag, "url": rel.get("html_url"),
        }
    except urllib.error.HTTPError as e:
        if e.code != 404:
            raise
    tags = _get_json(f"https://api.github.com/repos/{repo}/tags?per_page=1")
    if not tags:
        raise LookupError(f"{repo} has no releases or tags")
    tag = tags[0]["name"]
    return {"tag": tag, "version": versions.clean(tag), "published": None, "name": tag,
            "url": f"https://github.com/{repo}/releases/tag/{tag}"}


def fetch_live() -> tuple[dict[str, dict], dict[str, str]]:
    """(releases by tool key, error text by tool key) from GitHub, in parallel."""
    def one(tool):
        try:
            return tool.key, fetch_one(tool.repo), None
        except Exception as e:  # network down, rate limited, repo renamed: report per tool, keep going
            return tool.key, None, f"{type(e).__name__}: {e}"

    with ThreadPoolExecutor(max_workers=8) as pool:
        results = list(pool.map(one, registry.TOOLS))
    return ({k: r for k, r, _ in results if r}, {k: err for k, _, err in results if err})


def _age_seconds(iso: Optional[str]) -> float:
    try:
        return (datetime.datetime.now(datetime.timezone.utc) - datetime.datetime.fromisoformat(iso)).total_seconds()
    except (TypeError, ValueError):
        return float("inf")


def latest_releases(refresh: bool = False) -> dict:
    """{"source": "live" | "cache" | "snapshot", "fetched_at", "releases": {key: {...}}, "errors": {key: str}}.

    A live query runs when `refresh` is set or the cache is older than the TTL.
    Tools that fail to fetch keep their previous value (from the cache, else the
    snapshot) and are listed in `errors`, so one bad repo never blanks the page."""
    snapshot = _read(SNAPSHOT_PATH) or {"fetched_at": None, "releases": {}}
    cache = _read(CACHE_PATH)
    base = cache if cache and _age_seconds(cache.get("fetched_at")) != float("inf") else snapshot
    base_source = "cache" if base is cache else "snapshot"

    if not refresh and cache and _age_seconds(cache.get("fetched_at")) < CACHE_TTL_SECONDS:
        return {"source": "cache", "fetched_at": cache["fetched_at"], "releases": cache["releases"], "errors": {}}

    live, errors = fetch_live()
    if not live:
        return {"source": base_source, "fetched_at": base.get("fetched_at"), "releases": base["releases"], "errors": errors}

    merged = {**snapshot["releases"], **(cache or {}).get("releases", {}), **live}
    fetched_at = _now()
    try:
        CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
        CACHE_PATH.write_text(json.dumps({"fetched_at": fetched_at, "releases": merged}, indent=2) + "\n")
    except OSError:
        pass  # a read-only home directory only costs the caching
    return {"source": "live", "fetched_at": fetched_at, "releases": merged, "errors": errors}
