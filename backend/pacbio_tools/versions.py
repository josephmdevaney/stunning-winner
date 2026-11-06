"""Version strings from tool output ("5.0.0", "sawfish 2.0.3", "0.2.7-8ee9b5b",
"v3.5.0") compared to a release tag."""

from __future__ import annotations

import re
from typing import Optional

_VERSION = re.compile(r"(\d+)\.(\d+)(?:\.(\d+))?")


def parse(text: Optional[str]) -> Optional[tuple[int, int, int]]:
    """(major, minor, patch) from the first x.y[.z] in `text`, or None."""
    if not text:
        return None
    m = _VERSION.search(str(text))
    if not m:
        return None
    return int(m.group(1)), int(m.group(2)), int(m.group(3) or 0)


def clean(text: Optional[str]) -> Optional[str]:
    """The bare x.y.z (without a leading tool name or "v", or a trailing build id)."""
    v = parse(text)
    return ".".join(str(n) for n in v) if v else None


def gap(observed: Optional[str], latest: Optional[str]) -> Optional[str]:
    """How far `observed` is behind `latest`: "major", "minor", "patch", "none"
    when equal, "ahead" when newer than the latest release we know of, or None
    when either can't be parsed."""
    o, l = parse(observed), parse(latest)
    if o is None or l is None:
        return None
    if o == l:
        return "none"
    if o > l:
        return "ahead"
    if o[0] != l[0]:
        return "major"
    if o[1] != l[1]:
        return "minor"
    return "patch"
