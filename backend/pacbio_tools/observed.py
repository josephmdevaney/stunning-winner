"""Which version of each tool produced the data this suite loads.

Every reader here takes a path or a string and returns the version text the tool
itself wrote (VCF header lines, a recorded caller string), never a guess. Where
the suite doesn't read a file that records a tool's version, the tool simply has
no observation (shown as "not recorded"), rather than an assumed one.
"""

from __future__ import annotations

import gzip
import re
from collections import defaultdict
from pathlib import Path
from typing import Callable, Iterable, Optional

_SOURCE_LINE = re.compile(r'^##source="?([^"\n]+)"?\s*$')


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if str(path).endswith(".gz") else open(path, encoding="utf-8", errors="replace")


def vcf_header_value(path: Path, pattern: re.Pattern) -> Optional[str]:
    """Group 1 of the first `##` header line matching `pattern`, or None."""
    try:
        with _open_text(path) as f:
            for line in f:
                if not line.startswith("##"):
                    break
                m = pattern.match(line.rstrip("\n"))
                if m:
                    return m.group(1).strip()
    except OSError:
        return None
    return None


def sawfish_version(vcf: Path) -> Optional[str]:
    """`##source="sawfish 2.0.3"` -> "sawfish 2.0.3" (None for another tool's VCF)."""
    src = vcf_header_value(vcf, _SOURCE_LINE)
    return src if src and src.lower().startswith("sawfish") else None


def mitorsaw_version(caller: Optional[str]) -> Optional[str]:
    """A recorded caller string such as "mitorsaw 0.2.7-8ee9b5b" -> "0.2.7-8ee9b5b"."""
    m = re.match(r"^\s*mitorsaw\s+(\S+)", caller or "", flags=re.IGNORECASE)
    return m.group(1) if m else None


Reading = tuple[str, Optional[str], str]  # (group label, raw version text or None, sample id)


def group(readings: Iterable[Reading]) -> list[dict]:
    """One observation per (group label, version): how many samples carry it and a
    few example IDs. Samples whose version couldn't be read are counted separately
    so a partly-recorded cohort doesn't look fully recorded."""
    buckets: dict[tuple[str, Optional[str]], list[str]] = defaultdict(list)
    for label, raw, sample in readings:
        buckets[(label, raw)].append(sample)
    out = []
    for (label, raw), samples in sorted(buckets.items(), key=lambda kv: (kv[0][0], kv[0][1] or "")):
        out.append({
            "label": label,
            "raw": raw,
            "count": len(samples),
            "examples": sorted(samples)[:3],
        })
    return out


def read_all(label: str, items: dict[str, Path], reader: Callable[[Path], Optional[str]]) -> list[Reading]:
    return [(label, reader(path), sample) for sample, path in items.items()]
