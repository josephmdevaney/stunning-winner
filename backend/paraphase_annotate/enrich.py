"""Turns interpret_paraphase.py's own `Pathogenic_variants_annotated` cells
into structured, individually-verified variant records.

interpret_paraphase.py already calls Ensembl VEP and NCBI ClinVar for every
candidate variant in the whole VCF (hundreds per region, rate-limited to
~3/sec -- the reason a run can take tens of minutes), but (a) it doesn't keep
VEP's `hgvsc`/`hgvsp` fields, and (b) the ClinVar description it keeps is
whatever ClinVar record it happened to find first, not distinguished from a
genuinely definitive Pathogenic/Likely pathogenic call (e.g. "Conflicting
classifications of pathogenicity" prints in the same undifferentiated way a
clean "Pathogenic" would, in its raw output). Re-running its own pipeline just
to capture more fields would mean re-paying that entire rate-limited cost.
Instead, this module re-queries VEP (with `hgvs=1`) and ClinVar's esummary
(by the ClinVar ID interpret_paraphase.py already found) only for the small
number of variants that already made it into a `Pathogenic_variants_annotated`
cell -- a handful per sample, not the whole candidate set -- so it's fast, and
doesn't touch interpret_paraphase.py's own source at all.
"""

from __future__ import annotations

import re
import time
import urllib.parse
from functools import lru_cache
from typing import Optional, TypedDict

import requests

# One variant=info entry per comma-separated segment, each info made of one or
# more `:`-joined parts (see interpret_paraphase.py's write_region_output:
# `f"{variant}={':'.join(info_parts)}"`, comma-joined across variants).
_ENTRY_RE = re.compile(r"^([^=]+)=(.+)$")
_CLINVAR_PART_RE = re.compile(r"^ClinVar-id(\d+);(.*);$")
_TRANSCRIPT_PART_RE = re.compile(r"^(ENST\d+);(.*);([A-Z_]+);$")

# ACMG's own definitive terms -- everything else (Conflicting classifications
# of pathogenicity, Uncertain significance, Benign, ...) is real ClinVar data
# too, just not a clean pathogenic call, so it must not read as one.
_DEFINITIVE_CLASSIFICATIONS = {"pathogenic", "likely pathogenic", "pathogenic/likely pathogenic"}

_MIN_REQUEST_INTERVAL_S = 0.34  # match interpret_paraphase.py's own rate limit
_last_request_at = 0.0


class VariantRecord(TypedDict):
    variant: str
    transcript: Optional[str]
    consequence_terms: list[str]
    impact: Optional[str]
    hgvsc: Optional[str]
    hgvsp: Optional[str]
    clinvar_id: Optional[str]
    clinvar_classification: Optional[str]
    clinvar_review_status: Optional[str]
    clinvar_condition: Optional[str]
    is_definitive_pathogenic: bool


def _rate_limited_get(url: str) -> Optional[dict]:
    global _last_request_at
    elapsed = time.monotonic() - _last_request_at
    if elapsed < _MIN_REQUEST_INTERVAL_S:
        time.sleep(_MIN_REQUEST_INTERVAL_S - elapsed)
    _last_request_at = time.monotonic()
    try:
        resp = requests.get(url, headers={"Content-Type": "application/json"}, timeout=10)
        resp.raise_for_status()
        return resp.json()
    except requests.RequestException:
        return None


@lru_cache(maxsize=4096)
def hgvs_for(variant: str, transcript_id: str) -> tuple[Optional[str], Optional[str]]:
    """(hgvsc, hgvsp) for one variant on one transcript -- (None, None) if VEP
    has nothing for that pair (e.g. a non-coding consequence has no hgvsp) or
    the request failed."""
    url = (
        f"https://rest.ensembl.org/vep/human/hgvs/{urllib.parse.quote(variant)}"
        "?content-type=application/json&hgvs=1"
    )
    data = _rate_limited_get(url)
    if not data:
        return None, None
    for tc in data[0].get("transcript_consequences", []):
        if tc.get("transcript_id") == transcript_id:
            return tc.get("hgvsc"), tc.get("hgvsp")
    return None, None


@lru_cache(maxsize=4096)
def clinvar_classification_for(clinvar_id: str) -> tuple[Optional[str], Optional[str], Optional[str]]:
    """(classification, review_status, condition) straight from ClinVar's own
    current record for this ID -- (None, None, None) if the lookup failed.
    `condition` is the first trait name that isn't ClinVar's own generic
    "not specified"/"not provided" placeholders, if any."""
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={clinvar_id}&retmode=json"
    data = _rate_limited_get(url)
    if not data:
        return None, None, None
    record = data.get("result", {}).get(clinvar_id)
    if not record:
        return None, None, None
    germline = record.get("germline_classification") or {}
    classification = germline.get("description")
    review_status = germline.get("review_status")
    condition = None
    for trait in germline.get("trait_set", []):
        name = trait.get("trait_name")
        if name and name.lower() not in ("not specified", "not provided"):
            condition = name
            break
    return classification, review_status, condition


def _parse_entry(entry: str) -> Optional[VariantRecord]:
    m = _ENTRY_RE.match(entry)
    if not m:
        return None
    variant, info = m.group(1), m.group(2)

    record: VariantRecord = {
        "variant": variant, "transcript": None, "consequence_terms": [], "impact": None,
        "hgvsc": None, "hgvsp": None, "clinvar_id": None, "clinvar_classification": None,
        "clinvar_review_status": None, "clinvar_condition": None, "is_definitive_pathogenic": False,
    }
    for part in info.split(":"):
        cm = _CLINVAR_PART_RE.match(part)
        if cm:
            record["clinvar_id"] = cm.group(1)
            continue
        tm = _TRANSCRIPT_PART_RE.match(part)
        if tm:
            record["transcript"] = tm.group(1)
            record["consequence_terms"] = tm.group(2).split(";")
            record["impact"] = tm.group(3)

    if record["clinvar_id"]:
        classification, review_status, condition = clinvar_classification_for(record["clinvar_id"])
        record["clinvar_classification"] = classification
        record["clinvar_review_status"] = review_status
        record["clinvar_condition"] = condition
        record["is_definitive_pathogenic"] = (classification or "").strip().lower() in _DEFINITIVE_CLASSIFICATIONS
    # else: no ClinVar entry at all -- this variant was flagged purely on
    # VEP's own HIGH-impact call (see interpret_paraphase.py's
    # process_vcf_file), a predicted consequence, not a ClinVar-confirmed
    # pathogenic classification. is_definitive_pathogenic stays False: the
    # user's own ask was specifically "check against ClinVar to see if the
    # variant is path or likely path", and claiming "definitive" here would
    # overstate a computational prediction as clinical curation.

    if record["transcript"]:
        hgvsc, hgvsp = hgvs_for(variant, record["transcript"])
        record["hgvsc"] = hgvsc
        record["hgvsp"] = hgvsp

    return record


def parse_pathogenic_variants(cell: str) -> list[VariantRecord]:
    """One structured, ClinVar-and-VEP-verified record per variant in a
    Pathogenic_variants_annotated cell. Empty list for cells with no real
    variant (empty, "None annotated", "known_nonfunctional")."""
    if not cell or cell in ("None annotated", "known_nonfunctional"):
        return []
    records = [_parse_entry(entry) for entry in cell.split(",")]
    return [r for r in records if r is not None]
