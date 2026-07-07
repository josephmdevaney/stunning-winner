"""Interpret PacBio Paraphase's per-gene output for the carrier panel's non-repeat
targets: CNV, gene-conversion, and structural-variant genes a TRGT VCF can't
express (see `carrier_screen.non_repeat_targets.json`'s own `mechanism` field for
why each one needs a different tool).

Every finding here is read directly off a field Paraphase (or, for F8, PacBio's
dedicated intron-1/22 inversion caller) itself already asserts -- `total_cn`,
`gene_cn`, `genotype`, `fusions_called`, `sv_called`, or the inversion caller's own
GT/SUPP -- nothing here re-derives a pathogenic call from raw haplotype sequences.
Where the tool itself couldn't attribute a copy-number change to a specific gene
(`gene_cn` is null, as with this panel's real CYP21A2/TNXB region), that's reported
as an unresolved/atypical finding rather than guessed as benign or pathogenic.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from .strchive import STATUS_COLORS, STATUS_INTERMEDIATE, STATUS_NORMAL, STATUS_PATHOGENIC, STATUS_UNCERTAIN

# Expected total copy number for the genes whose Paraphase call is a plain total_cn
# comparison, and whether that count spans a gene+pseudogene pair (GBA/GBAP1,
# CYP21A2/CYP21A1P -- 4 total copies in a normal diploid sample) or just the gene
# itself (HBB, RPGR -- no same-size pseudogene, so 2).
_EXPECTED_TOTAL_CN = {"gba": 4, "cyp21": 4, "hbb": 2, "rpgr": 2}
_HAS_PSEUDOGENE_PAIR = {"gba": True, "cyp21": True, "hbb": False, "rpgr": False}

# Maps a carrier_screen.non_repeat_targets.json `gene` value to its Paraphase JSON
# key. TNXB is deliberately absent: PacBio's own carrier pipeline has no caller for
# it yet (see that target's own "tool": "Not yet implemented"), so it always stays
# "Not applicable" regardless of what Paraphase data is available.
GENE_TO_PARAPHASE_KEY = {
    "GBA": "gba",
    "HBB": "hbb",
    "HBA1/2": "hba",
    "SMN1/2": "smn1",
    "CYP21A2": "cyp21",
    "RPGR": "rpgr",
    "F8": "f8",  # handled via the dedicated inversion VCF, not the generic JSON
}


@dataclass(frozen=True)
class ParaphaseFinding:
    status: str
    status_color: str
    note: str


def load_paraphase_calls(json_path: str | Path) -> dict[str, dict]:
    with open(json_path, encoding="utf-8") as f:
        return json.load(f)


def load_f8_inversion_calls(vcf_path: str | Path) -> dict[str, dict]:
    """Parse PacBio's dedicated F8 intron-1/22 inversion caller VCF.

    One record per breakpoint (IDs `f8inv1`, `f8inv22`); GT=0 with SUPP=0 means no
    reads support the inversion breakpoint at that probe, GT=1 means the inversion
    was called.
    """
    calls: dict[str, dict] = {}
    with open(vcf_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            record_id, _ref, _alt, _qual, _filter, info, fmt, sample_field = (
                fields[2], fields[3], fields[4], fields[5], fields[6], fields[7], fields[8], fields[9],
            )
            fmt_d = dict(zip(fmt.split(":"), sample_field.split(":")))
            info_d = dict(kv.split("=", 1) for kv in info.split(";") if "=" in kv)
            calls[record_id] = {
                "gt": fmt_d.get("GT"),
                "supp": int(fmt_d["SUPP"]) if fmt_d.get("SUPP", ".") != "." else None,
                "evidence": info_d.get("EVIDENCE"),
            }
    return calls


def _finding(status: str, note: str) -> ParaphaseFinding:
    return ParaphaseFinding(status=status, status_color=STATUS_COLORS[status], note=note)


def _classify_hba(call: dict) -> ParaphaseFinding:
    genotype = call.get("genotype")
    if genotype == "aa/aa":
        return _finding(STATUS_NORMAL, "Paraphase: genotype aa/aa (4 of 4 alpha-globin copies present, no deletion)")
    if genotype:
        return _finding(
            STATUS_UNCERTAIN,
            f"Paraphase reported genotype {genotype!r} (not the normal aa/aa) -- clinical review recommended",
        )
    return _finding(STATUS_UNCERTAIN, "Paraphase did not report a resolved HBA1/HBA2 genotype")


def _classify_smn1(call: dict) -> ParaphaseFinding:
    smn1_cn = call.get("smn1_cn")
    if smn1_cn is None:
        return _finding(STATUS_UNCERTAIN, "Paraphase did not report an SMN1 copy number")
    if smn1_cn == 0:
        return _finding(
            STATUS_PATHOGENIC,
            "Paraphase: 0 SMN1 copies detected -- consistent with spinal muscular atrophy, refer for confirmation",
        )
    if smn1_cn == 1:
        return _finding(STATUS_INTERMEDIATE, "Paraphase: 1 SMN1 copy detected -- heterozygous SMA carrier")
    return _finding(STATUS_NORMAL, f"Paraphase: {smn1_cn} SMN1 copies detected (normal)")


def _classify_generic(gene_key: str, call: dict) -> ParaphaseFinding:
    expected = _EXPECTED_TOTAL_CN[gene_key]
    fusions = call.get("fusions_called")
    if fusions:
        return _finding(STATUS_PATHOGENIC, f"Paraphase flagged a fusion/gene-conversion event: {fusions}")

    gene_cn = call.get("gene_cn")
    if gene_cn is not None:
        if gene_cn == expected:
            return _finding(STATUS_NORMAL, f"Paraphase: {gene_cn} gene copies detected (normal)")
        return _finding(STATUS_UNCERTAIN, f"Paraphase: {gene_cn} gene copies detected (expected {expected})")

    total_cn = call.get("total_cn")
    if total_cn is None:
        return _finding(STATUS_UNCERTAIN, "Paraphase did not report a copy number for this region")
    unit = "gene+pseudogene copies" if _HAS_PSEUDOGENE_PAIR[gene_key] else "gene copies"
    if total_cn == expected:
        return _finding(STATUS_NORMAL, f"Paraphase: {total_cn} total {unit} detected (normal)")
    attribution = (
        "could not attribute the difference to a specific gene vs. pseudogene copy"
        if _HAS_PSEUDOGENE_PAIR[gene_key]
        else "could not attribute the difference further"
    )
    return _finding(
        STATUS_UNCERTAIN,
        f"Paraphase: {total_cn} total {unit} detected (expected {expected}) but {attribution} "
        "-- clinical review recommended",
    )


def _classify_f8(inv_calls: dict[str, dict]) -> ParaphaseFinding:
    inverted_at = [
        record_id for record_id, c in inv_calls.items() if c.get("gt") not in (None, "0") or (c.get("supp") or 0) > 0
    ]
    if inverted_at:
        return _finding(
            STATUS_PATHOGENIC,
            f"Intron 1/22 inversion breakpoint detected at {', '.join(sorted(inverted_at))} "
            "(pathogenic structural rearrangement causing Hemophilia A)",
        )
    if inv_calls:
        return _finding(
            STATUS_NORMAL,
            "No intron 1/22 inversion breakpoint detected (dedicated F8 inversion caller, both probes reference)",
        )
    return _finding(STATUS_UNCERTAIN, "No F8 inversion caller output available for this sample")


def classify_gene(
    gene: str, paraphase_calls: dict[str, dict], f8_inv_calls: Optional[dict[str, dict]]
) -> Optional[ParaphaseFinding]:
    """Build a finding for one carrier-panel non-repeat gene, or None if no real
    data is available for it (caller should keep the generic "Not applicable"
    placeholder in that case)."""
    key = GENE_TO_PARAPHASE_KEY.get(gene)
    if key is None:
        return None

    if key == "f8":
        return _classify_f8(f8_inv_calls or {})

    call = paraphase_calls.get(key)
    if call is None:
        return None
    if key == "hba":
        return _classify_hba(call)
    if key == "smn1":
        return _classify_smn1(call)
    return _classify_generic(key, call)
