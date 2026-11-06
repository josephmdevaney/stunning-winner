"""Interpret PacBio Paraphase's per-gene output for the carrier panel's non-repeat
targets: CNV, gene-conversion, and structural-variant genes a TRGT VCF can't
express (see `carrier_screen.non_repeat_targets.json`'s own `mechanism` field for
why each one needs a different tool).

Every finding here is read directly off a field Paraphase (or, for F8, PacBio's
dedicated intron-1/22 inversion caller; or, for HBA1/2's large deletions, the
`sawfish` structural-variant caller; or, for named pathogenic variants, the
"havanno" per-haplotype variant annotation) itself already asserts -- `total_cn`,
`gene_cn`, `genotype`, `fusions_called`, `sv_called`, the inversion caller's own
GT/SUPP, sawfish's own GT/SVTYPE, or havanno's own `variants.pathogenic` list --
nothing here re-derives a pathogenic call from raw haplotype sequences. Where the
tool itself couldn't attribute a copy-number change to a specific gene (`gene_cn`
is null, as with this panel's real CYP21A2/TNXB region), that's reported as an
unresolved/atypical finding rather than guessed as benign or pathogenic.
"""

from __future__ import annotations

import gzip
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .strchive import STATUS_COLORS, STATUS_INTERMEDIATE, STATUS_NORMAL, STATUS_PATHOGENIC, STATUS_UNCERTAIN


def _open_text(path: str | Path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt", encoding="utf-8")

# "Full expansion" (STATUS_PATHOGENIC) is STR-specific terminology -- none of this
# module's findings are repeat expansions (they're SNVs/indels, gene conversions,
# CNVs, or a structural inversion), so pathogenic findings here use ClinVar's own
# clinical-significance term instead. Same severity/color as STATUS_PATHOGENIC.
STATUS_PATHOGENIC_VARIANT = "Pathogenic"
_STATUS_COLORS = {**STATUS_COLORS, STATUS_PATHOGENIC_VARIANT: STATUS_COLORS[STATUS_PATHOGENIC]}

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
    zygosity: Optional[str] = None
    # Named variants this finding is built from (only populated for havanno-derived
    # findings -- CYP21A2/GBA's named pathogenic variant lists). Empty otherwise.
    variants: tuple[str, ...] = ()
    # Other raw fields the underlying tool asserted, kept verbatim for a "what did
    # the tool actually say" detail view (e.g. HBA's `genotype`/`sv_called`,
    # total_cn, GT/SUPP) -- never a re-derived or summarized value.
    raw_fields: dict[str, str] = field(default_factory=dict)


def load_paraphase_calls(json_path: str | Path) -> dict[str, dict]:
    with open(json_path, encoding="utf-8") as f:
        return json.load(f)


def load_havanno_annotations(json_path: str | Path) -> dict[str, dict]:
    """Parse PacBio's "havanno" per-haplotype variant annotation JSON: for each
    Paraphase gene key, a dict of haplotype-phase-call-id -> named pathogenic (and
    pseudogene) variants asserted on that call (e.g. Gauchian GBA1 calls, the
    curated CYP21A2 variant list). Returns the `annotations` section keyed by
    Paraphase gene key (`cyp21`, `gba`, ...), or `{}` if unavailable."""
    with open(json_path, encoding="utf-8") as f:
        data = json.load(f)
    return data.get("annotations", {})


def load_sawfish_calls(vcf_path: str | Path) -> list[dict]:
    """Parse sawfish structural-variant VCF records. Used as a fallback for
    HBA1/HBA2's large deletions, which Paraphase's own copy-number caller can
    fail to resolve to a genotype (see `carrier_screen.non_repeat_targets.json`'s
    HBA1/2 entry: "small deletion via Paraphase, large deletion via sawfish")."""
    calls: list[dict] = []
    with _open_text(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, _id, _ref, _alt, _qual, filter_, info, fmt, sample_field = fields[:10]
            if filter_ != "PASS":
                continue
            info_d = dict(kv.split("=", 1) for kv in info.split(";") if "=" in kv)
            fmt_d = dict(zip(fmt.split(":"), sample_field.split(":")))
            calls.append({
                "chrom": chrom,
                "pos": int(pos),
                "end": int(info_d["END"]) if "END" in info_d else int(pos),
                "svtype": info_d.get("SVTYPE"),
                "svlen": int(info_d["SVLEN"]) if "SVLEN" in info_d else None,
                "gt": fmt_d.get("GT"),
            })
    return calls


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


def _finding(
    status: str,
    note: str,
    zygosity: Optional[str] = None,
    variants: tuple[str, ...] = (),
    raw_fields: Optional[dict[str, str]] = None,
) -> ParaphaseFinding:
    return ParaphaseFinding(
        status=status, status_color=_STATUS_COLORS[status], note=note, zygosity=zygosity,
        variants=variants, raw_fields=raw_fields or {},
    )


def _gt_zygosity(gt: Optional[str]) -> Optional[str]:
    """Zygosity directly off a VCF-style GT string -- hemizygous single-allele
    (`1`), diploid homozygous (`1/1`), diploid heterozygous (`0/1`), phased
    (`0|1`), or None if reference/missing/unparseable."""
    if not gt:
        return None
    alleles = gt.replace("|", "/").split("/")
    if any(a == "." for a in alleles):
        return None
    if not any(a != "0" for a in alleles):
        return None
    if len(alleles) == 1:
        return "Hemizygous"
    return "Homozygous" if len(set(alleles)) == 1 else "Heterozygous"


def _classify_hba(
    call: dict,
    havanno_gene_annotations: Optional[dict[str, dict]] = None,
    sawfish_calls: Optional[list[dict]] = None,
) -> ParaphaseFinding:
    havanno_finding = _classify_havanno_gene(havanno_gene_annotations or {}, call)
    if havanno_finding is not None:
        return havanno_finding

    genotype = call.get("genotype")
    sv_called = call.get("sv_called") or {}
    raw_fields = {"genotype": genotype} if genotype else {}
    if sv_called:
        raw_fields["sv_called"] = ", ".join(f"{k}: {v}" for k, v in sv_called.items())

    if genotype == "aa/aa":
        return _finding(
            STATUS_NORMAL, "Paraphase: genotype aa/aa (4 of 4 alpha-globin copies present, no deletion)",
            zygosity="Homozygous", raw_fields=raw_fields,
        )
    if genotype:
        return _finding(
            STATUS_UNCERTAIN,
            f"Paraphase reported genotype {genotype!r} (not the normal aa/aa) -- clinical review recommended",
            raw_fields=raw_fields,
        )
    dels = [c for c in (sawfish_calls or []) if c.get("svtype") == "DEL" and _gt_zygosity(c.get("gt"))]
    if dels:
        parts = ", ".join(f"{c['chrom']}:{c['pos']}-{c['end']} ({c['svlen']} bp, GT={c['gt']})" for c in dels)
        zygosities = {_gt_zygosity(c["gt"]) for c in dels}
        zygosity = zygosities.pop() if len(zygosities) == 1 else None
        return _finding(
            STATUS_UNCERTAIN,
            "Paraphase did not report a resolved HBA1/HBA2 genotype, and sawfish found a deletion overlapping "
            f"the region ({parts}) that cannot be attributed to HBA1 vs. HBA2 from this data -- clinical review "
            "recommended",
            zygosity=zygosity,
            raw_fields={**raw_fields, "sawfish_deletions": parts},
        )
    return _finding(STATUS_UNCERTAIN, "Paraphase did not report a resolved HBA1/HBA2 genotype", raw_fields=raw_fields)


def _classify_smn1(call: dict) -> ParaphaseFinding:
    raw_fields = {
        k: str(call[k]) for k in ("smn1_cn", "smn2_cn", "smn_del78_cn") if call.get(k) is not None
    }
    smn1_cn = call.get("smn1_cn")
    if smn1_cn is None:
        return _finding(STATUS_UNCERTAIN, "Paraphase did not report an SMN1 copy number", raw_fields=raw_fields)
    if smn1_cn == 0:
        return _finding(
            STATUS_PATHOGENIC_VARIANT,
            "Paraphase: 0 SMN1 copies detected -- consistent with spinal muscular atrophy, refer for confirmation",
            raw_fields=raw_fields,
        )
    if smn1_cn == 1:
        return _finding(
            STATUS_INTERMEDIATE, "Paraphase: 1 SMN1 copy detected -- heterozygous SMA carrier",
            raw_fields=raw_fields,
        )
    return _finding(STATUS_NORMAL, f"Paraphase: {smn1_cn} SMN1 copies detected (normal)", raw_fields=raw_fields)


def _classify_generic(gene_key: str, call: dict) -> ParaphaseFinding:
    expected = _EXPECTED_TOTAL_CN[gene_key]
    fusions = call.get("fusions_called")
    if fusions:
        return _finding(
            STATUS_PATHOGENIC_VARIANT, f"Paraphase flagged a fusion/gene-conversion event: {fusions}",
            raw_fields={"fusions_called": str(fusions)},
        )

    gene_cn = call.get("gene_cn")
    if gene_cn is not None:
        raw_fields = {"gene_cn": str(gene_cn)}
        if gene_cn == expected:
            return _finding(STATUS_NORMAL, f"Paraphase: {gene_cn} gene copies detected (normal)", raw_fields=raw_fields)
        return _finding(
            STATUS_UNCERTAIN, f"Paraphase: {gene_cn} gene copies detected (expected {expected})",
            raw_fields=raw_fields,
        )

    total_cn = call.get("total_cn")
    if total_cn is None:
        return _finding(STATUS_UNCERTAIN, "Paraphase did not report a copy number for this region")
    raw_fields = {"total_cn": str(total_cn)}
    unit = "gene+pseudogene copies" if _HAS_PSEUDOGENE_PAIR[gene_key] else "gene copies"
    if total_cn == expected:
        return _finding(STATUS_NORMAL, f"Paraphase: {total_cn} total {unit} detected (normal)", raw_fields=raw_fields)
    attribution = (
        "could not attribute the difference to a specific gene vs. pseudogene copy"
        if _HAS_PSEUDOGENE_PAIR[gene_key]
        else "could not attribute the difference further"
    )
    return _finding(
        STATUS_UNCERTAIN,
        f"Paraphase: {total_cn} total {unit} detected (expected {expected}) but {attribution} "
        "-- clinical review recommended",
        raw_fields=raw_fields,
    )


def _havanno_zygosity(entries_with_variants: dict, call: Optional[dict]) -> Optional[str]:
    """Approximate zygosity by weighing how many of the gene's total physical
    copies (Paraphase's own `total_cn`, adjusted for `two_copy_haplotypes` --
    see docs/vcf.md: "haplotypes present in two copies are reported as separate
    columns... labelled with the '_cp2' suffix") carry the pathogenic variant
    set. Paraphase couldn't phase this locus into 2 chromosomes for either real
    sample this project has (`alleles_final`/`haplotype_links` come back empty),
    so this is the most precise zygosity call available without guessing which
    haplotype-phase calls pair onto which chromosome."""
    if not call:
        return None
    total_copies = call.get("total_cn")
    if not total_copies:
        return None
    two_copy_haplotypes = set(call.get("two_copy_haplotypes") or [])
    pathogenic_copies = sum(2 if hap_id in two_copy_haplotypes else 1 for hap_id in entries_with_variants)
    if pathogenic_copies == total_copies:
        return "Homozygous"
    if 0 < pathogenic_copies < total_copies:
        return "Heterozygous"
    return None


def _classify_havanno_gene(gene_annotations: dict[str, dict], call: Optional[dict] = None) -> Optional[ParaphaseFinding]:
    """Build a Pathogenic finding from havanno's named-variant annotations, or
    None if no haplotype-phase call for this gene carries any named pathogenic
    variant (caller should then fall back to the total_cn-based classifier).

    A single physical haplotype is often represented by several redundant
    haplotype-phase calls (e.g. `cyp21_hap2`/`hap5`/`hap7` all carrying the same
    variant list); those are collapsed by reporting the distinct variant sets
    actually seen, not the raw call count.
    """
    entries_with_variants = {
        hap_id: tuple(sorted(set(entry.get("variants", {}).get("pathogenic", []))))
        for hap_id, entry in gene_annotations.items()
        if entry.get("num_pathogenic_variants", 0) > 0
    }
    if not entries_with_variants:
        return None

    distinct_variant_sets = sorted(set(entries_with_variants.values()), key=lambda s: (-len(s), s))
    union_variants = sorted({v for variants in entries_with_variants.values() for v in variants})
    n_calls = len(entries_with_variants)
    n_total = len(gene_annotations)
    set_note = (
        "1 distinct variant set"
        if len(distinct_variant_sets) == 1
        else f"{len(distinct_variant_sets)} distinct (overlapping) variant sets"
    )
    return _finding(
        STATUS_PATHOGENIC_VARIANT,
        f"Named pathogenic variant annotation flagged {len(union_variants)} variant(s) "
        f"({set_note}) on {n_calls} of {n_total} haplotype-phase calls: {', '.join(union_variants)}",
        zygosity=_havanno_zygosity(entries_with_variants, call),
        variants=tuple(union_variants),
        raw_fields={"haplotype_phase_calls": f"{n_calls} of {n_total}"},
    )


def _classify_f8(inv_calls: dict[str, dict]) -> ParaphaseFinding:
    raw_fields = {
        record_id: f"GT={c.get('gt')}, SUPP={c.get('supp')}" for record_id, c in sorted(inv_calls.items())
    }
    inverted = {
        record_id: c for record_id, c in inv_calls.items()
        if _gt_zygosity(c.get("gt")) or (c.get("supp") or 0) > 0
    }
    if inverted:
        zygosities = {_gt_zygosity(c.get("gt")) for c in inverted.values()} - {None}
        zygosity = zygosities.pop() if len(zygosities) == 1 else None
        return _finding(
            STATUS_PATHOGENIC_VARIANT,
            f"Intron 1/22 inversion breakpoint detected at {', '.join(sorted(inverted))} "
            "(pathogenic structural rearrangement causing Hemophilia A)",
            zygosity=zygosity,
            raw_fields=raw_fields,
        )
    if inv_calls:
        return _finding(
            STATUS_NORMAL,
            "No intron 1/22 inversion breakpoint detected (dedicated F8 inversion caller, both probes reference)",
            raw_fields=raw_fields,
        )
    return _finding(STATUS_UNCERTAIN, "No F8 inversion caller output available for this sample")


def classify_gene(
    gene: str,
    paraphase_calls: dict[str, dict],
    f8_inv_calls: Optional[dict[str, dict]],
    havanno_annotations: Optional[dict[str, dict]] = None,
    sawfish_calls: Optional[list[dict]] = None,
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
        return _classify_hba(call, (havanno_annotations or {}).get("hba"), sawfish_calls)
    if key == "smn1":
        return _classify_smn1(call)

    # GBA and CYP21A2 each have a same-size pseudogene pair: havanno's per-haplotype named-variant
    # pattern-match can't tell a real pathogenic variant on the functional gene copy from the
    # pseudogene's own normal sequence (which carries many of the same "differentiating site"
    # alleles its curated list treats as pathogenic). Trusting it flagged ~50% of the real carrier
    # cohort, against a true carrier frequency of a few percent -- so skip it for those genes and
    # fall through to the more conservative copy-number/fusion reasoning in _classify_generic.
    if havanno_annotations and not _HAS_PSEUDOGENE_PAIR.get(key):
        havanno_finding = _classify_havanno_gene(havanno_annotations.get(key, {}), call)
        if havanno_finding is not None:
            return havanno_finding

    return _classify_generic(key, call)
