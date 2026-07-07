"""Parser for PacBio TRGT (Tandem Repeat Genotyping Tool) VCF output.

Reference: https://github.com/PacificBiosciences/trgt/blob/main/docs/vcf_files.md

TRGT writes one record per repeat region with:

INFO
    TRID    identifier of the repeat region from the input catalog
    END     ending position of the repeat region
    MOTIFS  comma-separated motif units observed at the locus
    STRUC   expected locus structure (deprecated)

FORMAT (per sample)
    GT      diploid genotype, allele indices into REF(0)/ALT(1..) in *declared* order
    AL      allele length in bp, one entry per declared allele (REF+ALT), in allele-index
            order -- NOT in GT/phase display order
    ALLR    length confidence interval ("low-high") per declared allele, same order as AL
    SD      number of spanning reads supporting each declared allele, same order as AL
    MC      motif copy count(s) per declared allele, same order as AL; a single allele's
            value may itself be underscore-joined sub-counts for compound repeat structures
            (e.g. "17_2"), which we sum to a total repeat-unit count
    MS      run-length segmentation string per allele (kept raw, not parsed further)
    AP      allele purity (fraction of motif-consistent bases), same order as AL
    AM      mean 5mCpG methylation level per allele ("." when unavailable), same order as AL
    PS      phase set identifier, present when GT is phased

Critically, AL/ALLR/SD/MC/MS/AP/AM are indexed by *allele index* (0=REF, 1=first ALT, ...),
not by haplotype/phase position -- so reading off "the two called alleles" requires looking
up GT's two allele indices into these arrays, not just taking their first two values.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional

from .motifs import is_cyclic_rotation


def _open_text(path: str | Path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt", encoding="utf-8")


def _parse_info(info: str) -> dict[str, str]:
    fields: dict[str, str] = {}
    if info in ("", "."):
        return fields
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            fields[k] = v
        else:
            fields[kv] = ""
    return fields


def _parse_gt(gt: str) -> tuple[Optional[int], Optional[int], bool, bool]:
    """Parse GT into (hap1_idx, hap2_idx, phased, hemizygous).

    A GT with no `/` or `|` separator (e.g. `1`, not `1/1`) is a real, valid TRGT
    output -- trgt's own hemizygous calling mode for chrX/chrY loci in a male
    sample emits exactly one allele index and a single-slot FORMAT array, not a
    two-allele diploid call. Treating that single value as `hap1_idx` with
    `hap2_idx=None` (rather than padding a second, nonexistent haplotype) is what
    `hemizygous=True` signals to `TrgtCall._hap_values`.
    """
    phased = "|" in gt
    sep = "|" if phased else ("/" if "/" in gt else None)
    if sep is None:
        value = None if gt == "." else int(gt)
        return value, None, False, True
    parts = gt.split(sep)
    alleles = tuple(None if p == "." else int(p) for p in parts)
    while len(alleles) < 2:
        alleles = alleles + (None,)
    return alleles[0], alleles[1], phased, False


def _parse_int_list(raw: str) -> list[Optional[int]]:
    if raw in ("", "."):
        return []
    return [None if v == "." else int(v) for v in raw.split(",")]


def _parse_range_list(raw: str) -> list[Optional[tuple[int, int]]]:
    if raw in ("", "."):
        return []
    out: list[Optional[tuple[int, int]]] = []
    for v in raw.split(","):
        if v == "." or "-" not in v:
            out.append(None)
            continue
        lo, hi = v.split("-", 1)
        out.append((int(lo), int(hi)))
    return out


def _parse_float_list(raw: str) -> list[Optional[float]]:
    if raw in ("", "."):
        return []
    return [None if v == "." else float(v) for v in raw.split(",")]


def _parse_mc_segments(raw: str) -> list[Optional[list[int]]]:
    """Parse MC into per-allele lists of per-STRUC-segment repeat-unit counts.

    For a compound locus like HTT's `(CAG)nCAACAG(CCG)n` (MOTIFS=CAG,CCG), MC reports
    one underscore-joined count per STRUC segment, e.g. "18_8" = 18 CAG units + 8 CCG
    units. Keeping these separate (rather than summing) matters clinically: STRchive's
    thresholds are almost always about one specific motif's count, not the sum across
    an unrelated flanking repeat.
    """
    if raw in ("", "."):
        return []
    out: list[Optional[list[int]]] = []
    for v in raw.split(","):
        if v == ".":
            out.append(None)
            continue
        out.append([int(part) for part in v.split("_") if part != "."])
    return out


@dataclass
class TrgtCall:
    chrom: str
    pos: int
    end: int
    trid: str
    motifs: list[str]
    struc: str
    ref: str
    alt: list[str]
    hap1_idx: Optional[int]
    hap2_idx: Optional[int]
    phased: bool
    hemizygous: bool
    allele_length: list[Optional[int]]
    allele_length_range: list[Optional[tuple[int, int]]]
    allele_support: list[Optional[int]]
    motif_count: list[Optional[list[int]]]
    allele_purity: list[Optional[float]]
    allele_methylation: list[Optional[float]]

    def _hap_values(self, per_allele: list) -> tuple[Optional[object], Optional[object]]:
        """Map a FORMAT array to (hap1_value, hap2_value).

        TRGT's AL/ALLR/SD/MC/MS/AP/AM fields are NOT sized to cover every declared
        REF+ALT allele (Number=R) -- they always have exactly one slot per haplotype
        (2 for a diploid call), ordered by ascending *original* allele index of the
        two haplotypes actually called, regardless of how many ALT alleles exist at
        the site. E.g. GT=1/2 with no haplotype calling REF yields a 2-element array
        for indices {1,2} at positions [0,1] -- not a 3-element REF+ALT1+ALT2 array.
        For a homozygous call (both haplotypes the same index), the two slots are
        still independent per-haplotype-copy values (e.g. SD legitimately differs
        between them), not a single deduplicated value.
        """
        if self.hemizygous:
            # Exactly one physical haplotype was called (male chrX/chrY), so the
            # FORMAT array has exactly one slot -- there's no "rank" to compute.
            if self.hap1_idx is None or not per_allele:
                return None, None
            return per_allele[0], None

        if self.hap1_idx is None or self.hap2_idx is None:
            return None, None
        lo = min(self.hap1_idx, self.hap2_idx)

        def get(idx: int, other_idx: int, is_first_of_tie: bool):
            # Rank: the lower original allele-index gets slot 0, the higher gets slot 1.
            # When both haplotypes share the same index, preserve haplotype order
            # (hap1 -> slot 0, hap2 -> slot 1) since each still has its own real value.
            if idx < other_idx or (idx == other_idx and is_first_of_tie):
                rank = 0
            else:
                rank = 1
            if rank >= len(per_allele):
                return None
            return per_allele[rank]

        return (
            get(self.hap1_idx, self.hap2_idx, True),
            get(self.hap2_idx, self.hap1_idx, False),
        )

    @property
    def genotyped(self) -> bool:
        if self.hemizygous:
            return self.hap1_idx is not None
        return self.hap1_idx is not None and self.hap2_idx is not None

    @property
    def allele_lengths_bp(self) -> tuple[Optional[int], Optional[int]]:
        return self._hap_values(self.allele_length)

    @property
    def allele_length_ranges(self) -> tuple[Optional[tuple[int, int]], Optional[tuple[int, int]]]:
        return self._hap_values(self.allele_length_range)

    @property
    def allele_depths(self) -> tuple[Optional[int], Optional[int]]:
        return self._hap_values(self.allele_support)

    def repeat_units(self, target_motif: Optional[str] = None) -> tuple[Optional[int], Optional[int]]:
        """Repeat-unit counts per haplotype for a specific motif.

        Compound loci (e.g. HTT's `(CAG)nCAACAG(CCG)n`) report one MC sub-count per
        STRUC segment. When `target_motif` is given, this picks out that segment's
        count specifically (matching by cyclic rotation, since motif reading-frame is
        arbitrary -- STRchive may record "AAG" where the VCF's own MOTIFS says "GAA").
        Without a target motif (e.g. for generic, non-disease catalog regions), falls
        back to summing all segments. Falls back further to AL / motif length if MC is
        entirely absent.
        """
        segs1, segs2 = self._hap_values(self.motif_count)

        def pick(segs: Optional[list[int]]) -> Optional[int]:
            if segs is None:
                return None
            if target_motif is not None and len(self.motifs) == len(segs):
                for i, m in enumerate(self.motifs):
                    if is_cyclic_rotation(target_motif, m):
                        return segs[i]
            return sum(segs)

        u1, u2 = pick(segs1), pick(segs2)
        if u1 is not None and u2 is not None:
            return u1, u2

        al1, al2 = self.allele_lengths_bp
        motif_len = len(target_motif) if target_motif else (len(self.motifs[0]) if self.motifs and self.motifs[0] else None)
        if motif_len:
            u1 = u1 if u1 is not None else (round(al1 / motif_len) if al1 is not None else None)
            u2 = u2 if u2 is not None else (round(al2 / motif_len) if al2 is not None else None)
        return u1, u2

    def dominant_motifs(self) -> tuple[Optional[str], Optional[str]]:
        """The motif with the largest MC segment count for each haplotype.

        For loci like RFC1/FAME genes where multiple candidate motifs can occupy the
        same repeat tract, pathogenicity hinges on *which* motif expanded, not just how
        long the tract is. Returns None for a haplotype whose MC data is missing, has no
        nonzero segment, or doesn't line up 1:1 with the declared MOTIFS list.
        """
        segs1, segs2 = self._hap_values(self.motif_count)

        def pick(segs: Optional[list[int]]) -> Optional[str]:
            if not segs or len(segs) != len(self.motifs):
                return None
            best_idx = max(range(len(segs)), key=lambda i: segs[i])
            if segs[best_idx] <= 0:
                return None
            return self.motifs[best_idx]

        return pick(segs1), pick(segs2)

    @property
    def allele_methylation_pct(self) -> tuple[Optional[float], Optional[float]]:
        return self._hap_values(self.allele_methylation)

    @property
    def allele_purity_frac(self) -> tuple[Optional[float], Optional[float]]:
        return self._hap_values(self.allele_purity)

    @property
    def total_depth(self) -> Optional[int]:
        """Sum of supporting reads across both haplotypes.

        SD reports one independent value per haplotype copy -- even a homozygous call
        (both haplotypes the same allele) gets two distinct real read-support counts,
        not one shared value -- so this is a plain sum of both slots.
        """
        d1, d2 = self.allele_depths
        vals = [d for d in (d1, d2) if d is not None]
        return sum(vals) if vals else None


def sample_names(path: str | Path) -> list[str]:
    with _open_text(path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                return cols[9:]
    return []


def iter_trgt_records(path: str | Path, sample: str | None = None) -> Iterator[TrgtCall]:
    """Yield one TrgtCall per VCF record for the given sample (first sample if omitted)."""
    with _open_text(path) as f:
        sample_idx = 0
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                samples = cols[9:]
                if sample is not None:
                    if sample not in samples:
                        raise ValueError(f"sample {sample!r} not found in VCF (have: {samples})")
                    sample_idx = samples.index(sample)
                continue

            fields = line.rstrip("\n").split("\t")
            chrom, pos, _id, ref, alt, _qual, _filter, info, fmt = fields[:9]
            sample_field = fields[9 + sample_idx]

            info_d = _parse_info(info)
            fmt_keys = fmt.split(":")
            fmt_vals = sample_field.split(":")
            fmt_d = dict(zip(fmt_keys, fmt_vals))

            gt = fmt_d.get("GT", "./.")
            hap1_idx, hap2_idx, phased, hemizygous = _parse_gt(gt)

            yield TrgtCall(
                chrom=chrom,
                pos=int(pos),
                end=int(info_d.get("END", pos)),
                trid=info_d.get("TRID", ""),
                hemizygous=hemizygous,
                motifs=(info_d.get("MOTIFS", "").split(",") if info_d.get("MOTIFS") else []),
                struc=info_d.get("STRUC", ""),
                ref=ref,
                alt=([] if alt in (".", "") else alt.split(",")),
                hap1_idx=hap1_idx,
                hap2_idx=hap2_idx,
                phased=phased,
                allele_length=_parse_int_list(fmt_d.get("AL", "")),
                allele_length_range=_parse_range_list(fmt_d.get("ALLR", "")),
                allele_support=_parse_int_list(fmt_d.get("SD", "")),
                motif_count=_parse_mc_segments(fmt_d.get("MC", "")),
                allele_purity=_parse_float_list(fmt_d.get("AP", "")),
                allele_methylation=_parse_float_list(fmt_d.get("AM", "")),
            )


class TrgtIndex:
    """Chromosome-bucketed, start-sorted index of TRGT calls for coordinate/TRID lookup."""

    def __init__(self, calls: list[TrgtCall]):
        self.by_trid: dict[str, TrgtCall] = {}
        self._by_chrom: dict[str, list[TrgtCall]] = {}
        for call in calls:
            if call.trid:
                self.by_trid[call.trid] = call
            self._by_chrom.setdefault(call.chrom, []).append(call)
        for chrom_calls in self._by_chrom.values():
            chrom_calls.sort(key=lambda c: c.pos)
        self._starts: dict[str, list[int]] = {
            chrom: [c.pos for c in chrom_calls] for chrom, chrom_calls in self._by_chrom.items()
        }

    def find_overlap(self, chrom: str, start: int, stop: int) -> Optional[TrgtCall]:
        """Find a call overlapping [start, stop) on chrom, if any, by binary search on pos."""
        import bisect

        starts = self._starts.get(chrom)
        if not starts:
            return None
        calls = self._by_chrom[chrom]
        i = bisect.bisect_right(starts, stop)
        # Scan backward from the first call starting at or before `stop` looking for overlap;
        # repeat regions in these catalogs are non-overlapping and sorted, so a short
        # backward scan is sufficient.
        for j in range(i - 1, max(i - 5, -1), -1):
            call = calls[j]
            if call.end < start:
                break
            if call.pos <= stop and call.end >= start:
                return call
        return None

    @classmethod
    def from_vcf(cls, path: str | Path, sample: str | None = None) -> "TrgtIndex":
        return cls(list(iter_trgt_records(path, sample=sample)))
