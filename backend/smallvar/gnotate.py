"""Reader for slivar "gnotate" frequency zips, such as the gnomAD and HPRC resources in PacBio's
HiFi WDL resource bundle (``static_resources/GRCh38/slivar_gnotate/``).

A gnotate zip holds, per chromosome directory (``sli.var/<chrom>/``):

  gnotate-variant.bin   sorted little-endian uint64, one per variant
  gnotate-<field>.bin   float32, parallel to variant.bin (AF, AC and nhomalt are all float32)
  long-alleles.txt      variants whose alleles are too long to pack, with their values inline

Layout of the uint64 (worked out from the data and checked against a call set: 98% of
HG002's chr22 SNVs and indels are found, none of a shuffled-allele control):

  bits 36-63   position, 0-based (VCF POS - 1)
  bits 32-35   len(ref)
  bits 28-31   len(alt)
  bits  2-27   ref + alt bases, 2 bits each (A,C,G,T = 0..3), first base most significant
  bits  0-1    reserved (0)

Alleles with len(ref) + len(alt) > 13 do not fit and are looked up in long-alleles.txt, whose
rows are ``position (0-based) / ref / alt / filter / v1|v2|v3`` with the values in the
order of fields.txt. A variant absent from the file is simply not in it; for the PacBio
gnomAD file that means not seen at a PASS site in gnomAD v3.0 genomes.

The arrays are deflate-compressed inside the zip, so a chromosome is decompressed whole
(chr1 is about 550 MB for the variant array and one value array) and chromosomes are
processed one at a time.
"""

from __future__ import annotations

import bisect
import os
import sys
import zipfile
from array import array
from typing import Dict, List, Optional, Sequence, Tuple

from .normalize import canon_chrom

MAX_PACKED = 13  # len(ref) + len(alt) that still fits in the 26 base bits
_BASE = {"A": 0, "C": 1, "G": 2, "T": 3}


def encode_allele(ref: str, alt: str) -> Optional[int]:
    """The low 36 bits of a variant's code, or None if the alleles don't pack (too long, or
    not plain ACGT)."""
    n = len(ref) + len(alt)
    if n > MAX_PACKED or not ref or not alt:
        return None
    v = 0
    for ch in ref:
        b = _BASE.get(ch)
        if b is None:
            return None
        v = (v << 2) | b
    for ch in alt:
        b = _BASE.get(ch)
        if b is None:
            return None
        v = (v << 2) | b
    return (len(ref) << 32) | (len(alt) << 28) | (v << 2)


class Gnotate:
    def __init__(self, path: str):
        self.path = str(path)
        self._zip = zipfile.ZipFile(self.path)
        names = set(self._zip.namelist())
        self.chroms = self._read_lines("sli.var/chroms.txt")
        self.fields = self._read_lines("sli.var/fields.txt")
        self.args = " ".join(self._read_lines("sli.var/args.txt"))
        self.created = " ".join(self._read_lines("sli.var/message.txt"))
        self._dirs = {canon_chrom(c): c for c in self.chroms if f"sli.var/{c}/gnotate-variant.bin" in names}
        st = os.stat(self.path)
        self.stamp = f"{os.path.basename(self.path)}:{st.st_size}:{st.st_mtime_ns}"

    def _read_lines(self, name: str) -> List[str]:
        try:
            return [l for l in self._zip.read(name).decode().splitlines() if l.strip()]
        except KeyError:
            return []

    def _array(self, chrom_dir: str, name: str, typecode: str) -> array:
        a = array(typecode)
        a.frombytes(self._zip.read(f"sli.var/{chrom_dir}/{name}"))
        if sys.byteorder == "big":  # the files are little-endian
            a.byteswap()
        return a

    def has_chrom(self, chrom: str) -> bool:
        return canon_chrom(chrom) in self._dirs

    def lookup(
        self, chrom: str, positions: Sequence[int], refs: Sequence[str], alts: Sequence[str]
    ) -> Tuple[List[bool], Dict[str, List[float]]]:
        """Look up variants on one chromosome. ``positions`` are VCF POS (1-based) and the
        alleles are as written in the VCF. Returns (found flags, {field: values}); values are
        only meaningful where the flag is True."""
        n = len(positions)
        found = [False] * n
        values = {f: [0.0] * n for f in self.fields}
        d = self._dirs.get(canon_chrom(chrom))
        if d is None or n == 0:
            return found, values

        keyed: List[Tuple[int, int]] = []  # (row, packed key)
        long_rows: List[int] = []
        for i in range(n):
            ref, alt = refs[i].upper(), alts[i].upper()
            code = encode_allele(ref, alt)
            if code is None:
                if len(ref) + len(alt) > MAX_PACKED:
                    long_rows.append(i)
                continue
            keyed.append((i, ((positions[i] - 1) << 36) | code))

        if keyed:
            variant = self._array(d, "gnotate-variant.bin", "Q")
            size = len(variant)
            hits: List[Tuple[int, int]] = []
            for i, key in keyed:
                at = bisect.bisect_left(variant, key)
                if at < size and variant[at] == key:
                    hits.append((i, at))
            del variant
            if hits:
                for f in self.fields:
                    arr = self._array(d, f"gnotate-{f}.bin", "f")
                    column = values[f]
                    for i, at in hits:
                        column[i] = arr[at]
                    del arr
                for i, _ in hits:
                    found[i] = True

        if long_rows:
            wanted = {(positions[i] - 1, refs[i].upper(), alts[i].upper()): i for i in long_rows}
            with self._zip.open(f"sli.var/{d}/long-alleles.txt") as fh:
                next(fh, None)  # header
                for raw in fh:
                    p, r, a, _flt, vals = raw.decode().rstrip("\n").split("\t")
                    i = wanted.get((int(p), r, a))
                    if i is not None:
                        for f, v in zip(self.fields, vals.split("|")):
                            values[f][i] = float(v)
                        found[i] = True
        return found, values
