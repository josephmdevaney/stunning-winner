"""A minimal reader for UCSC bigBed files: enough to dump every record, with no external tool.

Layout (little-endian files, as UCSC writes them): a 64-byte header, zoom headers, the autoSql
definition, a chromosome B+ tree (name -> id and size), then the data blocks indexed by an R-tree.
Each data block is zlib-compressed when the header's uncompress-buffer size is non-zero; its records are
``chromId, start, end`` (uint32 each) followed by the remaining fields as one tab-separated string.

This is not a general bigBed library (no range queries, no zoom levels, no big-endian files); it reads
the whole file, which is what building the MANE gene table needs."""

from __future__ import annotations

import struct
import zlib
from typing import Dict, Iterator, List, Tuple

BIGBED_MAGIC = 0x8789F2EB
BPT_MAGIC = 0x78CA8C91
CIRTREE_MAGIC = 0x2468ACE0


class BigBed:
    def __init__(self, data: bytes):
        self.b = data
        if struct.unpack_from("<I", data, 0)[0] != BIGBED_MAGIC:
            raise ValueError("not a little-endian bigBed file")
        (_ver, _zooms, self.chrom_tree, self.full_data, self.full_index, self.field_count, self.defined_fields,
         self.autosql_off, _summary, self.uncompress_buf) = struct.unpack_from("<HHQQQHHQQI", data, 4)
        self.autosql = ""
        if self.autosql_off:
            self.autosql = data[self.autosql_off:data.index(b"\0", self.autosql_off)].decode()
        self.chroms = self._chromosomes()

    def field_names(self) -> List[str]:
        """Column names from the autoSql definition, in order."""
        names = []
        for line in self.autosql.splitlines()[3:]:
            line = line.strip()
            if line.startswith(")") or not line.endswith(";") and ";" not in line:
                continue
            decl = line.split(";")[0].split()
            if len(decl) >= 2:
                names.append(decl[-1])
        return names

    def _chromosomes(self) -> Dict[int, Tuple[str, int]]:
        b, off = self.b, self.chrom_tree
        magic, _block, key_size, val_size, _count, _res = struct.unpack_from("<IIIIQQ", b, off)
        if magic != BPT_MAGIC:
            raise ValueError("bad chromosome tree")
        out: Dict[int, Tuple[str, int]] = {}

        def node(at: int) -> None:
            is_leaf, _r, count = struct.unpack_from("<BBH", b, at)
            at += 4
            for _ in range(count):
                key = b[at:at + key_size].split(b"\0", 1)[0].decode()
                at += key_size
                if is_leaf:
                    chrom_id, size = struct.unpack_from("<II", b, at)
                    out[chrom_id] = (key, size)
                else:
                    node(struct.unpack_from("<Q", b, at)[0])
                at += val_size if is_leaf else 8

        node(off + 32)
        return out

    def _leaves(self) -> Iterator[Tuple[int, int]]:
        b = self.b
        if struct.unpack_from("<I", b, self.full_index)[0] != CIRTREE_MAGIC:
            raise ValueError("bad R-tree index")

        def node(at: int) -> Iterator[Tuple[int, int]]:
            is_leaf, _r, count = struct.unpack_from("<BBH", b, at)
            at += 4
            for _ in range(count):
                if is_leaf:
                    _sc, _sb, _ec, _eb, data_off, data_size = struct.unpack_from("<IIIIQQ", b, at)
                    at += 32
                    yield data_off, data_size
                else:
                    child = struct.unpack_from("<Q", b, at + 16)[0]
                    at += 24
                    yield from node(child)

        yield from node(self.full_index + 48)

    def records(self) -> Iterator[Tuple[str, int, int, List[str]]]:
        """(chrom, start, end, [remaining fields]) for every record, in index order."""
        b = self.b
        for off, size in self._leaves():
            block = b[off:off + size]
            if self.uncompress_buf:
                block = zlib.decompress(block)
            at = 0
            while at < len(block):
                chrom_id, start, end = struct.unpack_from("<III", block, at)
                at += 12
                stop = block.index(b"\0", at)
                rest = block[at:stop].decode().split("\t")
                at = stop + 1
                yield self.chroms[chrom_id][0], start, end, rest
