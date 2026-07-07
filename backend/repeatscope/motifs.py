"""Shared helpers for comparing STR motif sequences."""

from __future__ import annotations


def is_cyclic_rotation(a: str, b: str) -> bool:
    """True if `a` is a cyclic rotation of `b` (e.g. AAG/GAA/AGA are all rotations of
    each other) -- STR motif notation's reading frame is arbitrary, so two sources
    (e.g. STRchive and a VCF's own MOTIFS field) may record the "same" repeat unit
    starting at a different base."""
    return len(a) == len(b) and a in (b + b)
