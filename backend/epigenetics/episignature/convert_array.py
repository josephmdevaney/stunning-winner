#!/usr/bin/env python3
"""Turn EpiSignaLR's array-derived reference data (a pickled dict of pandas tables) into plain JSON.

    python -I convert_array.py <no_strand_all_points_dict.pickle> <collapsed_coordinates_sorted.bed> <out.json>

The pickle holds one table per condition: a row per signature CpG (chromosome, 1-based position), a column per
reference profile (the median array methylation of each of 34 conditions and of controls), and a closing
`Start_forward` column. It is read the way convert.py reads a model: pandas is never imported, anything the file
asks for outside a short allow-list stops the conversion, and each pandas class is an empty stand-in that keeps
the state it is given. The tables are put back together from the stored blocks by hand.

The CpG order in each table is what pairs a genome's values with the reference rows, so it is checked against the
coordinates EpiSignaLR intersects a pileup with (`collapsed_coordinates_sorted.bed`): the two must list the same
CpGs in the same order, one position apart (1-based table, 0-based bed). Self-contained, like convert.py, so it
runs under `python -I`.
"""

from __future__ import annotations

import hashlib
import json
import math
import pickle
import sys
import warnings

import numpy as np

ALLOWED = {
    ("numpy", "dtype"),
    ("numpy", "ndarray"),
    ("numpy._core.multiarray", "_reconstruct"),
    ("numpy._core.multiarray", "scalar"),
    ("numpy.core.multiarray", "_reconstruct"),
    ("numpy.core.multiarray", "scalar"),
    ("builtins", "slice"),
}
EXPECTED_CONDITIONS = 34
LAST_COLUMN = "Start_forward"


class Stub:
    """Stands in for a pandas class: keeps its constructor arguments and any state it is handed."""

    qualified = ""

    def __init__(self, *args, **kwargs):
        self.args = args

    def __setstate__(self, state):
        self.state = state


def _unpickle_block(values, placement, ndim):
    """Stand-in for pandas' block constructor: just the stored values and which columns they belong to."""
    return {"values": values, "placement": placement}


class RestrictedUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if (module, name) in ALLOWED:
            return super().find_class(module, name)
        if (module, name) == ("pandas._libs.internals", "_unpickle_block"):
            return _unpickle_block
        if module.startswith("pandas."):
            return type(name, (Stub,), {"qualified": f"{module}.{name}"})
        raise pickle.UnpicklingError(f"refusing to load global {module}.{name}")


def _number(v):
    if isinstance(v, str):
        return None if v in ("", "NA", "nan") else float(v)
    f = float(v)
    return None if math.isnan(f) else f


def table(df) -> tuple[list[str], dict[str, list]]:
    """(column names, {column: values}) from a stand-in DataFrame."""
    blocks, axes = df.state["_mgr"].args
    columns = [str(c) for c in axes[0].args[1]["data"]]
    out: dict[str, list] = {}
    for block in blocks:
        placement = block["placement"]
        index = list(range(*placement.indices(len(columns)))) if isinstance(placement, slice) else list(placement)
        for k, j in enumerate(index):
            out[columns[j]] = list(block["values"][k])
    return columns, out


def read_bed(path: str) -> dict[str, list[tuple[str, int]]]:
    """condition -> its CpGs (chromosome, 0-based start) in file order; a row may list several conditions."""
    out: dict[str, list[tuple[str, int]]] = {}
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            chrom, start, _end, names = line.split()[:4]
            for name in names.split(","):
                out.setdefault(name.strip(), []).append((chrom, int(start)))
    return out


def convert(pickle_path: str, bed_path: str) -> dict:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with open(pickle_path, "rb") as f:
            frames = RestrictedUnpickler(f).load()
    if not isinstance(frames, dict) or len(frames) != EXPECTED_CONDITIONS:
        raise ValueError(f"expected a dict of {EXPECTED_CONDITIONS} conditions")
    bed = read_bed(bed_path)

    conditions = {}
    for key, df in frames.items():
        # EpiSignaLR renames this condition's key (not its column) before scoring; the bed already uses the new name.
        name = "MRXCJS" if key == "MRXSCJ" else key
        columns, cols = table(df)
        profiles = columns[2:-1]
        if columns[:2] != ["Chr", "Pos"] or columns[-1] != LAST_COLUMN or key not in profiles:
            raise ValueError(f"{key}: unexpected columns {columns[:3]}…{columns[-1]}")
        coords = [(str(c), int(p) - 1) for c, p in zip(cols["Chr"], cols["Pos"])]
        if coords != bed.get(name):
            raise ValueError(f"{key}: the reference rows are not the bed's CpGs in the bed's order")
        if [int(v) for v in cols[LAST_COLUMN]] != [s for _, s in coords]:
            raise ValueError(f"{key}: {LAST_COLUMN} disagrees with Pos")
        conditions[name] = {
            "positive_profile": key,                 # the reference column the SVM treats as the positive example
            "coords": [list(c) for c in coords],
            "profiles": profiles,
            "values": [[_number(cols[p][i]) for p in profiles] for i in range(len(coords))],
        }
    extra = set(bed) - set(conditions)
    if extra:
        raise ValueError(f"bed lists conditions with no reference table: {sorted(extra)}")
    with open(pickle_path, "rb") as f:
        digest = hashlib.sha256(f.read()).hexdigest()
    return {"format": 1, "source": {"file": pickle_path.rsplit("/", 1)[-1], "sha256": digest}, "conditions": conditions}


def main(argv: list[str]) -> int:
    if len(argv) != 4:
        print(__doc__.split("\n\n")[0], file=sys.stderr)
        return 2
    data = convert(argv[1], argv[2])
    with open(argv[3], "w") as f:
        json.dump(data, f, separators=(",", ":"))
    n_cpgs = len({tuple(c) for v in data["conditions"].values() for c in v["coords"]})
    print(f"{len(data['conditions'])} conditions, {n_cpgs} distinct CpGs -> {argv[3]}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
