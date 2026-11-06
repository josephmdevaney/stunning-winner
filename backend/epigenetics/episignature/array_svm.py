#!/usr/bin/env python3
"""Score one genome against EpiSignaLR's 34 array-derived classifiers.

    python -I array_svm.py <array_references.json>      # a JSON object {"chrom:start": fraction | null} on stdin

Writes {"scores": {condition: SVM decision value}, "sklearn": version} as JSON on stdout. Needs numpy and scikit-learn, which is why it is a
separate process from the API: it runs only the converted reference data (convert_array.py) and a genome's values.

It repeats what EpiSignaLR's `run_SVM.py` does for each condition, which is not a fixed model: the SVM is trained
anew for every genome. The condition's reference table (one median methylation profile per condition, plus
controls) gets the genome's values as an extra row; every missing value, in the reference or the genome, is
replaced by that CpG's mean over all rows; a linear SVC with the condition's own profile as the single positive
example, weighted 10 to 1, is fitted to the reference rows; and the genome's decision value is the score.
Because the genome's values enter the column means, they reach the fit through any reference gap, so the fit
cannot be done once and stored.
"""

from __future__ import annotations

import json
import sys
import warnings

import numpy as np
import sklearn
from sklearn.svm import SVC

POSITIVE_WEIGHT = 10


def score(references: dict, genome: dict) -> dict[str, float]:
    out = {}
    for name, cond in references["conditions"].items():
        reference = np.array([[np.nan if v is None else v for v in row] for row in cond["values"]], dtype=float).T
        sample = np.array([[np.nan if genome.get(f"{c}:{s}") is None else genome[f"{c}:{s}"] for c, s in cond["coords"]]], dtype=float)
        x = np.vstack([reference, sample])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fill = np.nanmean(x, axis=0)
        if np.isnan(fill).any():
            raise ValueError(f"{name}: a CpG has no value in any row")
        x = np.where(np.isnan(x), fill, x)
        y = [1 if p == cond["positive_profile"] else 0 for p in cond["profiles"]]
        svm = SVC(kernel="linear", class_weight={0: 1, 1: POSITIVE_WEIGHT})
        svm.fit(x[:-1], y)
        out[name] = float(svm.decision_function(x[-1:])[0])
    return out


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__.split("\n\n")[0], file=sys.stderr)
        return 2
    with open(argv[1]) as f:
        references = json.load(f)
    json.dump({"scores": score(references, json.load(sys.stdin)), "sklearn": sklearn.__version__}, sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
