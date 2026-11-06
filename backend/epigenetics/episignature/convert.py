#!/usr/bin/env python3
"""Turn an EpiSignaLR model (a joblib pickle of a scikit-learn pipeline) into plain JSON.

    python -I convert.py <enet_svm_model.joblib> <model.json>

A pickle can run arbitrary code when it is loaded, so this never loads it the usual way and never imports
scikit-learn. Every global the file asks for is checked against a short allow-list: numpy's array
reconstruction and joblib's array wrapper are real, the scikit-learn classes (and the authors' own
`ElasticNetSelector`) become empty stand-ins that only keep the state the file hands them, and anything else
-- `os.system`, `builtins.eval`, any other module -- stops the conversion. Only numpy and joblib are needed.

The pipeline is StandardScaler -> ElasticNetSelector -> SVC(kernel="linear", probability=True). With a linear
kernel the SVM is a weight vector over the standardised selected features plus an intercept, and the
probability is a Platt sigmoid of its decision value, so what is written out is everything scoring needs and
the API process scores in plain Python.

This file is deliberately self-contained (it runs under `python -I`, which does not put its own folder on the
import path), and it reads the pickle but writes nothing but the JSON.
"""

from __future__ import annotations

import hashlib
import json
import pickle
import sys

import numpy as np
from joblib.numpy_pickle import NumpyUnpickler

ALLOWED = {
    ("numpy", "dtype"),
    ("numpy", "ndarray"),
    ("numpy._core.multiarray", "_reconstruct"),
    ("numpy._core.multiarray", "scalar"),
    ("numpy.core.multiarray", "_reconstruct"),
    ("numpy.core.multiarray", "scalar"),
    ("joblib.numpy_pickle", "NumpyArrayWrapper"),
}
STUBBED_PREFIXES = ("sklearn.", "__main__")


class Stub:
    """Stands in for a class defined elsewhere: keeps only the state the pickle assigns it."""

    qualified = ""

    def __init__(self, *args, **kwargs):
        pass

    def __setstate__(self, state):
        self.state = state


class RestrictedUnpickler(NumpyUnpickler):
    def find_class(self, module, name):
        if (module, name) in ALLOWED:
            return super().find_class(module, name)
        if module == "__main__" or module.startswith(STUBBED_PREFIXES):
            return type(name, (Stub,), {"qualified": f"{module}.{name}"})
        raise pickle.UnpicklingError(f"refusing to load global {module}.{name}")


def load(path: str):
    with open(path, "rb") as f:
        return RestrictedUnpickler(path, f, True).load()


def _step(pipeline, name: str, expected: str) -> dict:
    steps = dict(pipeline.state["steps"])
    if name not in steps:
        raise ValueError(f"pipeline has no {name!r} step (found {sorted(steps)})")
    obj = steps[name]
    if not obj.qualified.endswith(expected):
        raise ValueError(f"step {name!r} is {obj.qualified}, expected …{expected}")
    return obj.state


def convert(path: str) -> dict:
    pipe = load(path)
    if not isinstance(pipe, Stub) or not pipe.qualified.endswith("Pipeline"):
        raise ValueError(f"{path} does not hold a scikit-learn Pipeline")
    scaler = _step(pipe, "scaler", "StandardScaler")
    select = _step(pipe, "select", "ElasticNetSelector")
    svm = _step(pipe, "svm", "SVC")

    if svm["kernel"] != "linear" or not svm["probability"] or svm.get("_sparse"):
        raise ValueError("only a dense linear SVC with probability=True is supported")
    if list(svm["classes_"]) != [0, 1]:
        raise ValueError(f"expected classes [0, 1], found {list(svm['classes_'])}")

    support = [bool(v) for v in select["support_"]]
    dual = np.asarray(svm["dual_coef_"], dtype=float)                 # (1, n_support_vectors)
    vectors = np.asarray(svm["support_vectors_"], dtype=float)        # (n_support_vectors, n_selected)
    n_selected = sum(support)
    if vectors.shape[1] != n_selected or dual.shape[0] != 1:
        raise ValueError("support vectors do not match the selected features")
    if len(scaler["mean_"]) != len(support):
        raise ValueError("scaler and selector disagree on the number of features")

    with open(path, "rb") as f:
        digest = hashlib.sha256(f.read()).hexdigest()
    return {
        "format": 1,
        "source": {"file": path.rsplit("/", 1)[-1], "sha256": digest},
        "n_features": len(support),
        "mean": [float(v) for v in scaler["mean_"]],
        "scale": [float(v) for v in scaler["scale_"]],
        "support": support,
        # f(x) = weights . z + intercept, with z the standardised selected features; sklearn's public sign.
        "weights": [float(v) for v in (dual @ vectors)[0]],
        "intercept": float(svm["intercept_"][0]),
        "prob_a": float(svm["_probA"][0]),
        "prob_b": float(svm["_probB"][0]),
        "svm": {"kernel": "linear", "C": float(svm["C"]), "n_support_vectors": int(vectors.shape[0])},
    }


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(__doc__.split("\n\n")[0], file=sys.stderr)
        return 2
    model = convert(argv[1])
    with open(argv[2], "w") as f:
        json.dump(model, f, indent=1)
        f.write("\n")
    print(f"{argv[1]}: {model['n_features']} features, {sum(model['support'])} selected -> {argv[2]}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
