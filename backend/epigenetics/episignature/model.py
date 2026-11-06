"""An EpiSignaLR episignature model: scoring one sample's DMR methylation, in plain Python.

The authors' pipeline is StandardScaler -> ElasticNet feature selection -> linear SVM with Platt-scaled
probabilities, trained on a handful of long-read genomes per condition (convert.py turns the pickle into JSON).
For a sample with methylation fractions x at the DMR panel:

    z_j  = (x_j - mean_j) / scale_j            the selected features only
    f    = sum_j w_j z_j + b                    the SVM decision value; each term is one DMR's contribution
    P    = Platt sigmoid of f                   the case probability the authors report

A DMR with no data is replaced by the training mean (z = 0), as EpiSignaLR does. A probability of at least
0.5 is a positive call and 0.4-0.6 is flagged indeterminate (both from the model's own metadata).
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from functools import lru_cache
from typing import Optional

from . import fetched

MIN_PROBABILITY = 1e-7      # libsvm clamps the Platt probability away from 0 and 1


def _sigmoid_predict(decision: float, a: float, b: float) -> float:
    fapb = decision * a + b
    # written to avoid cancelling 1 - p for large |fapb|, as libsvm does
    return math.exp(-fapb) / (1.0 + math.exp(-fapb)) if fapb >= 0 else 1.0 / (1.0 + math.exp(fapb))


def _multiclass_probability(r: list[list[float]]) -> list[float]:
    """libsvm's `multiclass_probability` (Wu, Lin and Weng's second method) from pairwise probabilities `r`.

    scikit-learn runs it even for two classes, where it stops at a tolerance of 0.005 / k and so returns
    probabilities within about 0.005 of the plain Platt sigmoid. It is ported as written so that scores match
    the ones the EpiSignaLR authors' scikit-learn pipeline reports.
    """
    k = len(r)
    p = [1.0 / k] * k
    q = [[0.0] * k for _ in range(k)]
    for t in range(k):
        for j in range(t):
            q[t][t] += r[j][t] * r[j][t]
            q[t][j] = q[j][t]
        for j in range(t + 1, k):
            q[t][t] += r[j][t] * r[j][t]
            q[t][j] = -r[j][t] * r[t][j]
    eps = 0.005 / k
    for _ in range(max(100, k)):
        qp = [sum(q[t][j] * p[j] for j in range(k)) for t in range(k)]
        pqp = sum(p[t] * qp[t] for t in range(k))
        if max(abs(qp[t] - pqp) for t in range(k)) < eps:
            break
        for t in range(k):
            diff = (-qp[t] + pqp) / q[t][t]
            p[t] += diff
            pqp = (pqp + diff * (diff * q[t][t] + 2 * qp[t])) / (1 + diff) / (1 + diff)
            for j in range(k):
                qp[j] = (qp[j] + diff * q[t][j]) / (1 + diff)
                p[j] /= 1 + diff
    return p


def case_probability(decision: float, prob_a: float, prob_b: float) -> float:
    """P(case) from the SVM decision value, as scikit-learn's SVC.predict_proba computes it.

    libsvm fits its sigmoid to the opposite sign of scikit-learn's public decision value for a binary problem,
    and the first class is control, so the case probability is the second of the two.
    """
    first = min(max(_sigmoid_predict(-decision, prob_a, prob_b), MIN_PROBABILITY), 1.0 - MIN_PROBABILITY)
    return _multiclass_probability([[0.0, first], [1.0 - first, 0.0]])[1]


@dataclass(frozen=True)
class Contribution:
    dmr: str
    weight: float
    value: Optional[float]        # the sample's fraction, or None when imputed
    z: float
    term: float                   # weight * z, this DMR's share of the decision value
    imputed: bool


@dataclass(frozen=True)
class Score:
    probability: float
    decision: float
    intercept: float
    contributions: list[Contribution]

    @property
    def imputed(self) -> list[str]:
        return [c.dmr for c in self.contributions if c.imputed]


@dataclass(frozen=True)
class EpiModel:
    key: str
    label: str
    features: list[str]           # every DMR the scaler was fit on, in column order
    mean: list[float]
    scale: list[float]
    support: list[bool]
    weights: list[float]          # one per selected feature
    intercept: float
    prob_a: float
    prob_b: float
    metadata: dict

    @property
    def selected(self) -> list[str]:
        return [f for f, keep in zip(self.features, self.support) if keep]

    @property
    def threshold(self) -> float:
        return float(self.metadata["scoring"]["threshold"])

    @property
    def indeterminate_zone(self) -> tuple[float, float]:
        lo, hi = self.metadata["scoring"]["indeterminate_zone"]
        return float(lo), float(hi)

    def score(self, values: dict[str, Optional[float]]) -> Score:
        """`values` maps DMR name -> methylation fraction (0-1), or None / missing when the DMR has no data."""
        contributions = []
        decision = self.intercept
        w = iter(self.weights)
        for j, name in enumerate(self.features):
            if not self.support[j]:
                continue
            weight = next(w)
            x = values.get(name)
            imputed = x is None or (isinstance(x, float) and math.isnan(x))
            z = 0.0 if imputed else (x - self.mean[j]) / self.scale[j]
            term = weight * z
            decision += term
            contributions.append(Contribution(name, weight, None if imputed else x, z, term, imputed))
        return Score(case_probability(decision, self.prob_a, self.prob_b), decision, self.intercept, contributions)


def _read_training(path) -> tuple[list[str], list[dict]]:
    """(DMR columns in file order, rows) from enet_svm_training_data.csv: DMR_* columns, group_bin, sample_id."""
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        columns = [c for c in (reader.fieldnames or []) if c not in ("sample_id", "group_bin")]
        rows = [{**{c: (float(r[c]) if r[c] not in ("", "NA", "nan") else None) for c in columns}, "group": int(float(r["group_bin"]))}
                for r in reader]
    return columns, rows


@lru_cache(maxsize=None)
def training(key: str) -> dict:
    """The model's training cohort as `{features, rows: [{group, DMR_x: fraction, ...}]}` (group 1 = case).
    Sample ids are dropped: only the distributions are used."""
    label = fetched.CATALOG[key][0]
    columns, rows = _read_training(fetched.model_files(label)["training"])
    return {"features": columns, "rows": rows}


@lru_cache(maxsize=None)
def loo_scores(key: str) -> list[dict]:
    """The authors' leave-one-out case probabilities for their own training samples: `{group, score}`."""
    path = fetched.model_files(fetched.CATALOG[key][0])["loo"]
    with open(path, newline="") as f:
        return [{"group": int(float(r["true_label"])), "score": float(r["loo_prob_case"])} for r in csv.DictReader(f)]


@lru_cache(maxsize=None)
def load(key: str) -> EpiModel:
    """The converted model joined to its feature names. Raises ValueError if the pieces do not agree."""
    label = fetched.CATALOG[key][0]
    files = fetched.model_files(label)
    converted = json.loads(files["converted"].read_text())
    metadata = json.loads(files["metadata"].read_text())
    features = training(key)["features"]

    if len(features) != converted["n_features"]:
        raise ValueError(f"{label}: training data has {len(features)} DMRs, the model was fit on {converted['n_features']}")
    selected = [f for f, keep in zip(features, converted["support"]) if keep]
    expected = metadata["final_model"]["selected_dmrs"]
    if selected != expected:
        raise ValueError(f"{label}: the model selects {selected}, its metadata lists {expected}")
    published = metadata["final_model"]["svm_weights"]
    if any(abs(w - published[name]) > 1e-6 for name, w in zip(selected, converted["weights"])):
        raise ValueError(f"{label}: converted SVM weights differ from the weights in the model metadata")

    return EpiModel(
        key=key, label=label, features=features, mean=converted["mean"], scale=converted["scale"],
        support=converted["support"], weights=converted["weights"], intercept=converted["intercept"],
        prob_a=converted["prob_a"], prob_b=converted["prob_b"], metadata=metadata,
    )
