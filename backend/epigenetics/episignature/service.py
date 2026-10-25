"""Scoring a registered genome against the EpiSignaLR episignature models.

For each model: average the genome's methylation over the model's DMR panel (regions.py), score the selected
DMRs (model.py), and say how far to trust the answer. The result carries the evidence behind the score -- each
selected DMR's value beside the training cases and controls, its share of the decision value -- and the
caveats that apply to this genome's data, so a probability is never shown without them.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Optional

from .. import cpg
from ..samples import Sample
from . import model as models
from . import fetched, regions

# EpiSignaLR derives the models from count-mode pileups at reference CpG sites.
EXPECTED_PILEUP_MODE = "count"
EXPECTED_MODSITES_MODE = "reference"
# Imputing more than this share of the selected DMRs means the score mostly reflects the training mean.
MAX_IMPUTED_SHARE = 0.25

_cache: dict[tuple, dict] = {}


@lru_cache(maxsize=None)
def _panel(key: str) -> dict[str, regions.Dmr]:
    return regions.parse_panel(fetched.model_files(fetched.CATALOG[key][0])["dmrs"])


def _provenance_warnings(provenance: dict) -> list[dict]:
    out = []
    mode, sites = provenance.get("pileup_mode"), provenance.get("modsites_mode")
    if mode is None and sites is None:
        out.append({"id": "pileup_unrecorded", "severity": "caution",
                    "text": "The pileup does not record how it was made. EpiSignaLR's scores come from count-mode pileups at reference "
                            "CpG sites; a different mode shifts methylation values, and the scores with them."})
    elif mode != EXPECTED_PILEUP_MODE or sites != EXPECTED_MODSITES_MODE:
        out.append({"id": "pileup_mode", "severity": "caution",
                    "text": f"This pileup was made in {mode or 'unknown'} mode with {sites or 'unknown'} sites; EpiSignaLR uses "
                            f"{EXPECTED_PILEUP_MODE} mode with {EXPECTED_MODSITES_MODE} sites. Methylation values differ between the modes, "
                            "so these scores are not directly comparable with the published ones. Re-run pb-CpG-tools with "
                            f"--pileup-mode {EXPECTED_PILEUP_MODE} --modsites-mode {EXPECTED_MODSITES_MODE} for a like-for-like result."})
    return out


def _band(score: float, m: models.EpiModel) -> str:
    lo, hi = m.indeterminate_zone
    if lo <= score <= hi:
        return "indeterminate"
    return "positive" if score >= m.threshold else "negative"


def _one(key: str, combined, provenance: dict) -> dict:
    m = models.load(key)
    panel = _panel(key)
    values = regions.region_values(combined, panel)
    score = m.score({name: v.fraction for name, v in values.items()})

    cohort = models.training(key)
    features = []
    for c in score.contributions:
        dmr = panel.get(c.dmr)
        v = values.get(c.dmr)
        controls = [r[c.dmr] for r in cohort["rows"] if r["group"] == 0 and r[c.dmr] is not None]
        cases = [r[c.dmr] for r in cohort["rows"] if r["group"] == 1 and r[c.dmr] is not None]
        j = m.features.index(c.dmr)
        features.append({
            "dmr": c.dmr,
            "region": None if dmr is None else dict(zip(("chrom", "start", "end"), dmr.span)),
            "weight": c.weight, "z": c.z, "contribution": c.term,
            "value": c.value, "imputed": c.imputed,
            "n_cpgs": v.n_cpgs if v else 0, "n_cpgs_listed": dmr.n_listed if dmr else 0,
            "mean_coverage": v.mean_coverage if v else 0.0,
            "training_mean": m.mean[j],
            "controls": controls, "cases": cases,
        })

    imputed = score.imputed
    withheld = len(imputed) / len(features) > MAX_IMPUTED_SHARE
    band = "insufficient_data" if withheld else _band(score.probability, m)
    warnings = _provenance_warnings(provenance)
    if imputed:
        warnings.append({"id": "imputed", "severity": "caution" if withheld else "note",
                         "text": f"{len(imputed)} of {len(features)} selected DMRs had no covered CpG and were replaced by the training mean "
                                 f"({', '.join(imputed)})." + (" Too many: no call is made." if withheld else "")})

    perf = m.metadata["loo_performance"]
    cohort_meta = m.metadata["training_cohort"]
    lo, hi = m.indeterminate_zone
    return {
        "id": key, "name": fetched.CATALOG[key][1], "gene": fetched.CATALOG[key][2],
        "probability": None if withheld else score.probability,
        "decision": score.decision, "intercept": score.intercept,
        "call": band,
        "threshold": m.threshold, "indeterminate_zone": [lo, hi],
        "features": features,
        "dmrs": {"panel": len(m.features), "selected": len(features), "with_data": sum(1 for v in values.values() if v.fraction is not None),
                 "imputed_selected": len(imputed)},
        "training": {"cases": cohort_meta["n_cases"], "controls": cohort_meta["n_controls"]},
        "performance": {
            "design": m.metadata["evaluation"],
            "auc": perf["auc"], "sensitivity": perf["sensitivity"], "specificity": perf["specificity"],
            "balanced_accuracy": perf["balanced_accuracy"], "permutation_p": perf["permutation_p"],
            "confusion_matrix": perf["confusion_matrix"],
        },
        "loo": models.loo_scores(key),
        "warnings": warnings,
    }


def overview(sample: Sample) -> dict:
    combined = sample.paths["combined"]
    stamp = (sample.key, str(combined), combined.stat().st_mtime_ns)
    if stamp in _cache:
        return _cache[stamp]
    provenance = cpg.header(combined)
    out = {
        "sample": sample.key, "label": sample.label,
        "source": {k: fetched.info()[k] for k in ("repo", "commit", "fetched_at", "source", "licence")},
        "provenance": provenance,
        "models": [_one(key, combined, provenance) for key in fetched.CATALOG],
    }
    _cache[stamp] = out
    return out
