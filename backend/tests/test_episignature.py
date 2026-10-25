"""EpiSignaLR episignature scoring (epigenetics/episignature/) and its route.

The maths and the region averaging are tested on synthetic models and CpGs. Tests that need the fetched
EpiSignaLR files (backend/scripts/setup_episignalr.py) skip when they are not there; those use the files'
own training rows and probabilities recorded from the authors' scikit-learn pipeline (1.8.0)."""

import os
import pickle
import subprocess
import sys

import pytest
from fastapi.testclient import TestClient

from epigenetics import cpg, router, samples
from epigenetics.cpg import CpG
from epigenetics.episignature import model, regions, service, vendor

needs_vendor = pytest.mark.skipif(not vendor.installed(), reason="EpiSignaLR files not fetched (backend/scripts/setup_episignalr.py)")

KABUKI_A, KABUKI_B = -2.852093863675567, 0.25868962896184056


# ---- the probability, against scikit-learn ----------------------------------------------------------------------

# (decision value, SVC.predict_proba[:, 1]) from scikit-learn 1.8.0 for the Kabuki model's sigmoid.
SKLEARN_PROBABILITY = [
    (-3.0, 0.000249115909), (-1.0, 0.069564149419), (-0.0973, 0.5),
    (0.0, 0.564463863166), (0.5, 0.843634529522), (1.5, 0.989407171832),
]


@pytest.mark.parametrize("decision, expected", SKLEARN_PROBABILITY)
def test_case_probability_reproduces_scikit_learns(decision, expected):
    assert model.case_probability(decision, KABUKI_A, KABUKI_B) == pytest.approx(expected, abs=1e-9)


def test_case_probability_is_increasing_and_stays_off_zero_and_one():
    values = [model.case_probability(f / 4, KABUKI_A, KABUKI_B) for f in range(-60, 61)]
    assert values == sorted(values)
    assert 0 < values[0] and values[-1] < 1


def test_two_class_refinement_stays_within_its_tolerance_of_the_plain_sigmoid():
    # libsvm stops at 0.005 / k, so the refined value is within about 0.005 of 1 - sigmoid
    for first in (0.05, 0.3, 0.5, 0.8, 0.97):
        p = model._multiclass_probability([[0.0, first], [1 - first, 0.0]])
        assert sum(p) == pytest.approx(1.0) and abs(p[0] - first) < 0.005


# ---- scoring a sample -------------------------------------------------------------------------------------------

def _model(**over):
    fields = dict(
        key="syn", label="Syn", features=["DMR_1", "DMR_2", "DMR_3"], mean=[0.5, 0.4, 0.2], scale=[0.1, 0.2, 0.05],
        support=[True, False, True], weights=[-0.5, 0.25], intercept=-0.1, prob_a=-2.0, prob_b=0.2,
        metadata={"scoring": {"threshold": 0.5, "indeterminate_zone": [0.4, 0.6]}},
    )
    fields.update(over)
    return model.EpiModel(**fields)


def test_selected_features_follow_the_support_mask():
    assert _model().selected == ["DMR_1", "DMR_3"]
    assert _model().threshold == 0.5 and _model().indeterminate_zone == (0.4, 0.6)


def test_decision_is_the_intercept_plus_each_selected_dmrs_contribution():
    s = _model().score({"DMR_1": 0.6, "DMR_2": 0.9, "DMR_3": 0.3})       # DMR_2 is not selected: ignored
    assert [c.dmr for c in s.contributions] == ["DMR_1", "DMR_3"]
    z1, z3 = (0.6 - 0.5) / 0.1, (0.3 - 0.2) / 0.05
    assert s.contributions[0].z == pytest.approx(z1) and s.contributions[1].z == pytest.approx(z3)
    assert s.contributions[0].term == pytest.approx(-0.5 * z1)
    assert s.decision == pytest.approx(-0.1 + -0.5 * z1 + 0.25 * z3)
    assert s.decision == pytest.approx(s.intercept + sum(c.term for c in s.contributions))
    assert s.probability == pytest.approx(model.case_probability(s.decision, -2.0, 0.2))


def test_a_dmr_with_no_data_is_the_training_mean_and_contributes_nothing():
    s = _model().score({"DMR_1": 0.6, "DMR_3": None})
    assert s.imputed == ["DMR_3"] and s.contributions[1].z == 0 and s.contributions[1].term == 0
    assert s.contributions[1].value is None
    same = _model().score({"DMR_1": 0.6, "DMR_3": 0.2})                  # 0.2 is that DMR's training mean
    assert s.decision == pytest.approx(same.decision) and same.imputed == []
    assert _model().score({"DMR_1": float("nan"), "DMR_3": 0.2}).imputed == ["DMR_1"]
    assert _model().score({}).decision == pytest.approx(-0.1)            # nothing known: only the intercept


# ---- region averaging -------------------------------------------------------------------------------------------

def test_parse_panel_groups_cpgs_by_dmr_in_file_order(tmp_path):
    bed = tmp_path / "p.bed"
    bed.write_text("chr1\t100\t101\tDMR_B\nchr1\t120\t121\tDMR_B\nchr2\t5\t6\tDMR_A\n")
    panel = regions.parse_panel(bed)
    assert list(panel) == ["DMR_B", "DMR_A"]
    assert panel["DMR_B"].sites == {"chr1": (100, 120)} and panel["DMR_B"].n_listed == 2
    assert panel["DMR_B"].span == ("chr1", 100, 121)


def _cpg(pos, cov, mod):
    return CpG(pos, 100.0 * mod / cov, cov, mod, cov - mod)


def test_region_value_is_total_modified_over_total_coverage_of_the_listed_cpgs(monkeypatch, tmp_path):
    bed = tmp_path / "p.bed"
    bed.write_text("chr1\t100\t101\tD\nchr1\t110\t111\tD\nchr1\t120\t121\tD\nchr1\t130\t131\tD\n")
    panel = regions.parse_panel(bed)
    pileup = [_cpg(100, 10, 10), _cpg(105, 99, 0),           # 105 is not listed: ignored
              _cpg(110, 30, 0), _cpg(120, 1, 1)]             # 120 has a single read: below the minimum
    monkeypatch.setattr(cpg, "fetch", lambda b, chrom, s, e: [c for c in pileup if s <= c.pos < e])
    v = regions.region_values(tmp_path / "x.bed.gz", panel)["D"]
    assert v.n_cpgs == 2 and v.mean_coverage == 20.0
    assert v.fraction == pytest.approx(10 / 40)                # pooled, so the deep CpG counts for more than a plain mean (0.5)


def test_a_dmr_without_a_covered_cpg_has_no_value(monkeypatch, tmp_path):
    bed = tmp_path / "p.bed"
    bed.write_text("chr1\t100\t101\tD\n")
    monkeypatch.setattr(cpg, "fetch", lambda *a: [])
    v = regions.region_values(tmp_path / "x.bed.gz", regions.parse_panel(bed))["D"]
    assert v.fraction is None and v.n_cpgs == 0


# ---- the converter: reads a pickle without trusting it ----------------------------------------------------------

def _python_with_numpy_and_joblib():
    for py in (sys.executable, str(vendor.VENDOR_DIR / ".convert-venv" / "bin" / "python")):
        if os.path.exists(py) and subprocess.run([py, "-I", "-c", "import numpy, joblib"], capture_output=True).returncode == 0:
            return py
    return None


needs_converter_python = pytest.mark.skipif(_python_with_numpy_and_joblib() is None, reason="no Python with numpy and joblib")
CONVERT = vendor.VENDOR_DIR.parent / "convert.py"

WRITE_FAKE_MODEL = '''
import sys, types, numpy as np, joblib

def fake(module, name, **state):
    parts = module.split(".")
    for i in range(1, len(parts) + 1):
        sys.modules.setdefault(".".join(parts[:i]), types.ModuleType(".".join(parts[:i])))
    cls = type(name, (), {})
    cls.__module__ = module
    setattr(sys.modules[module], name, cls)
    obj = cls()
    obj.__dict__.update(state)
    return obj

scaler = fake("sklearn.preprocessing._data", "StandardScaler", mean_=np.array([.5, .4, .2]), scale_=np.array([.1, .2, .05]))
select = fake("__main__", "ElasticNetSelector", support_=np.array([True, False, True]))
svm = fake("sklearn.svm._classes", "SVC", kernel="linear", probability=True, _sparse=False, C=0.01,
           classes_=np.array([0, 1]), support_vectors_=np.array([[1., 2.], [3., -1.], [0., 4.]]),
           dual_coef_=np.array([[0.5, -0.25, 0.125]]), intercept_=np.array([-0.3]),
           _probA=np.array([-2.5]), _probB=np.array([0.1]))
pipe = fake("sklearn.pipeline", "Pipeline", steps=[("scaler", scaler), ("select", select), ("svm", svm)])
joblib.dump(pipe, sys.argv[1])
'''


@needs_converter_python
def test_converter_turns_a_pipeline_into_scoring_numbers_without_scikit_learn(tmp_path):
    py = _python_with_numpy_and_joblib()
    pkl, out = tmp_path / "m.joblib", tmp_path / "m.json"
    subprocess.run([py, "-c", WRITE_FAKE_MODEL, str(pkl)], check=True, cwd=tmp_path)
    done = subprocess.run([py, "-I", str(CONVERT), str(pkl), str(out)], capture_output=True, text=True, cwd=tmp_path)
    assert done.returncode == 0, done.stderr
    import json
    m = json.loads(out.read_text())
    assert m["n_features"] == 3 and m["support"] == [True, False, True]
    assert m["mean"] == [0.5, 0.4, 0.2] and m["scale"] == [0.1, 0.2, 0.05]
    # weights = dual_coef_ @ support_vectors_ = [.5*1 - .25*3 + 0, .5*2 + .25*1 + .125*4]
    assert m["weights"] == pytest.approx([-0.25, 1.75])
    assert m["intercept"] == -0.3 and m["prob_a"] == -2.5 and m["prob_b"] == 0.1


class _Evil:
    def __reduce__(self):
        return (os.system, ("echo pwned > pwned.txt",))


@needs_converter_python
def test_converter_refuses_a_pickle_that_asks_for_anything_else_and_runs_nothing(tmp_path):
    py = _python_with_numpy_and_joblib()
    (tmp_path / "evil.joblib").write_bytes(pickle.dumps(_Evil()))
    done = subprocess.run([py, "-I", str(CONVERT), str(tmp_path / "evil.joblib"), str(tmp_path / "out.json")],
                          capture_output=True, text=True, cwd=tmp_path)
    assert done.returncode != 0 and "refusing to load global" in done.stderr
    assert not (tmp_path / "pwned.txt").exists() and not (tmp_path / "out.json").exists()


# ---- the fetched models -----------------------------------------------------------------------------------------

# Probability the authors' scikit-learn pipeline gives the first and last row of each model's own training data.
SKLEARN_TRAINING_ROWS = {
    "kabuki": (0.7556585488, 0.0725952239),
    "kleefstra": (0.9228247966, 0.116368815),
    "setd1b": (0.9121430531, 0.0787635412),
}


@needs_vendor
@pytest.mark.parametrize("key", list(vendor.CATALOG))
def test_fetched_model_scores_its_training_rows_as_scikit_learn_does(key):
    m = model.load(key)
    rows = model.training(key)["rows"]
    for row, expected in zip((rows[0], rows[-1]), SKLEARN_TRAINING_ROWS[key]):
        assert m.score({d: row[d] for d in m.features}).probability == pytest.approx(expected, abs=1e-8)


@needs_vendor
@pytest.mark.parametrize("key", list(vendor.CATALOG))
def test_fetched_model_is_consistent_with_its_panel_and_metadata(key):
    m = model.load(key)
    panel = regions.parse_panel(vendor.model_files(vendor.CATALOG[key][0])["dmrs"])
    assert set(m.features) == set(panel) and len(m.selected) == m.metadata["final_model"]["n_features"]
    assert len(model.loo_scores(key)) == m.metadata["training_cohort"]["n_samples"]
    assert m.threshold == 0.5 and m.indeterminate_zone == (0.4, 0.6)


# ---- the service ------------------------------------------------------------------------------------------------

def _fake_pileup(monkeypatch, fraction=0.5, mode="count", sites="reference"):
    def fake_fetch(bed, chrom, start, end):
        mod = round(fraction * 30)
        return [CpG(p, 100.0 * mod / 30, 30, mod, 30 - mod) for p in range(start, end)]
    monkeypatch.setattr(cpg, "fetch", fake_fetch)
    monkeypatch.setattr(cpg, "header", lambda bed: {"pb_cpg_tools_version": "3.0.0", "pileup_mode": mode, "modsites_mode": sites, "min_coverage": 4})


@pytest.fixture
def pileup_sample(tmp_path):
    bed = tmp_path / "S.combined.bed.gz"
    bed.write_bytes(b"")
    service._cache.clear()
    return samples.Sample(key="SYN", label="Synthetic", paths={"combined": bed})


@needs_vendor
def test_service_reports_each_model_with_its_evidence(monkeypatch, pileup_sample):
    _fake_pileup(monkeypatch)
    out = service.overview(pileup_sample)
    assert [m["id"] for m in out["models"]] == list(vendor.CATALOG)
    assert out["source"]["commit"] and out["provenance"]["pileup_mode"] == "count"
    kabuki = next(m for m in out["models"] if m["id"] == "kabuki")
    ref = model.load("kabuki")
    assert kabuki["probability"] == pytest.approx(ref.score({d: 0.5 for d in ref.features}).probability)
    assert kabuki["dmrs"]["selected"] == 13 and kabuki["dmrs"]["imputed_selected"] == 0
    f = kabuki["features"][0]
    assert f["value"] == 0.5 and f["controls"] and f["cases"] and f["region"]["chrom"].startswith("chr")
    assert sum(x["contribution"] for x in kabuki["features"]) + kabuki["intercept"] == pytest.approx(kabuki["decision"])
    assert kabuki["training"] == {"cases": 9, "controls": 10} and kabuki["performance"]["specificity"] == 0.5
    assert kabuki["warnings"] == []                      # count/reference pileup, nothing imputed


@needs_vendor
def test_a_pileup_made_in_another_mode_is_flagged_on_every_model(monkeypatch, pileup_sample):
    _fake_pileup(monkeypatch, mode="model", sites="denovo")
    for m in service.overview(pileup_sample)["models"]:
        assert [w["id"] for w in m["warnings"]] == ["pileup_mode"]
        assert "model mode with denovo sites" in m["warnings"][0]["text"]


@needs_vendor
def test_no_call_is_made_when_too_many_selected_dmrs_have_no_data(monkeypatch, pileup_sample):
    _fake_pileup(monkeypatch)
    monkeypatch.setattr(cpg, "fetch", lambda *a: [])
    kabuki = next(m for m in service.overview(pileup_sample)["models"] if m["id"] == "kabuki")
    assert kabuki["probability"] is None and kabuki["call"] == "insufficient_data"
    assert kabuki["dmrs"]["imputed_selected"] == 13
    assert any(w["id"] == "imputed" and w["severity"] == "caution" for w in kabuki["warnings"])


def test_call_bands_follow_the_models_threshold_and_indeterminate_zone():
    m = _model()
    assert [service._band(p, m) for p in (0.1, 0.39, 0.4, 0.55, 0.6, 0.61, 0.99)] == \
        ["negative", "negative", "indeterminate", "indeterminate", "indeterminate", "positive", "positive"]


def test_unrecorded_pileup_mode_is_a_caution():
    ids = lambda prov: [w["id"] for w in service._provenance_warnings(prov)]
    assert ids({"pileup_mode": None, "modsites_mode": None}) == ["pileup_unrecorded"]
    assert ids({"pileup_mode": "count", "modsites_mode": "reference"}) == []
    assert ids({"pileup_mode": "count", "modsites_mode": "denovo"}) == ["pileup_mode"]


# ---- route ------------------------------------------------------------------------------------------------------

@pytest.fixture
def client(tmp_path, monkeypatch):
    for name in ("S.combined.bed.gz", "S.combined.bed.gz.tbi"):
        (tmp_path / name).write_bytes(b"")
    cfg = tmp_path / "cfg.json"
    cfg.write_text('{"samples": {"SYN": {"dir": "%s", "prefix": "S", "label": "Synthetic"}}}' % tmp_path)
    monkeypatch.setenv("EPIGENETICS_CONFIG", str(cfg))
    from repeatscope.api.main import app
    with TestClient(app) as c:
        yield c


def test_route_says_how_to_install_when_the_models_are_not_fetched(client, monkeypatch):
    monkeypatch.setattr(router.episignature_vendor, "installed", lambda: False)
    body = client.get("/api/epigenetics/episignatures").json()
    assert body["installed"] is False and body["models"] == []
    assert body["setup_command"] == "python backend/scripts/setup_episignalr.py"
    assert {a["id"] for a in body["available"]} == {"setd1b", "kleefstra", "kabuki"}


def test_route_returns_the_service_result(client, monkeypatch):
    monkeypatch.setattr(router.episignature_vendor, "installed", lambda: True)
    monkeypatch.setattr(router.episignature_service, "overview", lambda s: {"sample": s.key, "models": [{"id": "kabuki"}]})
    body = client.get("/api/epigenetics/episignatures", params={"sample": "SYN"}).json()
    assert body == {"installed": True, "sample": "SYN", "models": [{"id": "kabuki"}]}
    assert client.get("/api/epigenetics/episignatures", params={"sample": "nope"}).status_code == 404


def test_route_reports_a_missing_tabix_and_inconsistent_files_as_unavailable(client, monkeypatch):
    monkeypatch.setattr(router.episignature_vendor, "installed", lambda: True)

    def broken(error):
        def overview(s):
            raise error
        return overview

    monkeypatch.setattr(router.episignature_service, "overview", broken(cpg.TabixError("tabix not found")))
    r = client.get("/api/epigenetics/episignatures")
    assert r.status_code == 503 and "tabix" in r.json()["detail"]
    monkeypatch.setattr(router.episignature_service, "overview", broken(ValueError("weights differ")))
    r = client.get("/api/epigenetics/episignatures")
    assert r.status_code == 503 and "inconsistent" in r.json()["detail"]
