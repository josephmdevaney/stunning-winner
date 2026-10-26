"""EpiSignaLR episignature scoring (epigenetics/episignature/) and its route.

The maths and the region averaging are tested on synthetic models and CpGs. Tests that need the fetched
EpiSignaLR files (backend/scripts/setup_episignalr.py) skip when they are not there; those use the files'
own training rows and probabilities recorded from the authors' scikit-learn pipeline (1.8.0)."""

import json
import os
import pickle
import subprocess
import sys

import pytest
from fastapi.testclient import TestClient

from epigenetics import cpg, router, samples
from epigenetics.cpg import CpG
from epigenetics.episignature import arrays, fetched, model, regions, service

needs_vendor = pytest.mark.skipif(not fetched.installed(), reason="EpiSignaLR files not fetched (backend/scripts/setup_episignalr.py)")

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

def _python_with(*modules):
    """An interpreter that can import `modules`: this one, or the one setup made."""
    for py in (sys.executable, str(fetched.VENDOR_DIR / ".venv" / "bin" / "python")):
        if os.path.exists(py) and subprocess.run([py, "-I", "-c", "import " + ", ".join(modules)], capture_output=True).returncode == 0:
            return py
    return None


def _python_with_numpy_and_joblib():
    return _python_with("numpy", "joblib")


needs_converter_python = pytest.mark.skipif(_python_with_numpy_and_joblib() is None, reason="no Python with numpy and joblib")
CONVERT = fetched.VENDOR_DIR.parent / "convert.py"

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


# ---- the array-derived classifiers: reading the pickle without trusting it -----------------------------------------

CONVERT_ARRAY = fetched.VENDOR_DIR.parent / "convert_array.py"
ARRAY_SVM = fetched.VENDOR_DIR.parent / "array_svm.py"
needs_numpy = pytest.mark.skipif(_python_with("numpy") is None, reason="no Python with numpy")
needs_scikit_learn = pytest.mark.skipif(_python_with("numpy", "sklearn") is None, reason="no Python with numpy and scikit-learn")

# Two conditions in a pickled stand-in for the pandas layout: an object block holding strings (with an "NA"), an
# integer block and a float block, placed by index array and by slice, as pandas stores them.
WRITE_FAKE_FRAMES = '''
import pickle, sys, types
import numpy as np

def put(module, name, obj):
    parts = module.split(".")
    for i in range(1, len(parts) + 1):
        sys.modules.setdefault(".".join(parts[:i]), types.ModuleType(".".join(parts[:i])))
    setattr(sys.modules[module], name, obj)
    obj.__module__, obj.__qualname__ = module, name
    return obj

def _unpickle_block(values, placement, ndim): raise RuntimeError("only referenced when writing")
put("pandas._libs.internals", "_unpickle_block", _unpickle_block)
def _new_Index(cls, d): raise RuntimeError("only referenced when writing")
put("pandas.core.indexes.base", "_new_Index", _new_Index)
Index = put("pandas.core.indexes.base", "Index", type("Index", (), {}))

class Block:
    def __init__(self, values, placement): self.values, self.placement = values, placement
    def __reduce__(self): return (_unpickle_block, (self.values, self.placement, 2))
class Names:
    def __init__(self, names): self.names = names
    def __reduce__(self): return (_new_Index, (Index, {"data": np.array(self.names, dtype=object), "name": None}))
class BlockManager:
    def __init__(self, blocks, axes): self.blocks, self.axes = blocks, axes
    def __reduce__(self): return (BlockManager, (tuple(self.blocks), self.axes))
put("pandas.core.internals.managers", "BlockManager", BlockManager)
class DataFrame: pass
put("pandas.core.frame", "DataFrame", DataFrame)

def frame(key, coords, control, own):
    # columns: Chr Pos Control <key> B Start_forward; "B" is an object column with an NA, like the real tables'
    names = ["Chr", "Pos", "Control", key, "B", "Start_forward"]
    blocks = [
        Block(np.array([[c for c, _ in coords], ["0.3", "NA"]], dtype=object), np.array([0, 4])),
        Block(np.array([[p + 1 for _, p in coords], [p for _, p in coords]], dtype=np.int32), np.array([1, 5])),
        Block(np.array([control, own], dtype=float), slice(2, 4, 1)),
    ]
    df = DataFrame()
    df.__dict__["_mgr"] = BlockManager(blocks, [Names(names), Names([0, 1])])
    return df

frames = {
    "A": frame("A", [("1", 99), ("1", 199)], [0.1, 0.2], [0.8, 0.9]),
    "MRXSCJ": frame("MRXSCJ", [("1", 199), ("2", 49)], [0.4, 0.5], [0.6, 0.7]),
}
pickle.dump(frames, open(sys.argv[1], "wb"))
'''

RUN_CONVERT_ARRAY = '''
import json, sys
sys.path.insert(0, sys.argv[1])
import convert_array
convert_array.EXPECTED_CONDITIONS = 2
print(json.dumps(convert_array.convert(sys.argv[2], sys.argv[3])))
'''


@needs_numpy
def test_array_converter_rebuilds_the_tables_and_checks_their_order_against_the_bed(tmp_path):
    py = _python_with("numpy")
    pkl, bed = tmp_path / "refs.pickle", tmp_path / "coords.bed"
    subprocess.run([py, "-c", WRITE_FAKE_FRAMES, str(pkl)], check=True, cwd=tmp_path)
    bed.write_text("1\t99\t100\tA\n1\t199\t200\tA,MRXCJS\n2\t49\t50\tMRXCJS\n")
    done = subprocess.run([py, "-c", RUN_CONVERT_ARRAY, str(CONVERT_ARRAY.parent), str(pkl), str(bed)], capture_output=True, text=True, cwd=tmp_path)
    assert done.returncode == 0, done.stderr
    out = json.loads(done.stdout)["conditions"]
    assert sorted(out) == ["A", "MRXCJS"]                         # the key is renamed as EpiSignaLR renames it
    a = out["A"]
    assert a["coords"] == [["1", 99], ["1", 199]] and a["profiles"] == ["Control", "A", "B"] and a["positive_profile"] == "A"
    assert a["values"] == [[0.1, 0.8, 0.3], [0.2, 0.9, None]]    # "NA" becomes null
    assert out["MRXCJS"]["positive_profile"] == "MRXSCJ"          # the reference column keeps its original name

    bed.write_text("1\t199\t200\tA,MRXCJS\n1\t99\t100\tA\n2\t49\t50\tMRXCJS\n")           # the same CpGs in another order
    done = subprocess.run([py, "-c", RUN_CONVERT_ARRAY, str(CONVERT_ARRAY.parent), str(pkl), str(bed)], capture_output=True, text=True, cwd=tmp_path)
    assert done.returncode != 0 and "bed's order" in done.stderr


@needs_numpy
def test_array_converter_refuses_a_pickle_that_asks_for_anything_else(tmp_path):
    py = _python_with("numpy")
    (tmp_path / "evil.pickle").write_bytes(pickle.dumps(_Evil()))
    (tmp_path / "c.bed").write_text("")
    done = subprocess.run([py, "-I", str(CONVERT_ARRAY), str(tmp_path / "evil.pickle"), str(tmp_path / "c.bed"), str(tmp_path / "o.json")],
                          capture_output=True, text=True, cwd=tmp_path)
    assert done.returncode != 0 and "refusing to load global" in done.stderr
    assert not (tmp_path / "pwned.txt").exists() and not (tmp_path / "o.json").exists()


# ---- the array-derived classifiers: scoring -------------------------------------------------------------------------

SYNTHETIC_REFERENCES = {
    "format": 1,
    "conditions": {
        "X": {
            "positive_profile": "X", "profiles": ["Control", "X", "Y"], "coords": [["1", 10], ["1", 20], ["1", 30]],
            "values": [[0.2, 0.8, 0.5], [0.3, 0.9, 0.4], [0.4, 0.7, 0.6]],
        },
    },
}


def _score_with_svm(genome, tmp_path):
    refs = tmp_path / "refs.json"
    refs.write_text(json.dumps(SYNTHETIC_REFERENCES))
    done = subprocess.run([_python_with("numpy", "sklearn"), "-I", str(ARRAY_SVM), str(refs)], input=json.dumps(genome),
                          capture_output=True, text=True, cwd=tmp_path)
    assert done.returncode == 0, done.stderr
    return json.loads(done.stdout)


@needs_scikit_learn
def test_array_scorer_puts_a_genome_like_the_condition_above_zero_and_like_control_below(tmp_path):
    like_x = _score_with_svm({"1:10": 0.8, "1:20": 0.9, "1:30": 0.7}, tmp_path)
    like_control = _score_with_svm({"1:10": 0.2, "1:20": 0.3, "1:30": 0.4}, tmp_path)
    assert like_x["scores"]["X"] > 0 > like_control["scores"]["X"] and like_x["sklearn"]


@needs_scikit_learn
def test_array_scorer_fills_a_missing_value_with_the_cpgs_mean_over_all_rows(tmp_path):
    # the reference has no gaps, so a genome value that is missing is the mean of that CpG's three reference rows
    rows = SYNTHETIC_REFERENCES["conditions"]["X"]["values"]
    means = {f"1:{p}": sum(rows[i]) / 3 for i, p in enumerate((10, 20, 30))}
    assert _score_with_svm({}, tmp_path)["scores"]["X"] == pytest.approx(_score_with_svm(means, tmp_path)["scores"]["X"])


needs_array = pytest.mark.skipif(not fetched.array_installed(), reason="array-derived classifiers not set up (backend/scripts/setup_episignalr.py)")

# Scores for the genome below from EpiSignaLR's own run_SVM.py (scikit-learn 1.8.0).
SYNTHETIC_GENOME_SCORES = {
    "Kabuki": 0.722136928, "HVDAS_C": -1.082398094, "EEOC": -1.095668339, "MRD51": -1.141009014,
    "ICF2_3_4": -1.144200229, "CdLS": -1.14448125,
}


def _kabuki_genome():
    """The Kabuki median profile at Kabuki's CpGs with every 7th one missing, and nothing elsewhere."""
    kab = arrays.references()["conditions"]["Kabuki"]
    j = kab["profiles"].index("Kabuki")
    genome = {f"{c}:{s}": None for c, s in arrays.sites()}
    for i, ((c, s), row) in enumerate(zip(kab["coords"], kab["values"])):
        genome[f"{c}:{s}"] = None if i % 7 == 0 or row[j] is None else row[j]
    return genome


@needs_array
def test_array_scores_match_episignalrs_own_script():
    scores = arrays.run_scorer(_kabuki_genome())["scores"]
    assert len(scores) == 34
    for name, expected in SYNTHETIC_GENOME_SCORES.items():
        assert scores[name] == pytest.approx(expected, abs=1e-8)


@needs_array
def test_array_reference_data_covers_the_34_conditions_and_3642_cpgs():
    assert set(arrays.references()["conditions"]) == set(arrays.conditions()) and len(arrays.conditions()) == 34
    assert len(arrays.sites()) == 3642


@needs_array
def test_array_overview_shapes_scores_and_calls(monkeypatch, pileup_sample):
    present = {tuple(k.split(":")): v for k, v in _kabuki_genome().items() if v is not None}
    monkeypatch.setattr(cpg, "fetch_sites", lambda bed, sites: {(c, int(s)): CpG(int(s), v * 100, 30, 0, 0) for (c, s), v in present.items()})
    monkeypatch.setattr(cpg, "header", lambda bed: {"pileup_mode": "model", "modsites_mode": "denovo"})
    arrays._cache.clear()
    out = arrays.overview(pileup_sample)
    assert out["call"] == "one" and out["matches"] == ["Kabuki"] and out["threshold"] == 0.35
    top = out["conditions"][0]
    assert top["id"] == "Kabuki" and top["matches"] and top["score"] == pytest.approx(0.722136928, abs=1e-8)
    assert top["name"] == "Kabuki syndrome 1 and 2" and top["genes"] == "KMT2D, KDM6A"
    assert (top["n_cpgs"], top["n_found"]) == (153, 131) and out["cpgs"] == {"total": 3642, "found": 131}
    assert [c["score"] for c in out["conditions"]] == sorted((c["score"] for c in out["conditions"]), reverse=True)
    detail = top["cpgs"]
    assert len(detail["coords"]) == len(detail["control"]) == len(detail["condition"]) == len(detail["sample"]) == 153
    assert detail["coords"][0] == ["1", 935928] and detail["sample"][0] is None     # 1-based position; every 7th is missing
    assert [w["id"] for w in out["warnings"]] == ["pileup_mode"] and out["scorer"]["reference_scikit_learn"] == "1.8.0"


@needs_array
def test_array_call_is_control_when_nothing_reaches_the_threshold(monkeypatch, pileup_sample):
    monkeypatch.setattr(cpg, "fetch_sites", lambda bed, sites: {})
    monkeypatch.setattr(cpg, "header", lambda bed: {"pileup_mode": "count", "modsites_mode": "reference"})
    arrays._cache.clear()
    out = arrays.overview(pileup_sample)
    assert out["call"] == "control" and out["matches"] == [] and out["warnings"] == [] and out["cpgs"]["found"] == 0


# ---- reading signature CpGs from a pileup ----------------------------------------------------------------------------

def test_fetch_sites_reads_all_sites_in_one_tabix_call_and_matches_the_files_chromosome_names(monkeypatch):
    seen = {}

    def fake_run(args):
        seen["args"] = args
        seen["regions"] = open(args[1]).read()
        return "chr2\t49\t50\t55.5\tTotal\t30\t17\t13\t55.0\nchr1\t99\t100\t5.0\tTotal\t20\t1\t19\t5.0\n"

    monkeypatch.setattr(cpg, "_run", fake_run)
    monkeypatch.setattr(cpg, "contigs", lambda bed: ["chr1", "chr2"])
    got = cpg.fetch_sites("p.bed.gz", [("2", 49), ("1", 99), ("1", 500)])
    assert seen["args"][0] == "-R" and seen["args"][2] == "p.bed.gz"
    assert seen["regions"] == "chr1\t99\t100\nchr1\t500\t501\nchr2\t49\t50\n"           # sorted, "chr" added, 0-based BED
    assert set(got) == {("1", 99), ("2", 49)} and got[("2", 49)].mod_score == 55.5 and got[("1", 99)].coverage == 20


def test_fetch_sites_leaves_names_alone_when_the_file_has_no_chr_prefix(monkeypatch):
    monkeypatch.setattr(cpg, "_run", lambda args: "1\t99\t100\t5.0\tTotal\t20\t1\t19\t5.0\n")
    monkeypatch.setattr(cpg, "contigs", lambda bed: ["1", "2"])
    assert set(cpg.fetch_sites("p.bed.gz", [("1", 99)])) == {("1", 99)}


# ---- the fetched models -----------------------------------------------------------------------------------------

# Probability the authors' scikit-learn pipeline gives the first and last row of each model's own training data.
SKLEARN_TRAINING_ROWS = {
    "kabuki": (0.7556585488, 0.0725952239),
    "kleefstra": (0.9228247966, 0.116368815),
    "setd1b": (0.9121430531, 0.0787635412),
}


@needs_vendor
@pytest.mark.parametrize("key", list(fetched.CATALOG))
def test_fetched_model_scores_its_training_rows_as_scikit_learn_does(key):
    m = model.load(key)
    rows = model.training(key)["rows"]
    for row, expected in zip((rows[0], rows[-1]), SKLEARN_TRAINING_ROWS[key]):
        assert m.score({d: row[d] for d in m.features}).probability == pytest.approx(expected, abs=1e-8)


@needs_vendor
@pytest.mark.parametrize("key", list(fetched.CATALOG))
def test_fetched_model_is_consistent_with_its_panel_and_metadata(key):
    m = model.load(key)
    panel = regions.parse_panel(fetched.model_files(fetched.CATALOG[key][0])["dmrs"])
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
    assert [m["id"] for m in out["models"]] == list(fetched.CATALOG)
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
    monkeypatch.setattr(router.episignature_fetched, "installed", lambda: False)
    body = client.get("/api/epigenetics/episignatures").json()
    assert body["installed"] is False and body["models"] == []
    assert body["setup_command"] == "python backend/scripts/setup_episignalr.py"
    assert {a["id"] for a in body["available"]} == {"setd1b", "kleefstra", "kabuki"}


def test_route_returns_the_service_result(client, monkeypatch):
    monkeypatch.setattr(router.episignature_fetched, "installed", lambda: True)
    monkeypatch.setattr(router.episignature_service, "overview", lambda s: {"sample": s.key, "models": [{"id": "kabuki"}]})
    body = client.get("/api/epigenetics/episignatures", params={"sample": "SYN"}).json()
    assert body == {"installed": True, "sample": "SYN", "models": [{"id": "kabuki"}]}
    assert client.get("/api/epigenetics/episignatures", params={"sample": "nope"}).status_code == 404


def test_route_reports_a_missing_tabix_and_inconsistent_files_as_unavailable(client, monkeypatch):
    monkeypatch.setattr(router.episignature_fetched, "installed", lambda: True)

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


def test_array_route_says_how_to_install_when_the_files_are_not_there(client, monkeypatch):
    monkeypatch.setattr(router.episignature_fetched, "array_installed", lambda: False)
    body = client.get("/api/epigenetics/episignatures/array").json()
    assert body["installed"] is False and body["conditions"] == [] and body["setup_command"].endswith("setup_episignalr.py")


def test_array_route_returns_the_scores_and_reports_failures_as_unavailable(client, monkeypatch):
    monkeypatch.setattr(router.episignature_fetched, "array_installed", lambda: True)
    monkeypatch.setattr(router.episignature_arrays, "overview", lambda s: {"sample": s.key, "call": "control", "conditions": []})
    assert client.get("/api/epigenetics/episignatures/array").json() == {"installed": True, "sample": "SYN", "call": "control", "conditions": []}

    def broken(error):
        def overview(s):
            raise error
        return overview

    monkeypatch.setattr(router.episignature_arrays, "overview", broken(arrays.ScoringError("scorer failed")))
    r = client.get("/api/epigenetics/episignatures/array")
    assert r.status_code == 503 and r.json()["detail"] == "scorer failed"
    monkeypatch.setattr(router.episignature_arrays, "overview", broken(cpg.TabixError("tabix not found")))
    assert client.get("/api/epigenetics/episignatures/array").status_code == 503
