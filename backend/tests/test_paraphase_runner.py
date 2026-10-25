"""Paraphase Annotation runner: vendored-script versioning, input preparation and
run provenance (paraphase_annotate/runner.py). Uses a stand-in `interpret_paraphase.py`
and synthetic JSON in a temp directory -- no real sample data and no network."""

import json
import sys
import time

import pytest

from paraphase_annotate import runner, samples

SAMPLE = "S1"

V3_JSON = {
    "smn1": {
        "smn1_cn": 2, "smn2_cn": 0, "total_cn": None, "highest_total_cn": 2,
        "final_haplotypes": {"h1": "11"}, "region_depth": {"median": 60.0},
    },
    "AMY1A": {"total_cn": 4, "final_haplotypes": {}, "region_depth": {"median": 120.0}},
}

# Writes what the real script writes (two CSVs + its completion line) and records
# its own command line, so the test can see exactly how it was launched.
STAND_IN_SCRIPT = '''
import argparse, json, os, sys
ap = argparse.ArgumentParser()
ap.add_argument("--input"); ap.add_argument("--out"); ap.add_argument("--region"); ap.add_argument("--config")
a = ap.parse_args()
sid = os.path.basename(a.input).split(".paraphase.json")[0]
vcfs = os.path.join(os.path.dirname(a.input), sid + "_paraphase_vcfs")
json.dump({"argv": sys.argv[1:], "vcfs_exists": os.path.isdir(vcfs),
           "json_keys": sorted(json.load(open(a.input))["smn1"])}, open(os.path.join(a.out, "seen.json"), "w"))
open(os.path.join(a.out, sid + "_summary.csv"), "w").write("Region\\nsmn1\\n")
open(os.path.join(a.out, sid + "_detailed.csv"), "w").write("Region\\nsmn1\\n")
print("Completed annotating all regions")
'''


@pytest.fixture
def env(tmp_path, monkeypatch):
    # sample: a v3-shaped JSON with a sibling VCF directory (as Paraphase lays it out)
    data = tmp_path / "data"
    (data / f"{SAMPLE}_paraphase_vcfs").mkdir(parents=True)
    (data / f"{SAMPLE}_paraphase_vcfs" / f"{SAMPLE}_smn1.vcf").write_text(
        "##fileformat=VCFv4.2\n##paraphase_version=3.5.0\n##paraphase_command=paraphase --bam x.bam\n#CHROM\n"
    )
    (data / f"{SAMPLE}.paraphase.json").write_text(json.dumps(V3_JSON))
    monkeypatch.setattr(samples, "SAMPLES", {SAMPLE: str(data / f"{SAMPLE}.paraphase.json")})

    # vendored script: stand-in, tagged v4.1.0
    vendor = tmp_path / "vendor"
    (vendor / "data").mkdir(parents=True)
    script = vendor / "interpret_paraphase.py"
    script.write_text(STAND_IN_SCRIPT)
    (vendor / "VERSION.json").write_text(json.dumps({
        "paraphase_tag": "v4.1.0", "fetched_at": "2026-10-25T00:00:00+00:00", "source": "x", "sha256": {},
    }))
    monkeypatch.setattr(runner, "VENDOR_DIR", vendor)
    monkeypatch.setattr(runner, "VENDOR_SCRIPT", script)
    monkeypatch.setattr(runner, "VENDOR_CONFIG", vendor / "data" / "config_grch38.yaml")
    monkeypatch.setattr(runner, "VENDOR_VERSION_FILE", vendor / "VERSION.json")
    monkeypatch.setattr(runner, "RESULTS_ROOT", tmp_path / "results")
    monkeypatch.setattr(runner, "interpreter_for_script", lambda: sys.executable)
    runner._RUNNING.clear()
    return tmp_path


def _wait_done(sample, timeout=20):
    deadline = time.time() + timeout
    while time.time() < deadline:
        st = runner.status(sample)
        if st["state"] != "running":
            return st
        time.sleep(0.1)
    raise AssertionError("run did not finish")


def test_run_adapts_v3_json_for_the_v4_script_and_records_provenance(env):
    assert runner.start_run(SAMPLE) == {"state": "running"}
    assert _wait_done(SAMPLE) == {"state": "done"}

    seen = json.loads((runner.results_dir(SAMPLE) / "seen.json").read_text())
    assert seen["vcfs_exists"] is True  # the script found the VCFs through the symlink
    # the script was handed the normalized JSON: region-specific fields moved under region_specific_info
    assert "region_specific_info" in seen["json_keys"] and "smn1_cn" not in seen["json_keys"]
    assert str(runner.input_dir(SAMPLE)) in seen["argv"][seen["argv"].index("--input") + 1]

    meta = runner.run_meta(SAMPLE)
    assert meta["script_tag"] == "v4.1.0"
    assert meta["input_schema"] == "v3" and meta["script_schema"] == "v4" and meta["adapted_json"] is True
    assert meta["data_paraphase_version"] == "3.5.0"
    assert "finished_at" in meta  # set once the run is done
    # the default medical regions that the JSON doesn't have are reported, not hidden
    assert "rccx" in meta["regions_missing_from_json"] and "smn1" not in meta["regions_missing_from_json"]


def test_the_original_sample_directory_is_never_modified(env):
    before = sorted(p.name for p in (env / "data").rglob("*"))
    runner.start_run(SAMPLE)
    _wait_done(SAMPLE)
    assert sorted(p.name for p in (env / "data").rglob("*")) == before
    assert json.loads((env / "data" / f"{SAMPLE}.paraphase.json").read_text()) == V3_JSON


def test_v3_vendored_script_gets_the_json_unchanged(env, monkeypatch):
    (runner.VENDOR_VERSION_FILE).write_text(json.dumps({"paraphase_tag": "v3.5.0"}))
    runner.start_run(SAMPLE)
    _wait_done(SAMPLE)
    seen = json.loads((runner.results_dir(SAMPLE) / "seen.json").read_text())
    assert "smn1_cn" in seen["json_keys"]
    assert runner.run_meta(SAMPLE)["adapted_json"] is False


def test_v4_json_is_blocked_for_a_v3_vendored_script(env):
    v4 = {"smn1": {"total_cn": 2, "region_specific_info": {"smn1_cn": 2}}}
    (env / "data" / f"{SAMPLE}.paraphase.json").write_text(json.dumps(v4))
    runner.VENDOR_VERSION_FILE.write_text(json.dumps({"paraphase_tag": "v3.5.0"}))
    reason = runner.preflight(SAMPLE)
    assert "Paraphase v4+" in reason and "setup_paraphase_annotate.py" in reason
    assert runner.start_run(SAMPLE)["state"] == "blocked"


def test_blocked_when_no_suitable_python_is_found(env, monkeypatch):
    monkeypatch.setattr(runner, "interpreter_for_script", lambda: None)
    reason = runner.preflight(SAMPLE)
    assert "needs Python 3.12+" in reason and "PARAPHASE_ANNOTATE_PYTHON" in reason


def test_a_script_vendored_before_version_tracking_is_treated_as_v3_5_0(env):
    runner.VENDOR_VERSION_FILE.unlink()
    assert runner.vendor_info()["paraphase_tag"] == "v3.5.0"
    assert runner.script_major() == 3


def test_vendor_info_is_none_when_nothing_is_vendored(env):
    runner.VENDOR_SCRIPT.unlink()
    assert runner.vendor_info() is None


def test_run_meta_is_none_before_any_run(env):
    assert runner.run_meta(SAMPLE) is None


def test_a_rerun_archives_the_previous_results_instead_of_deleting_them(env):
    runner.start_run(SAMPLE)
    _wait_done(SAMPLE)
    old_summary = runner.summary_csv(SAMPLE).read_text()
    runner.summary_csv(SAMPLE).write_text("Region\nfrom the first run\n")

    runner.start_run(SAMPLE)
    _wait_done(SAMPLE)

    archives = sorted((runner.results_dir(SAMPLE) / "archive").iterdir())
    assert len(archives) == 1
    kept = {p.name for p in archives[0].iterdir()}
    assert {f"{SAMPLE}_summary.csv", f"{SAMPLE}_detailed.csv", "run.log", "run_meta.json"} <= kept
    assert (archives[0] / f"{SAMPLE}_summary.csv").read_text() == "Region\nfrom the first run\n"
    assert runner.summary_csv(SAMPLE).read_text() == old_summary  # the new run wrote fresh results


def test_the_first_run_has_nothing_to_archive(env):
    assert runner.archive_previous(SAMPLE) is None
    runner.start_run(SAMPLE)
    _wait_done(SAMPLE)
    assert not (runner.results_dir(SAMPLE) / "archive").exists()
