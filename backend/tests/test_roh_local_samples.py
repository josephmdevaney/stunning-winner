"""The ROH suite's local sample registry: a gitignored config naming this machine's own ROH files."""

import json

from roh import analysis


def test_only_samples_whose_file_exists_are_registered(tmp_path):
    present = tmp_path / "a.roh.txt"
    present.write_text("")
    config = tmp_path / "roh.local.json"
    config.write_text(json.dumps({"samples": {"A": str(present), "B": str(tmp_path / "missing.roh.txt")}}))
    assert analysis._local_samples(config) == {"A": str(present)}


def test_no_config_registers_nothing(tmp_path):
    assert analysis._local_samples(tmp_path / "absent.json") == {}


def test_a_corrupt_config_registers_nothing(tmp_path):
    config = tmp_path / "roh.local.json"
    config.write_text("{not json")
    assert analysis._local_samples(config) == {}


def test_a_config_without_samples_registers_nothing(tmp_path):
    config = tmp_path / "roh.local.json"
    config.write_text("{}")
    assert analysis._local_samples(config) == {}
