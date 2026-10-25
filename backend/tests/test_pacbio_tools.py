"""PacBio tool-versions report (pacbio_tools/). Synthetic headers and stubbed
GitHub responses only -- no network and no real sample data."""

import gzip
import json
import urllib.error

import pytest

from pacbio_tools import observed, registry, releases, report, versions

# ---- version parsing and comparison -----------------------------------------

@pytest.mark.parametrize("text, expected", [
    ("5.1.0", (5, 1, 0)),
    ("v3.5.0", (3, 5, 0)),
    ("sawfish 2.0.3", (2, 0, 3)),
    ("0.2.7-8ee9b5b", (0, 2, 7)),
    ("26.2.0", (26, 2, 0)),
    ("2.11", (2, 11, 0)),
    ("TRGT", None),
    ("", None),
    (None, None),
])
def test_parse(text, expected):
    assert versions.parse(text) == expected


def test_clean_strips_names_prefixes_and_build_ids():
    assert versions.clean("v4.1.0") == "4.1.0"
    assert versions.clean("mitorsaw 0.2.7-8ee9b5b") == "0.2.7"
    assert versions.clean("no version here") is None


@pytest.mark.parametrize("observed_v, latest_v, expected", [
    ("3.5.0", "4.1.0", "major"),
    ("0.2.7", "0.2.14", "patch"),       # 0.2.x vs 0.2.y: same major and minor
    ("2.0.3", "2.3.0", "minor"),
    ("5.1.0", "5.1.2", "patch"),
    ("5.1.0", "5.1.0", "none"),
    ("5.2.0", "5.1.0", "ahead"),
    ("garbage", "5.1.0", None),
    ("5.1.0", None, None),
])
def test_gap(observed_v, latest_v, expected):
    assert versions.gap(observed_v, latest_v) == expected


def test_minor_numbers_compare_as_integers_not_text():
    assert versions.gap("2.9.0", "2.11.0") == "minor"  # "9" > "11" as text; 9 < 11 as numbers


# ---- readers ----------------------------------------------------------------

def _write_vcf(path, header_lines):
    text = "##fileformat=VCFv4.2\n" + "".join(l + "\n" for l in header_lines) + "#CHROM\tPOS\n"
    if str(path).endswith(".gz"):
        with gzip.open(path, "wt") as f:
            f.write(text)
    else:
        path.write_text(text)


def test_sawfish_version_reads_the_source_header(tmp_path):
    p = tmp_path / "a.sv.vcf.gz"
    _write_vcf(p, ['##source="sawfish 2.0.3"', '##sawfish_cmdline="sawfish joint-call"'])
    assert observed.sawfish_version(p) == "sawfish 2.0.3"


def test_sawfish_version_ignores_another_tools_vcf(tmp_path):
    p = tmp_path / "x.vcf"
    _write_vcf(p, ['##source="DeepVariant 1.6"'])
    assert observed.sawfish_version(p) is None
    assert observed.sawfish_version(tmp_path / "missing.vcf") is None


def test_mitorsaw_version_from_a_recorded_caller_string():
    assert observed.mitorsaw_version("mitorsaw 0.2.7-8ee9b5b") == "0.2.7-8ee9b5b"
    assert observed.mitorsaw_version("MitoCaller v1.4 (synthetic)") is None
    assert observed.mitorsaw_version(None) is None


def test_group_counts_samples_per_label_and_version_and_keeps_unreadable_ones_apart():
    out = observed.group([
        ("Carrier", "5.0.0", "s1"), ("Carrier", "5.0.0", "s2"), ("Carrier", None, "s3"),
        ("Genome-wide", "5.1.0", "HG002"),
    ])
    by = {(o["label"], o["raw"]): o for o in out}
    assert by[("Carrier", "5.0.0")]["count"] == 2
    assert by[("Carrier", None)]["count"] == 1  # a partly-recorded cohort doesn't look fully recorded
    assert by[("Genome-wide", "5.1.0")]["examples"] == ["HG002"]


# ---- report -----------------------------------------------------------------

LATEST = {
    "source": "snapshot", "fetched_at": "2026-10-25T00:00:00+00:00", "errors": {},
    "releases": {
        "paraphase": {"tag": "v4.1.0", "version": "4.1.0", "published": "2026-09-25", "name": "x", "url": "u"},
        "trgt": {"tag": "v5.1.0", "version": "5.1.0", "published": "2026-06-10", "name": "x", "url": "u"},
        "sawfish": {"tag": "v2.3.0", "version": "2.3.0", "published": "2026-10-21", "name": "x", "url": "u"},
        "hiphase": {"tag": "v1.7.0", "version": "1.7.0", "published": "2026-05-27", "name": "x", "url": "u"},
    },
}


def _tool(rep, key):
    return next(t for t in rep["tools"] if t["key"] == key)


def test_report_statuses():
    rep = report.build({
        "paraphase": [{"label": "A", "raw": "3.5.0", "count": 2, "examples": []},
                      {"label": "script", "raw": "v4.1.0", "count": 1, "examples": []}],
        "trgt": [{"label": "G", "raw": "5.1.0", "count": 1, "examples": []}],
        "sawfish": [{"label": "C", "raw": None, "count": 3, "examples": []}],
    }, LATEST)
    assert _tool(rep, "paraphase")["status"] == "behind" and _tool(rep, "paraphase")["worst_gap"] == "major"
    assert _tool(rep, "trgt")["status"] == "current"
    assert _tool(rep, "sawfish")["status"] == "unknown"          # files exist, no version in them
    assert _tool(rep, "pbmm2")["status"] == "not_recorded"       # used, but nothing in the data records it
    assert _tool(rep, "hiphase")["status"] == "not_used"
    assert rep["summary"]["used"] == sum(t.used for t in registry.TOOLS)
    assert rep["summary"]["behind"] == 1 and rep["summary"]["current"] == 1


def test_report_attaches_the_cleaned_observed_version_and_the_gap():
    rep = report.build({"paraphase": [{"label": "A", "raw": "3.5.0", "count": 1, "examples": []}]}, LATEST)
    o = _tool(rep, "paraphase")["observed"][0]
    assert (o["version"], o["gap"]) == ("3.5.0", "major")


def test_a_tool_with_no_known_latest_release_has_no_gap():
    rep = report.build({"pbmm2": [{"label": "x", "raw": "26.1.0", "count": 1, "examples": []}]}, LATEST)
    assert _tool(rep, "pbmm2")["latest"] is None
    assert _tool(rep, "pbmm2")["observed"][0]["gap"] is None


# ---- releases ---------------------------------------------------------------

@pytest.fixture
def isolated(tmp_path, monkeypatch):
    snap = tmp_path / "snapshot.json"
    snap.write_text(json.dumps({"fetched_at": "2026-01-01T00:00:00+00:00", "releases": {
        "trgt": {"tag": "v5.0.0", "version": "5.0.0", "published": "2026-01-01", "name": "t", "url": "u"},
    }}))
    monkeypatch.setattr(releases, "SNAPSHOT_PATH", snap)
    monkeypatch.setattr(releases, "CACHE_PATH", tmp_path / "cache.json")
    return tmp_path


def _fake_live(monkeypatch, live, errors=None):
    monkeypatch.setattr(releases, "fetch_live", lambda: (live, errors or {}))


def test_offline_falls_back_to_the_snapshot_and_says_so(isolated, monkeypatch):
    _fake_live(monkeypatch, {}, {"trgt": "URLError: offline"})
    out = releases.latest_releases()
    assert out["source"] == "snapshot" and out["fetched_at"] == "2026-01-01T00:00:00+00:00"
    assert out["releases"]["trgt"]["version"] == "5.0.0"
    assert "trgt" in out["errors"]


def test_a_live_result_is_cached_and_overlays_the_snapshot(isolated, monkeypatch):
    new = {"tag": "v5.1.0", "version": "5.1.0", "published": "2026-06-10", "name": "t", "url": "u"}
    _fake_live(monkeypatch, {"trgt": new})
    out = releases.latest_releases(refresh=True)
    assert out["source"] == "live" and out["releases"]["trgt"]["version"] == "5.1.0"
    # a second call inside the TTL is served from the cache without another live query
    monkeypatch.setattr(releases, "fetch_live", lambda: pytest.fail("should not query GitHub again"))
    again = releases.latest_releases()
    assert again["source"] == "cache" and again["releases"]["trgt"]["version"] == "5.1.0"


def test_one_failing_tool_keeps_its_previous_release_and_is_reported(isolated, monkeypatch):
    sawfish = {"tag": "v2.3.0", "version": "2.3.0", "published": "2026-10-21", "name": "s", "url": "u"}
    _fake_live(monkeypatch, {"sawfish": sawfish}, {"trgt": "HTTPError: 403"})
    out = releases.latest_releases(refresh=True)
    assert out["releases"]["sawfish"]["version"] == "2.3.0"
    assert out["releases"]["trgt"]["version"] == "5.0.0"  # kept from the snapshot
    assert out["errors"] == {"trgt": "HTTPError: 403"}


def test_fetch_one_falls_back_to_the_newest_tag_when_a_repo_has_no_releases(monkeypatch):
    def fake_get(url):
        if url.endswith("/releases/latest"):
            raise urllib.error.HTTPError(url, 404, "Not Found", {}, None)
        return [{"name": "v1.2.3"}]
    monkeypatch.setattr(releases, "_get_json", fake_get)
    out = releases.fetch_one("PacificBiosciences/example")
    assert out["tag"] == "v1.2.3" and out["version"] == "1.2.3" and out["published"] is None


def test_every_registered_tool_is_in_the_committed_snapshot():
    snap = json.loads(releases.SNAPSHOT_PATH.read_text())
    assert set(snap["releases"]) >= {t.key for t in registry.TOOLS}


# ---- route ------------------------------------------------------------------

def test_versions_route_combines_observed_versions_with_the_newest_releases(monkeypatch):
    from fastapi.testclient import TestClient

    from pacbio_tools import router as tools_router
    from repeatscope.api.main import app

    monkeypatch.setattr(tools_router, "collect", lambda: {
        "trgt": [{"label": "Genome-wide", "raw": "5.0.0", "count": 1, "examples": ["HG002"]}],
    })
    monkeypatch.setattr(releases, "latest_releases", lambda refresh=False: {**LATEST, "refresh_seen": refresh})

    with TestClient(app) as client:
        body = client.get("/api/tools/versions").json()
        trgt = next(t for t in body["tools"] if t["key"] == "trgt")
        assert trgt["status"] == "behind" and trgt["worst_gap"] == "minor"
        assert trgt["observed"][0]["version"] == "5.0.0" and trgt["latest"]["tag"] == "v5.1.0"
        assert body["source"] == "snapshot" and body["summary"]["behind"] == 1
        assert client.get("/api/tools/versions", params={"refresh": "true"}).status_code == 200
