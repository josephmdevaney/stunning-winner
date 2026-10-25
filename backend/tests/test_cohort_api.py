import pytest
from fastapi.testclient import TestClient

from repeatscope.api.main import app


@pytest.fixture(scope="module")
def client():
    with TestClient(app) as c:
        yield c


@pytest.mark.parametrize("mode", ["carrier", "puretarget_v2"])
def test_cohort_aggregates_every_sample_in_the_mode(client, mode):
    cohort = client.get("/api/cohort", params={"mode": mode}).json()
    samples = client.get("/api/samples", params={"mode": mode}).json()["samples"]
    assert cohort["stats"]["samples"] == len(cohort["samples"]) == len(samples)
    assert cohort["stats"]["qc_pass"] + cohort["stats"]["qc_fail"] == len(samples)
    assert sum(b["count"] for b in cohort["depth_histogram"]) <= len(samples)
    assert cohort["stats"]["reportable_findings"] == sum(len(r["reportable"]) for r in cohort["samples"])


def test_cohort_reportable_matches_overview_pathogenic_count(client):
    cohort = client.get("/api/cohort", params={"mode": "puretarget_v2"}).json()
    for row in cohort["samples"][:5]:
        overview = client.get("/api/overview", params={"mode": "puretarget_v2", "sample": row["sample"]}).json()
        assert len(row["reportable"]) == overview["stats"]["pathogenic_findings"]


def test_cohort_unavailable_for_genome_wide(client):
    assert client.get("/api/cohort", params={"mode": "genome_wide"}).status_code == 404
