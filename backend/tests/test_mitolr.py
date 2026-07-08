"""Tests for MitoLR: seeded demo generation, real MitoRSaw VCF parsing, the
sample registry, and the API routes."""

from fastapi.testclient import TestClient

from mitolr import annotate, generate, qc, samples, vcf
from repeatscope.api.main import app

client = TestClient(app)

HG002_VCF = samples.SAMPLES["HG002"]["vcf"]


# --- seeded demo sample -----------------------------------------------------

def test_demo_variants_deterministic_and_counted():
    a = generate.build_variants()
    b = generate.build_variants()
    assert [v["id"] for v in a] == [v["id"] for v in b]
    assert len(a) == 79
    assert sum(1 for v in a if v["origin"] == "pathogenic") == 21
    assert sum(1 for v in a if v["origin"] == "novel") == 58


def test_demo_narrative_overrides():
    by_id = {v["id"]: v for v in generate.build_variants()}
    assert by_id["m3243AG"]["heteroplasmy"] == 38.0
    assert by_id["m3243AG"]["status"] == "Pathogenic"
    assert by_id["kss_common_del"]["heteroplasmy"] == 11.0


# --- real MitoRSaw VCF parsing ---------------------------------------------

def test_hg002_vcf_parses():
    vs = vcf.build_variants(HG002_VCF)
    assert len(vs) == 15
    assert [v["pos"] for v in vs] == sorted(v["pos"] for v in vs)
    # HG002 has no curated-pathogenic calls; all are novel/benign.
    assert {v["origin"] for v in vs} == {"novel"}


def test_hg002_heteroplasmy_from_vaf_and_depth_from_dp():
    by_pos = {v["pos"]: v for v in vcf.build_variants(HG002_VCF)}
    # m.263A>G is homoplasmic (VAF 1.0) at DP 1108.
    assert by_pos[263]["heteroplasmy"] == 100.0
    assert by_pos[263]["depth"] == 1108
    assert by_pos[263]["change"] == "A>G"
    # m.302 A>AC is a control-region insertion at VAF 0.9369.
    assert by_pos[302]["heteroplasmy"] == 93.7
    assert by_pos[302]["change"].startswith("ins")


def test_hg002_locus_annotation():
    by_pos = {v["pos"]: v for v in vcf.build_variants(HG002_VCF)}
    assert by_pos[750]["locus"] == "MT-RNR1"
    assert by_pos[4769]["locus"] == "MT-ND2"
    assert by_pos[15326]["locus"] == "MT-CYB"
    assert by_pos[263]["locus"] == "D-loop"
    assert by_pos[750]["consequence"] == "rRNA variant"


def test_pileup_ratio_matches_allele_depths():
    by_pos = {v["pos"]: v for v in vcf.build_variants(HG002_VCF)}
    v = by_pos[302]  # AD = 70 ref, 1040 alt
    reads = vcf.simulate_pileup(v)
    assert len(reads) == 48
    frac = sum(1 for r in reads if r["kind"] == "mut") / len(reads)
    assert abs(frac - 1040 / (70 + 1040)) < 0.03


def test_variant_notation_forms():
    assert annotate.variant_notation(263, "A", "G") == ("m.263A>G", "A>G", False)
    pos, change, is_del = annotate.variant_notation(302, "A", "AC")
    assert change == "insC" and not is_del
    pos, change, is_del = annotate.variant_notation(3106, "CN", "C")
    assert change == "del" and not is_del  # 1 bp deletion, not the large-arc kind


# --- sample registry / dispatch --------------------------------------------

def test_default_sample_is_hg002():
    assert samples.resolve(None) == "HG002"


def test_overview_dispatch_per_sample():
    hg = samples.overview("HG002")
    assert hg["stats"]["variants_called"] == 15
    assert hg["stats"]["pathogenic_findings"] == 0
    assert hg["haplogroup"]["name"] == "Unclassified"
    assert hg["caller"].startswith("mitorsaw")

    demo = samples.overview("MT-0417-P2")
    assert demo["stats"]["variants_called"] == 79
    assert demo["haplogroup"]["name"] == "H2a5b1"


def test_qc_uses_bam_when_available():
    # HG002 has a precomputed BAM-stats fixture -> rich, real QC.
    hg = samples.sample_qc("HG002")
    assert len(hg["read_length_histogram"]) == 18  # real, not empty
    assert len(hg["coverage_track"]) > 100  # continuous per-base track
    labels = [s["label"] for s in hg["stats"]]
    assert "Mean coverage" in labels and "chrM reads mapped" in labels
    assert hg["coverage_note"].startswith("Per-base read depth")

    demo = samples.sample_qc("MT-0417-P2")
    assert len(demo["read_length_histogram"]) == 18
    assert len(demo["coverage_track"]) == 60


def test_bam_qc_suppresses_trna_strand_balance():
    hg = samples.sample_qc("HG002")
    by_region = {r["region"]: r for r in hg["regions"]}
    # tiny tRNA locus -> strand balance shown as "—" (midpoint approximation).
    assert by_region["MT-TF"]["strandBalance"] == "—"
    # large rRNA/protein locus keeps its computed value.
    assert by_region["MT-RNR2"]["strandBalance"] != "—"


def test_vcf_qc_fallback_without_bam(monkeypatch):
    # If the BAM fixture is absent, HG002 falls back to VCF-derived QC.
    monkeypatch.setitem(samples.SAMPLES["HG002"], "bam_qc", "/nonexistent.json")
    samples._bam_qc.cache_clear()
    hg = samples.sample_qc("HG002")
    assert hg["read_length_histogram"] == []
    samples._bam_qc.cache_clear()


def test_snv_pileup_uses_bam_indel_uses_caller():
    snv = samples.variant_detail("HG002", "chrM_263_A_G")
    assert "aligned HiFi BAM" in snv["pileup_source"]
    assert all(r["kind"] == "mut" for r in snv["pileup_reads"])  # homoplasmic G

    indel = samples.variant_detail("HG002", "chrM_302_A_AC")
    assert "haplotype-aware" in indel["pileup_source"]


# --- API routes -------------------------------------------------------------

def test_samples_endpoint():
    body = client.get("/api/mito/samples").json()
    assert body["default"] == "HG002"
    keys = [s["key"] for s in body["samples"]]
    assert "HG002" in keys and "MT-0417-P2" in keys


def test_overview_endpoint_defaults_to_hg002():
    body = client.get("/api/mito/overview").json()
    assert body["stats"]["variants_called"] == 15


def test_overview_endpoint_demo_sample():
    body = client.get("/api/mito/overview?sample=MT-0417-P2").json()
    assert body["stats"]["variants_called"] == 79


def test_variant_endpoint_and_404():
    vs = client.get("/api/mito/overview").json()["variants"]
    vid = vs[0]["id"]
    r = client.get(f"/api/mito/variant/{vid}")
    assert r.status_code == 200
    assert r.json()["variant"]["id"] == vid
    assert client.get("/api/mito/variant/nope").status_code == 404
    assert client.get("/api/mito/overview?sample=bogus").status_code == 404


def test_repeatscope_routes_still_present():
    assert client.get("/api/modes").status_code == 200
