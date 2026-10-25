"""Epigenetics suite: CpG reading, pooled statistics, region calls, disorder patterns and MLID.

Synthetic data only. The one real-file test builds a tiny bgzipped, tabix-indexed bed in a temp directory
and is skipped when htslib's tabix/bgzip aren't installed."""

import json
import shutil
import subprocess

import pytest

from epigenetics import analysis, cpg, disorders, samples, service
from epigenetics.analysis import (
    STATUS_BORDERLINE, STATUS_GAIN, STATUS_IMPRINTED, STATUS_LOSS, STATUS_LOW_COVERAGE, STATUS_NO_DATA,
    Thresholds, TrackStats, classify, region_result, track_stats,
)
from epigenetics.cpg import CpG


def _cpgs(n, pct, cov=30, start=1000):
    """n CpGs, each with `pct` percent methylated at depth `cov`."""
    mod = round(cov * pct / 100)
    return [CpG(start + 10 * i, float(pct), cov, mod, cov - mod) for i in range(n)]


def stats(pct, n=20, cov=30):
    return TrackStats(n_cpgs=n, mean_coverage=float(cov), methylation=float(pct))


# ---- pooled statistics --------------------------------------------------------------------------------

def test_track_stats_pools_counts_so_deep_cpgs_weigh_more():
    shallow_high = CpG(1, 80.0, 10, 8, 2)
    deep_low = CpG(2, 1.0, 100, 1, 99)
    s = track_stats([shallow_high, deep_low])
    assert s.methylation == pytest.approx(100 * 9 / 110)            # pooled, not the 40.5% mean of the two scores
    assert s.n_cpgs == 2 and s.mean_coverage == 55.0


def test_track_stats_is_none_without_read_support():
    assert track_stats([]) is None
    assert track_stats([CpG(1, 0.0, 0, 0, 0)]) is None


# ---- region calls -------------------------------------------------------------------------------------

def test_a_normal_phased_imprint_is_allele_specific():
    r = classify(stats(48), stats(95), stats(3))
    assert (r["status"], r["phasing"], r["pattern"], r["allelic_difference"]) == (STATUS_IMPRINTED, "phased", "allele-specific", 92.0)


def test_phased_loss_and_gain_need_both_alleles_to_agree():
    assert classify(stats(4), stats(5), stats(3))["status"] == STATUS_LOSS
    assert classify(stats(4), stats(5), stats(3))["pattern"] == "both alleles unmethylated"
    assert classify(stats(96), stats(97), stats(95))["status"] == STATUS_GAIN
    assert classify(stats(50), stats(45), stats(60))["status"] == STATUS_BORDERLINE   # similar, intermediate


def test_allele_split_wins_over_the_pooled_value():
    # pooled 25% would read as "not hemimethylated" on its own; a 0 / 50 split is still an allelic imprint
    assert classify(stats(25), stats(0), stats(50))["status"] == STATUS_IMPRINTED


def test_unphased_calls_fall_back_to_bands_and_say_so():
    for pct, want in [(50, STATUS_IMPRINTED), (10, STATUS_LOSS), (90, STATUS_GAIN), (25, STATUS_BORDERLINE), (75, STATUS_BORDERLINE)]:
        r = classify(stats(pct), None, None)
        assert (r["status"], r["phasing"], r["pattern"]) == (want, "unphased", "unphased"), pct


def test_one_unusable_haplotype_means_unphased_with_a_note():
    r = classify(stats(50), stats(95), stats(3, n=2))   # hap2 has too few CpGs
    assert r["phasing"] == "unphased" and r["status"] == STATUS_IMPRINTED
    assert any("not usable" in n for n in r["notes"])


def test_thin_data_is_not_called():
    assert classify(stats(50, n=3), None, None)["status"] == STATUS_LOW_COVERAGE
    low = classify(stats(50, cov=5), None, None)
    assert low["status"] == STATUS_LOW_COVERAGE and "mean depth 5.0x" in low["notes"][0]
    assert classify(None, None, None)["status"] == STATUS_NO_DATA


def test_thresholds_are_configurable_and_reported():
    strict = Thresholds(loss_below=40.0)
    assert classify(stats(35), None, None, strict)["status"] == STATUS_LOSS
    assert classify(stats(35), None, None)["status"] == STATUS_IMPRINTED
    assert Thresholds().as_dict()["allelic_min_diff"] == 40.0


def test_region_result_carries_identity_and_the_call():
    region = {"id": "x", "name": "X:TSS-DMR", "chrom": "chr1", "category": "X", "primary_diagnostic": False,
              "methylated_allele": "maternal", "disorders": [], "disorder_notes": []}
    r = region_result(region, _cpgs(20, 50), _cpgs(20, 95), _cpgs(20, 3))
    assert r["status"] == STATUS_IMPRINTED and r["interpretable"] is False   # category X is shown, not interpreted
    assert r["hap1"]["methylation"] > 90 and r["combined"]["n_cpgs"] == 20


# ---- the region panel ---------------------------------------------------------------------------------

def test_panel_is_the_published_hg38_table():
    panel = service.panel()
    regions = panel["regions"]
    assert panel["source"]["pmcid"] == "PMC11295890" and panel["source"]["genome_build"] == "GRCh38"
    assert len(regions) == 48 and len({r["id"] for r in regions}) == 48
    assert {r["category"] for r in regions} == {"CA", "NC", "X"}
    for r in regions:
        assert r["chrom"].startswith("chr") and r["start"] < r["end"], r["name"]
        if "core" in r:
            assert r["core"]["chrom"] == r["chrom"] and r["core"]["start"] >= r["start"] - 1000


def test_every_disorder_rule_points_at_a_clinically_associated_dmr_in_the_panel():
    by_id = {r["id"]: r for r in service.panel()["regions"]}
    for rule in disorders.DISORDER_RULES:
        for dmr_id, direction, _ in rule["triggers"]:
            assert dmr_id in by_id, (rule["id"], dmr_id)
            assert by_id[dmr_id]["category"] == "CA", (rule["id"], dmr_id)
            assert direction in (STATUS_LOSS, STATUS_GAIN)


def test_the_conflicting_ig_dmr_coordinates_are_recorded_not_hidden():
    ig = next(r for r in service.panel()["regions"] if r["id"] == "meg3-dlk1-ig-dmr")
    assert "100,811,001" in ig["coordinate_notes"]


# ---- disorder patterns and MLID ----------------------------------------------------------------------

def _res(region_id, status, chrom="chr1", category="CA", pct=50.0):
    return {"id": region_id, "name": region_id.upper(), "chrom": chrom, "category": category, "status": status,
            "combined": {"methylation": pct, "n_cpgs": 20, "mean_coverage": 30.0}}


def _normal():
    return {r["id"]: _res(r["id"], STATUS_IMPRINTED, r["chrom"], r["category"]) for r in service.panel()["regions"]}


def _by_id(summaries):
    return {d["id"]: d for d in summaries}


def test_a_normal_genome_shows_no_disorder_pattern():
    out = _by_id(disorders.evaluate_disorders(_normal()))
    assert {d["status"] for d in out.values()} == {"not_observed"}


def test_snurf_gain_is_the_pws_pattern_and_loss_the_as_pattern():
    r = _normal()
    r["snurf-tss-dmr"]["status"] = STATUS_GAIN
    out = _by_id(disorders.evaluate_disorders(r))
    assert out["pws"]["status"] == "pattern_present" and out["as"]["status"] == "not_observed"
    r["snurf-tss-dmr"]["status"] = STATUS_LOSS
    out = _by_id(disorders.evaluate_disorders(r))
    assert out["as"]["status"] == "pattern_present" and out["pws"]["status"] == "not_observed"


def test_bws_and_srs_use_the_right_direction_at_ic1_and_ic2():
    r = _normal()
    r["kcnq1ot1-tss-dmr"]["status"] = STATUS_LOSS              # IC2 loss -> BWS
    assert _by_id(disorders.evaluate_disorders(r))["bws"]["status"] == "pattern_present"
    r = _normal()
    r["h19-igf2-ig-dmr-tss-dmr"]["status"] = STATUS_LOSS       # IC1 loss -> SRS, not BWS
    out = _by_id(disorders.evaluate_disorders(r))
    assert out["srs"]["status"] == "pattern_present" and out["bws"]["status"] == "not_observed"
    r = _normal()
    r["h19-igf2-ig-dmr-tss-dmr"]["status"] = STATUS_GAIN       # IC1 gain -> BWS
    assert _by_id(disorders.evaluate_disorders(r))["bws"]["status"] == "pattern_present"


def test_supportive_loci_alone_do_not_establish_a_pattern():
    r = _normal()
    r["mest-alt-tss-dmr"]["status"] = STATUS_GAIN              # maternal UPD7 pattern, IC1 normal
    out = _by_id(disorders.evaluate_disorders(r))["srs"]
    assert out["status"] == "supportive_only"


def test_missing_data_at_the_primary_locus_is_insufficient_not_normal():
    r = _normal()
    r["plagl1-alt-tss-dmr"]["status"] = STATUS_LOW_COVERAGE
    assert _by_id(disorders.evaluate_disorders(r))["tndm"]["status"] == "insufficient_data"
    assert _by_id(disorders.evaluate_disorders({}))["pws"]["status"] == "insufficient_data"


def test_the_14q32_disorders_have_opposite_directions():
    r = _normal()
    r["meg3-tss-dmr"]["status"] = STATUS_LOSS
    out = _by_id(disorders.evaluate_disorders(r))
    assert out["ts14"]["status"] == "pattern_present" and out["kos14"]["status"] == "not_observed"
    r["meg3-tss-dmr"]["status"] = STATUS_GAIN
    out = _by_id(disorders.evaluate_disorders(r))
    assert out["kos14"]["status"] == "pattern_present" and out["ts14"]["status"] == "not_observed"


def test_mlid_counts_chromosomes_not_dmrs():
    r = _normal()
    r["h19-igf2-ig-dmr-tss-dmr"]["status"] = STATUS_GAIN        # 11p15
    r["kcnq1ot1-tss-dmr"]["status"] = STATUS_LOSS               # also 11p15: one event (e.g. paternal UPD11)
    assert disorders.evaluate_mlid(r)["possible"] is False
    r["snurf-tss-dmr"]["status"] = STATUS_LOSS                  # 15q: a second chromosome
    m = disorders.evaluate_mlid(r)
    assert m["possible"] is True and m["chromosomes"] == ["chr11", "chr15"] and len(m["aberrant"]) == 3


def test_mlid_ignores_polymorphic_dmrs_and_treats_non_clinical_ones_as_research():
    r = _normal()
    r["snurf-tss-dmr"]["status"] = STATUS_LOSS
    for rid in ("gpr1-as-tss-dmr", "ppiel-ex1-dmr"):              # an X and an NC region
        r[rid]["status"] = STATUS_GAIN
    assert disorders.evaluate_mlid(r)["possible"] is False
    research = disorders.evaluate_mlid(r, clinical_only=False)
    assert research["possible"] is True and "research" in research["scope"]


# ---- the whole overview on a synthetic genome ---------------------------------------------------------

def _fake_genome(monkeypatch, override=None):
    """Patch cpg.fetch so every window looks like a normal imprint (hap1 95 / hap2 3), except `override`:
    {region id: (combined %, hap1 %, hap2 %)}."""
    override = override or {}
    ids_by_window = {}
    for r in service.panel()["regions"]:
        for w in (r["start"], r.get("core", {}).get("start")):
            if w is not None:
                ids_by_window[(r["chrom"], w)] = r["id"]

    def fake_fetch(bed, chrom, start, end):
        rid = ids_by_window.get((chrom, start))
        track = bed.name.split(".")[0]
        pct = override.get(rid, (48, 95, 3))
        value = {"combined": pct[0], "hap1": pct[1], "hap2": pct[2]}[track]
        return _cpgs(25, value, 30, start)

    monkeypatch.setattr(cpg, "fetch", fake_fetch)
    monkeypatch.setattr(cpg, "header", lambda bed: {"pb_cpg_tools_version": "3.0.0", "pileup_mode": "model", "modsites_mode": "denovo", "min_coverage": 4})


@pytest.fixture
def synthetic_sample(tmp_path):
    paths = {}
    for t in ("combined", "hap1", "hap2"):
        p = tmp_path / f"{t}.bed.gz"
        p.write_bytes(b"")
        paths[t] = p
    service._cache.clear()
    service._cpg_cache.clear()
    return samples.Sample(key="SYN", label="Synthetic", paths=paths)


def test_overview_of_a_normal_synthetic_genome(monkeypatch, synthetic_sample):
    _fake_genome(monkeypatch)
    o = service.overview(synthetic_sample)
    assert o["summary"]["regions"] == 48 and o["summary"]["by_status"] == {"imprinted": 48}
    assert {d["status"] for d in o["disorders"]} == {"not_observed"} and o["mlid"]["possible"] is False
    assert o["provenance"]["pb_cpg_tools_version"] == "3.0.0" and o["thresholds"]["loss_below"] == 20.0
    assert o["window"] == "core"
    snurf = next(r for r in o["regions"] if r["id"] == "snurf-tss-dmr")
    assert snurf["window"] == "core" and snurf["interval"]["start"] == 24954888 and snurf["dmr"]["start"] == 24954855


def test_overview_finds_a_pws_like_pattern_with_a_second_chromosome_for_mlid(monkeypatch, synthetic_sample):
    _fake_genome(monkeypatch, {"snurf-tss-dmr": (96, 96, 97), "plagl1-alt-tss-dmr": (4, 5, 3)})
    o = service.overview(synthetic_sample)
    d = {x["id"]: x for x in o["disorders"]}
    assert d["pws"]["status"] == "pattern_present" and d["tndm"]["status"] == "pattern_present"
    assert o["mlid"]["possible"] is True and o["mlid"]["chromosomes"] == ["chr15", "chr6"]
    assert o["summary"]["by_status"]["gain"] == 1 and o["summary"]["by_status"]["loss"] == 1


def test_changing_a_threshold_reclassifies_without_rereading_the_beds(monkeypatch, synthetic_sample):
    _fake_genome(monkeypatch)
    first = service.overview(synthetic_sample)
    reads = []
    monkeypatch.setattr(cpg, "fetch", lambda *a: reads.append(a) or [])
    stricter = service.overview(synthetic_sample, th=Thresholds(allelic_min_diff=99))
    assert reads == []                                           # CpGs came from the per-window cache
    assert first["summary"]["by_status"] == {"imprinted": 48}
    assert stricter["summary"]["by_status"] != first["summary"]["by_status"]
    assert stricter["thresholds"]["allelic_min_diff"] == 99.0


def test_window_can_be_switched_to_the_published_dmr(monkeypatch, synthetic_sample):
    _fake_genome(monkeypatch)
    o = service.overview(synthetic_sample, window="dmr")
    assert {r["window"] for r in o["regions"]} == {"dmr"}


def test_region_track_returns_cpg_level_points(monkeypatch, synthetic_sample):
    _fake_genome(monkeypatch)
    t = service.region_track(synthetic_sample, "snurf-tss-dmr", flank=500)
    assert t["name"] == "SNURF:TSS-DMR" and set(t["tracks"]) == {"combined", "hap1", "hap2"}
    assert t["analysis_window"]["which"] == "core" and t["view"]["start"] == 24954855 - 500
    assert service.region_track(synthetic_sample, "nope") is None


# ---- sample registry and the real file reader ---------------------------------------------------------

def test_registry_needs_a_combined_bed_and_its_index(tmp_path, monkeypatch):
    for name in ("S1.combined.bed.gz", "S1.combined.bed.gz.tbi", "S1.hap1.bed.gz", "S1.hap1.bed.gz.tbi", "S2.combined.bed.gz"):
        (tmp_path / name).write_bytes(b"")
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"samples": {
        "S1": {"dir": str(tmp_path), "prefix": "S1", "label": "Sample one"},
        "S2": {"dir": str(tmp_path), "prefix": "S2"},          # no .tbi: not registered
        "S3": {"dir": str(tmp_path), "prefix": "S3"},          # no files: not registered
    }}))
    monkeypatch.setenv("EPIGENETICS_CONFIG", str(cfg))
    got = samples.load()
    assert list(got) == ["S1"] and got["S1"].label == "Sample one"
    assert sorted(got["S1"].paths) == ["combined", "hap1"]      # hap2 absent: dropped, sample kept


def test_registry_without_a_config_is_empty(tmp_path, monkeypatch):
    monkeypatch.setenv("EPIGENETICS_CONFIG", str(tmp_path / "missing.json"))
    assert samples.load() == {}


@pytest.mark.skipif(not (shutil.which("tabix") and shutil.which("bgzip")), reason="htslib tabix/bgzip not installed")
def test_cpg_reader_on_a_real_tabix_indexed_bed(tmp_path):
    bed = tmp_path / "x.combined.bed"
    bed.write_text(
        "##pb-cpg-tools-version=3.0.0\n"
        "##cmdline=aligned_bam_to_cpg_scores --threads 4 --min-mapq 1 --min-coverage 4\n"
        "##pileup-mode=model\n##modsites-mode=denovo\n"
        "chr1\t100\t101\t80.0\tTotal\t10\t8\t2\t80.0\n"
        "chr1\t200\t201\t10.0\tTotal\t20\t2\t18\t10.0\n"
        "chr2\t50\t51\t55.5\tTotal\t30\t17\t13\t55.0\n"
    )
    subprocess.run(["bgzip", str(bed)], check=True)
    gz = tmp_path / "x.combined.bed.gz"
    subprocess.run(["tabix", "-p", "bed", str(gz)], check=True)

    assert [c.pos for c in cpg.fetch(gz, "chr1", 0, 1000)] == [100, 200]
    assert [c.pos for c in cpg.fetch(gz, "chr1", 150, 1000)] == [200]
    got = cpg.fetch(gz, "chr2", 0, 100)[0]
    assert (got.mod_score, got.coverage, got.est_mod, got.est_unmod) == (55.5, 30, 17, 13)
    assert cpg.fetch(gz, "chr1", 500, 400) == []
    assert cpg.header(gz) == {"pb_cpg_tools_version": "3.0.0", "pileup_mode": "model", "modsites_mode": "denovo", "min_coverage": 4}
    assert set(cpg.contigs(gz)) == {"chr1", "chr2"}
