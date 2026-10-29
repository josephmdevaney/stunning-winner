"""Tests for the Kivvi suite: output parsing, FSHD1 interpretation, the
auto-discovering sample registry and the API routes.

All inputs are small synthetic Kivvi-shaped files built in tmp_path -- no real
per-participant output is read or committed.
"""

import gzip
import json

import pytest
from fastapi.testclient import TestClient

from kivvi import parse, samples
from repeatscope.api.main import app

client = TestClient(app)


# --- parsing helpers --------------------------------------------------------

def test_unit_ids_drops_flanks():
    assert parse.unit_ids("LeftFlank-1-14-14-RightFlank") == ["1", "14", "14"]
    # Partial alleles have no left flank.
    assert parse.unit_ids("41-41-47-RightFlank") == ["41", "41", "47"]


def test_parse_allele_size_lower_bound():
    assert parse.parse_allele_size("28") == (28, False)
    assert parse.parse_allele_size(">14") == (14, True)
    # Kivvi v1.1.0 writes ">=N" for merged and partial alleles (the docs only show ">N").
    assert parse.parse_allele_size(">=12") == (12, True)
    assert parse.parse_allele_size(" >=10 ") == (10, True)


@pytest.mark.parametrize("raw", ["", "=12", "~12", "12+", ">", ">=", "1.5", "abc"])
def test_parse_allele_size_refuses_unknown_formats(raw):
    with pytest.raises(ValueError, match="allele_size"):
        parse.parse_allele_size(raw)


def test_split_chromosome():
    assert parse.split_chromosome("chr4:Group1.1") == ("chr4", "Group1.1")
    assert parse.split_chromosome("unknown") == ("unknown", None)


# --- FSHD1 interpretation ---------------------------------------------------

def _classify(size, chromosome="chr4", group="Group1.1", distal="qAIntactPolyA", lower=False):
    return parse.classify_d4z4_allele(size, lower, chromosome, group, distal)["status"]


def test_contracted_chr4_qa_intact_is_fshd1_range():
    assert _classify(7) == "pathogenic"


def test_contraction_boundary_is_ten_units():
    assert _classify(10) == "pathogenic"
    assert _classify(11) == "normal"


def test_contracted_4a166_group_is_not_pathogenic():
    assert _classify(7, group="Group2.1") == "normal"


@pytest.mark.parametrize("distal", ["qB", "qADisruptedPolyA"])
def test_contracted_non_intact_polya_is_not_pathogenic(distal):
    assert _classify(7, distal=distal) == "normal"


@pytest.mark.parametrize("chromosome,group", [("chr10", "Group2.1"), ("unknown", None)])
def test_contracted_qa_intact_off_chr4_needs_review(chromosome, group):
    assert _classify(7, chromosome=chromosome, group=group) == "review"


def test_partial_allele_lower_bound():
    # A lower bound above 10 rules contraction out; one at or below cannot.
    assert _classify(14, lower=True) == "normal"
    assert _classify(8, lower=True) == "indeterminate"


# --- shaping ----------------------------------------------------------------

def _allele_info(size, chromosome="chr4:Group1.1", distal="qAIntactPolyA", allele_type="assembled", name=None):
    n = int(str(size).lstrip(">="))
    return {
        "allele_name": name or "LeftFlank-" + "-".join(str(i) for i in range(1, n + 1)) + "-RightFlank",
        "chromosome": chromosome,
        "distal_haplotype": distal,
        "allele_type": allele_type,
        "allele_size": str(size),
        "methylation": 0.8,
        "ending_in_qal": False,
    }


def _d4z4_json(infos):
    return {"allele_info": infos, "additional": {"methylation_all_sites": 0.8, "qc_metrics": {"median_read_length": 18000.0}}}


def test_summary_pathogenic_outranks_everything():
    data = _d4z4_json([_allele_info(7), _allele_info(30), _allele_info(40, chromosome="chr10:Group2.1"), _allele_info(">9", chromosome="unknown", allele_type="partial")])
    summary = parse.shape_d4z4(data)["summary"]
    assert summary["status"] == "pathogenic"
    assert (summary["complete_alleles"], summary["partial_alleles"]) == (3, 1)


def test_summary_flags_unresolved_alleles():
    shaped = parse.shape_d4z4(_d4z4_json([_allele_info(30), _allele_info(40)]))
    assert shaped["summary"]["status"] == "normal"
    assert "2 of 4" in shaped["summary"]["notes"][0]


def test_merged_and_partial_alleles_with_gte_sizes():
    # The shape of a real v1.1.0 run: one assembled 4-unit chr4 allele plus merged / partial alleles sized ">=N".
    data = _d4z4_json(
        [
            _allele_info(4),
            _allele_info(">=12", allele_type="merged"),
            _allele_info(">=10", chromosome="chr10:Group2.1", distal="qADisruptedPolyA", allele_type="partial"),
        ]
    )
    alleles = parse.shape_d4z4(data)["alleles"]
    assert [(a["size"], a["size_label"], a["size_is_lower_bound"]) for a in alleles] == [(4, "4", False), (12, "≥12", True), (10, "≥10", True)]
    assert [a["allele_type"] for a in alleles] == ["assembled", "merged", "partial"]
    assert alleles[0]["interpretation"]["status"] == "pathogenic"
    assert alleles[1]["interpretation"]["status"] == "normal"  # at least 12 units: above the 1-10 range


def test_summary_with_no_alleles_is_indeterminate():
    shaped = parse.shape_d4z4(_d4z4_json([]))
    assert shaped["summary"]["status"] == "indeterminate"


def test_per_unit_methylation_attached_by_allele_name():
    info = _allele_info(3)
    data = _d4z4_json([info])
    data["methylation"] = {"all_allele_ends": {"median_methylation_per_unit": {info["allele_name"]: "0.5,0.6,0.7"}}}
    assert parse.shape_d4z4(data)["alleles"][0]["per_unit_methylation"] == [0.5, 0.6, 0.7]


def test_kiv2_copy_number_from_complete_alleles():
    short, long_ = "LeftFlank-1-2-3-RightFlank", "LeftFlank-1-2-3-4-5-RightFlank"
    data = {"allele_cn": "3/5", "complete_alleles": [short, long_], "supporting_reads": {short: ["r1"], long_: ["r1", "r2"]}}
    shaped = parse.shape_kiv2(data)
    assert [(a["copy_number"], a["supporting_reads"]) for a in shaped["alleles"]] == [(3, 1), (5, 2)]
    assert shaped["shortest_copy_number"] == 3
    assert shaped["notes"] == []


def test_kiv2_single_allele_is_noted():
    shaped = parse.shape_kiv2({"complete_alleles": ["LeftFlank-1-2-RightFlank"]})
    assert "1 of 2" in shaped["notes"][0]


def test_kiv2_no_alleles():
    shaped = parse.shape_kiv2({"allele_cn": "", "complete_alleles": []})
    assert shaped["alleles"] == [] and shaped["shortest_copy_number"] is None


# --- FSHD2 evidence ---------------------------------------------------------

def test_permissive_means_chr4_qa_intact_polya_outside_4a166():
    assert parse.is_permissive("chr4", "Group1.1", "qAIntactPolyA")
    assert parse.is_permissive("chr4", None, "qAIntactPolyA")
    assert not parse.is_permissive("chr4", "Group2.1", "qAIntactPolyA")      # 4A166
    assert not parse.is_permissive("chr10", "Group1.1", "qAIntactPolyA")
    assert not parse.is_permissive("chr4", "Group1.1", "qADisruptedPolyA")
    assert not parse.is_permissive("chr4", "Group1.1", "qB")


@pytest.mark.parametrize(
    "size, lower, expected",
    [(1, False, "contracted"), (10, False, "contracted"), (11, False, "fshd2_range"), (20, False, "fshd2_range"),
     (21, False, "above"), (21, True, "above"), (20, True, "unresolved"), (8, True, "unresolved")],
)
def test_size_range_boundaries(size, lower, expected):
    assert parse.size_range(size, lower) == expected


def _meth_allele(size, methylation, **kw):
    info = _allele_info(size, **kw)
    info["methylation"] = methylation
    return info


def test_fshd2_evidence_lays_out_each_allele_without_making_a_call():
    data = _d4z4_json([
        _meth_allele(14, 0.31),                                          # permissive, 11-20 units
        _meth_allele(9, 0.28),                                           # permissive, contracted
        _meth_allele(35, 0.40, chromosome="chr4:Group2.1"),             # 4A166: not permissive
        _meth_allele(">12", 0.35, chromosome="unknown", allele_type="partial", distal="qB"),
    ])
    data["additional"]["methylation_all_sites"] = 0.33
    fshd2 = parse.shape_d4z4(data)["fshd2"]
    assert fshd2["sample_methylation"] == 0.33
    assert [(a["index"], a["permissive"], a["size_range"]) for a in fshd2["alleles"]] == [
        (1, True, "fshd2_range"), (2, True, "contracted"), (3, False, "above"), (4, False, "unresolved"),
    ]
    assert fshd2["alleles"][3]["size_label"] == "≥12" and fshd2["alleles"][3]["methylation"] == 0.35
    assert fshd2["summary"] == {
        "permissive_alleles": 2, "permissive_contracted": 1, "permissive_fshd2_range": 1,
        "alleles_with_methylation": 4, "lowest_methylation": 0.28, "highest_methylation": 0.40,
    }
    # evidence only: no verdict anywhere in the block, and the FSHD1 summary is unchanged by it
    assert not {"status", "label", "call"} & set(fshd2) and not {"status", "label", "call"} & set(fshd2["summary"])
    assert parse.shape_d4z4(data)["summary"]["status"] == "pathogenic"


def test_fshd2_evidence_tolerates_alleles_without_methylation_and_no_alleles():
    info = _allele_info(14)
    info["methylation"] = None
    fshd2 = parse.shape_d4z4(_d4z4_json([info]))["fshd2"]
    assert fshd2["summary"]["alleles_with_methylation"] == 0 and fshd2["summary"]["lowest_methylation"] is None
    empty = parse.shape_d4z4(_d4z4_json([]))["fshd2"]
    assert empty["alleles"] == [] and empty["summary"]["permissive_alleles"] == 0


# --- VCF --------------------------------------------------------------------

VCF = "\n".join(
    [
        "##fileformat=VCFv4.2",
        "##kivviCommand=kivvi -b x.bam d4z4",
        "##kivviVersion=1.1.0",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tallele1\tallele2",
        "chr6\t160613775\t.\tG\tC\t.\tPASS\tRU=31:2.3:31:31,30:2.4:30:15\tGT\t0\t1",
        "chr6\t160613800\t.\tCG\tC\t.\tPASS\tRU=34:Unknown:1:1\tGT\t0\t0",
    ]
)


def test_parse_vcf():
    vcf = parse.parse_vcf(VCF)
    assert vcf["version"] == "1.1.0"
    assert vcf["allele_columns"] == ["allele1", "allele2"]
    first, second = vcf["records"]
    assert first["genotypes"] == ["0", "1"] and first["on_assembled_allele"] is True
    assert first["n_units"] == 2
    assert first["min_vaf"] == pytest.approx(0.5)
    assert first["units"][0] == {"unit": "31", "where": "2.3", "depth": 31, "reads": 31}
    # Variant on an unassembled unit: no allele carries it.
    assert second["on_assembled_allele"] is False


# --- registry + API ---------------------------------------------------------

@pytest.fixture
def results_root(tmp_path, monkeypatch):
    """A results root with one sample holding D4Z4 (with SVG) and no KIV2."""
    d4z4 = tmp_path / "S1" / "d4z4"
    d4z4.mkdir(parents=True)
    (d4z4 / "S1.kivvi.d4z4.json").write_text(json.dumps(_d4z4_json([_allele_info(7), _allele_info(30)])))
    (d4z4 / "S1.kivvi.d4z4.vcf").write_text(VCF)
    (d4z4 / "S1.kivvi.d4z4.svg").write_text("<svg xmlns='http://www.w3.org/2000/svg'/>")
    # A directory without Kivvi output is not a sample.
    (tmp_path / "scratch").mkdir()
    monkeypatch.setattr(samples, "RESULTS_ROOT", tmp_path)
    return tmp_path


def test_samples_are_discovered_from_results_root(results_root):
    body = client.get("/api/kivvi/samples").json()
    assert body == {"samples": [{"key": "S1", "label": "S1"}], "default": "S1"}


def test_no_results_root_yields_no_samples(tmp_path, monkeypatch):
    monkeypatch.setattr(samples, "RESULTS_ROOT", tmp_path / "missing")
    assert client.get("/api/kivvi/samples").json() == {"samples": [], "default": ""}


def test_overview_reports_missing_target_as_null(results_root):
    body = client.get("/api/kivvi/overview").json()
    assert body["targets"] == {"d4z4": True, "kiv2": False}
    assert body["kiv2"] is None
    assert body["d4z4"]["summary"]["status"] == "pathogenic"
    assert body["kivvi_version"] == "1.1.0"
    # No genome depth in the synthetic JSON -> bamlet note.
    assert body["input_note"]


def test_target_detail_includes_variants_and_svg_flag(results_root):
    body = client.get("/api/kivvi/target/d4z4", params={"sample": "S1"}).json()
    assert body["variants"]["n_total"] == 2 and body["variants"]["n_on_assembled"] == 1
    assert body["svg_available"] is True


def test_variant_columns_are_labelled_by_allele_in_vcf_order(results_root):
    # VCF columns follow complete_alleles order, which differs from allele_info order.
    big, small = _allele_info(30), _allele_info(7)
    data = _d4z4_json([small, big])  # allele_info order: small first
    data["complete_alleles"] = [big["allele_name"], small["allele_name"]]  # VCF order: big first
    (results_root / "S1" / "d4z4" / "S1.kivvi.d4z4.json").write_text(json.dumps(data))
    labels = client.get("/api/kivvi/target/d4z4").json()["variants"]["allele_labels"]
    assert labels == ["A1 · 30 · chr4", "A2 · 7 · chr4"]


def test_variant_labels_fall_back_to_vcf_columns(results_root):
    # The fixture JSON has no complete_alleles, so the VCF's own column names are kept.
    labels = client.get("/api/kivvi/target/d4z4").json()["variants"]["allele_labels"]
    assert labels == ["allele1", "allele2"]


def _swap_vcf_for_gz(results_root, text=VCF):
    """Replace the fixture's plain VCF with a bgzip-style one: several concatenated gzip members."""
    target = results_root / "S1" / "d4z4"
    (target / "S1.kivvi.d4z4.vcf").unlink()
    data = text.encode()
    mid = len(data) // 2
    (target / "S1.kivvi.d4z4.vcf.gz").write_bytes(gzip.compress(data[:mid]) + gzip.compress(data[mid:]))
    return target


def test_gzipped_vcf_is_read_directly(results_root):
    _swap_vcf_for_gz(results_root)
    body = client.get("/api/kivvi/target/d4z4", params={"sample": "S1"}).json()
    assert body["variants"]["n_total"] == 2 and body["variants"]["n_on_assembled"] == 1
    assert body["variants"]["allele_columns"] == ["allele1", "allele2"]
    assert client.get("/api/kivvi/overview").json()["kivvi_version"] == "1.1.0"


def test_plain_vcf_wins_when_both_exist(results_root):
    target = _swap_vcf_for_gz(results_root)
    plain = VCF.rsplit("\n", 1)[0]  # drop the second record
    (target / "S1.kivvi.d4z4.vcf").write_text(plain)
    assert client.get("/api/kivvi/target/d4z4").json()["variants"]["n_total"] == 1


def test_no_vcf_yields_no_variants(results_root):
    (results_root / "S1" / "d4z4" / "S1.kivvi.d4z4.vcf").unlink()
    variants = client.get("/api/kivvi/target/d4z4").json()["variants"]
    assert variants["n_total"] == 0 and variants["records"] == []


def test_svg_endpoint(results_root):
    ok = client.get("/api/kivvi/target/d4z4/svg")
    assert ok.status_code == 200 and ok.headers["content-type"] == "image/svg+xml"
    assert client.get("/api/kivvi/target/kiv2/svg").status_code == 404


def test_unknown_sample_and_target_are_404(results_root):
    assert client.get("/api/kivvi/overview", params={"sample": "nope"}).status_code == 404
    assert client.get("/api/kivvi/target/bogus").status_code == 404
    assert client.get("/api/kivvi/target/kiv2").status_code == 404  # known target, no output


def test_qc_lists_only_targets_with_output(results_root):
    qc = client.get("/api/kivvi/qc").json()
    assert [t["target"] for t in qc["target_qc"]] == ["d4z4"]
    assert qc["target_qc"][0]["alleles_resolved"] == 2
    assert qc["target_qc"][0]["alleles_expected"] == 4
    assert qc["target_qc"][0]["median_read_length"] == 18000.0


def test_target_detail_carries_the_fshd2_evidence(results_root):
    body = client.get("/api/kivvi/target/d4z4", params={"sample": "S1"}).json()
    assert body["fshd2"]["summary"]["permissive_alleles"] == 2
    assert body["fshd2"]["sample_methylation"] == 0.8
    assert [a["size_label"] for a in body["fshd2"]["alleles"]] == ["7", "30"]
