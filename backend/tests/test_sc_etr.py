"""Tests for the TenK10K sc-eTR annotation of the STRchive disease panel.

All inputs are small synthetic files in the formats of the paper's released tables
(Zenodo 10.5281/zenodo.15009519); none of the paper's data is in the repo (it is
CC BY-NC-ND 4.0)."""

import importlib.util
import json
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from repeatscope.api import main
from repeatscope.api.main import app
from repeatscope.motifs import is_same_repeat_unit
from repeatscope.sc_etr import SCHEMA_VERSION, load_sc_etr_annotations, locus_annotation

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "real_data" / "build_sc_etr_disease_annotation.py"
spec = importlib.util.spec_from_file_location("build_sc_etr_disease_annotation", SCRIPT)
build = importlib.util.module_from_spec(spec)
spec.loader.exec_module(build)

QVAL_HEADER = "gene_name\tgene_level_pval\tchr\tpos\tcoeff_meta_fixed\tse_meta_fixed\tpval_meta_fixed\tmotif\tqval"


def qval_row(gene, q, chroms, pos, motifs, coeffs, ses=None, pvals=None):
    n = len(chroms)
    ses = ses or [0.1] * n
    pvals = pvals or [1e-6] * n
    return "\t".join([gene, "0.0", repr(chroms), repr(pos), repr(coeffs), repr(ses), repr(pvals), repr(motifs), str(q)])


def write_qval(path, rows):
    path.write_text(QVAL_HEADER + "\n" + "\n".join(rows) + "\n")


@pytest.fixture
def qval_dir(tmp_path):
    d = tmp_path / "qval"
    d.mkdir()
    write_qval(
        d / "CD14_Mono_qval.tsv",
        [
            qval_row("ENSG1", 0.001, ["chr1"], [100], ["CTG"], [0.5], pvals=[1e-8]),
            # not significant: must be ignored
            qval_row("ENSG2", 0.5, ["chr1"], [100], ["CTG"], [0.9]),
            # tied leads for one gene
            qval_row("ENSG3", 0.01, ["chr1", "chr1"], [5000, 100], ["AC", "CTG"], [0.3, 0.4]),
            qval_row("ENSG4", 0.02, ["chr1"], [5010], ["A"], [-0.7]),
            qval_row("ENSG5", 0.03, ["chr2"], [100], ["GT"], [0.2]),
        ],
    )
    write_qval(d / "NK_qval.tsv", [qval_row("ENSG1", 0.02, ["chr1"], [100], ["CTG"], [-0.2])])
    # the permutation null: significant-looking rows that must not count
    write_qval(d / "CD4_TCM_permuted_qval.tsv", [qval_row("ENSG9", 1e-9, ["chr1"], [100], ["CTG"], [9.0])])
    return d


@pytest.fixture
def catalog_json(tmp_path):
    def locus(region, structure, locus_id, **extra):
        return {
            "ReferenceRegion": region,
            "LocusStructure": structure,
            "LocusId": locus_id,
            "Source": "synthetic",
            "ManeGeneName": None,
            "ManeGeneId": None,
            **extra,
        }

    loci = [
        locus("chr1:100-130", "(CTG)*", "1-101-130-CTG", VariantType="Repeat", ManeGeneName="GENEA", ManeGeneId="ENSG1.5"),
        # a multi-repeat locus: list-valued fields, one tested variant per component
        locus(
            ["chr1:5000-5010", "chr1:5010-5020"],
            "(AC)*(A)*",
            "1-5001-5010-AC",
            VariantType=["Repeat", "Repeat"],
            VariantId=["1-5001-5010-AC", "1-5001-5010-AC_A"],
            InterruptionBaseCount=[0, 0],
        ),
        locus("chr2:100-140", "(GT)*", "2-101-140-GT", VariantType="Repeat"),
    ]
    path = tmp_path / "catalog.json"
    path.write_text(json.dumps(loci, indent=2))
    return path


@pytest.fixture
def strchive_json(tmp_path):
    def locus(locus_id, gene, chrom, start, stop, motifs):
        return {
            "id": locus_id, "gene": gene, "chrom": chrom, "start_hg38": start, "stop_hg38": stop,
            "reference_motif_reference_orientation": motifs,
        }

    loci = [
        locus("TEST_A", "GENEA", "chr1", 100, 130, ["CAG"]),
        locus("TEST_B", "GENEB", "chr1", 5010, 5020, ["A"]),
        locus("TEST_FAR", "GENEC", "chr1", 6500, 6530, ["GA"]),
        locus("TEST_NONE", "GENED", "chr3", 100, 130, ["CCG"]),
    ]
    path = tmp_path / "strchive.json"
    path.write_text(json.dumps(loci))
    return path


@pytest.fixture
def built(qval_dir, catalog_json, strchive_json):
    leads, n_sig = build.read_leads(qval_dir)
    return build.build_annotations(
        leads, build.iter_catalog_repeats(catalog_json), build.load_panel(strchive_json), n_sig, built="2026-01-01"
    )


# ---------------------------------------------------------------- motifs
def test_same_repeat_unit_ignores_rotation_and_strand():
    for other in ("CAG", "AGC", "GCA", "CTG", "TGC", "GCT"):
        assert is_same_repeat_unit("CAG", other)
    assert not is_same_repeat_unit("CAG", "CCG")
    assert not is_same_repeat_unit("CAG", "CAGCAG")
    assert is_same_repeat_unit("gcc", "GGC")


# ---------------------------------------------------------------- reading the released tables
def test_read_leads_applies_cutoff_skips_permuted_and_splits_tied_leads(qval_dir):
    leads, n_sig = build.read_leads(qval_dir)
    assert n_sig == 5  # ENSG1, ENSG3, ENSG4, ENSG5 in CD14_Mono and ENSG1 in NK
    assert set(leads) == {
        ("chr1", 100, "CTG"), ("chr1", 5000, "AC"), ("chr1", 5010, "A"), ("chr2", 100, "GT"),
    }
    # the tied row gives one lead to ENSG3 at each of its two repeats
    assert "ENSG3" in leads[("chr1", 100, "CTG")] and "ENSG3" in leads[("chr1", 5000, "AC")]
    # ENSG2 (q=0.5) and the permutation null (ENSG9) are not leads
    assert "ENSG2" not in leads[("chr1", 100, "CTG")]
    assert all("ENSG9" not in by_gene for by_gene in leads.values())
    assert set(leads[("chr1", 100, "CTG")]["ENSG1"]) == {"CD14_Mono", "NK"}


def test_catalog_streaming_splits_multi_repeat_loci_into_components(catalog_json):
    repeats = list(build.iter_catalog_repeats(catalog_json))
    assert [(r["chrom"], r["start"], r["end"], r["motif"]) for r in repeats] == [
        ("chr1", 100, 130, "CTG"),
        ("chr1", 5000, 5010, "AC"),
        ("chr1", 5010, 5020, "A"),
        ("chr2", 100, 140, "GT"),
    ]
    # the gene id loses its version, and a null gene is empty
    assert (repeats[0]["gene_id"], repeats[0]["gene"]) == ("ENSG1", "GENEA")
    assert repeats[1]["gene_id"] == "" and repeats[1]["gene"] == ""


# ---------------------------------------------------------------- the join
def test_lead_at_the_disease_repeat(built):
    a = built["loci"]["TEST_A"]
    assert a["polymorphic"] is True
    (lead,) = a["leads"]
    assert (lead["chrom"], lead["start"], lead["stop"], lead["motif"]) == ("chr1", 100, 130, "CTG")
    assert lead["relation"] == "overlaps" and lead["distance_bp"] == 0
    assert lead["same_motif"] is True  # CTG is CAG read from the other strand

    first, second = lead["targets"]  # ordered by best q
    assert (first["gene_id"], first["gene"]) == ("ENSG1", "GENEA")
    assert [(e["cell_type"], e["coeff"], e["qval"]) for e in first["effects"]] == [
        ("CD14_Mono", 0.5, 0.001), ("NK", -0.2, 0.02),
    ]
    assert (second["gene_id"], second["gene"]) == ("ENSG3", None)  # no symbol in the catalog
    assert all(t["gene_id"] != "ENSG9" for t in lead["targets"])


def test_adjacent_repeat_is_kept_but_marked_as_a_different_repeat(built):
    b = built["loci"]["TEST_B"]
    overlapping, adjacent = b["leads"]  # the repeat at the locus sorts first
    assert (overlapping["relation"], overlapping["motif"], overlapping["same_motif"]) == ("overlaps", "A", True)
    # touches the locus (ends where it starts): adjacent at 0 bp, and not the disease motif
    assert (adjacent["relation"], adjacent["distance_bp"], adjacent["motif"]) == ("adjacent", 0, "AC")
    assert adjacent["same_motif"] is False
    # the lead 4.9 kb away at TEST_A's locus is out of range
    assert all(l["start"] >= 5000 for l in b["leads"])


def test_locus_without_a_lead_still_says_whether_it_was_polymorphic(built):
    far = built["loci"]["TEST_FAR"]  # the nearest lead is 1,480 bp away, past the 1 kb limit
    assert far == {"polymorphic": False, "leads": []}
    none = built["loci"]["TEST_NONE"]
    assert none == {"polymorphic": False, "leads": []}


def test_polymorphic_is_independent_of_leads(qval_dir, catalog_json, tmp_path):
    # a catalog repeat sits on the locus, but no lead does
    strchive = tmp_path / "s.json"
    strchive.write_text(json.dumps([{
        "id": "POLY_NO_LEAD", "gene": "G", "chrom": "chr2", "start_hg38": 110, "stop_hg38": 120,
        "reference_motif_reference_orientation": ["GT"],
    }]))
    leads, n_sig = build.read_leads(qval_dir)
    leads.pop(("chr2", 100, "GT"))
    out = build.build_annotations(leads, build.iter_catalog_repeats(catalog_json), build.load_panel(strchive), n_sig)
    assert out["loci"]["POLY_NO_LEAD"] == {"polymorphic": True, "leads": []}


def test_counts_and_source_are_recorded(built):
    src = built["source"]
    assert built["schema_version"] == SCHEMA_VERSION
    assert src["built"] == "2026-01-01" and src["genome_build"] == "GRCh38" and src["lead_only"] is True
    assert src["paper"]["doi"] == "10.1101/2024.11.02.621562" and src["paper"]["license"] == "CC BY-NC-ND 4.0"
    assert src["counts"] == {
        "catalog_repeats": 4,
        "significant_gene_cell_type_rows": 5,
        "lead_repeats": 4,
        "lead_repeats_placed": 4,
        "egenes": 4,
        "panel_loci": 4,
        "panel_loci_polymorphic": 2,
        "panel_loci_with_lead_at_repeat": 2,
        "panel_loci_with_lead_nearby_only": 0,
    }


def test_a_lead_with_no_matching_catalog_repeat_is_left_out_not_guessed(qval_dir, catalog_json, strchive_json):
    leads, n_sig = build.read_leads(qval_dir)
    leads[("chr1", 100, "GGGG")] = leads[("chr1", 100, "CTG")]  # right position, wrong motif
    out = build.build_annotations(leads, build.iter_catalog_repeats(catalog_json), build.load_panel(strchive_json), n_sig)
    assert out["source"]["counts"]["lead_repeats"] == 5
    assert out["source"]["counts"]["lead_repeats_placed"] == 4
    assert all(l["motif"] != "GGGG" for e in out["loci"].values() for l in e["leads"])


# ---------------------------------------------------------------- loading the file
def test_load_roundtrip_and_locus_lookup(built, tmp_path):
    path = tmp_path / "ann.json"
    path.write_text(json.dumps(built))
    data = load_sc_etr_annotations(path)
    ann = locus_annotation(data, "TEST_A")
    assert ann["available"] is True and ann["polymorphic"] is True
    assert ann["leads"] == built["loci"]["TEST_A"]["leads"]
    assert ann["source"]["paper"]["doi"] == "10.1101/2024.11.02.621562"
    assert locus_annotation(data, "NOT_A_LOCUS") is None


@pytest.mark.parametrize(
    "content",
    ['{"schema_version": 99, "loci": {}, "source": {}}', '{"schema_version": 1}', "[]", "{not json"],
)
def test_load_rejects_files_it_does_not_understand(tmp_path, content):
    path = tmp_path / "bad.json"
    path.write_text(content)
    with pytest.raises(ValueError):
        load_sc_etr_annotations(path)


# ---------------------------------------------------------------- API
def _locus(client, locus_id="HD_HTT"):
    resp = client.get(f"/api/locus/{locus_id}", params={"mode": "genome_wide"})
    assert resp.status_code == 200
    return resp.json()


def _annotation_file(tmp_path, loci):
    path = tmp_path / "ann.json"
    path.write_text(json.dumps({
        "schema_version": SCHEMA_VERSION, "source": {"paper": {"doi": "x"}}, "loci": loci,
    }))
    return path


def test_locus_detail_has_no_expression_effects_when_the_file_is_absent(tmp_path, monkeypatch):
    monkeypatch.setattr(main, "SC_ETR_ANNOTATIONS_JSON", str(tmp_path / "missing.json"))
    with TestClient(app) as client:
        assert _locus(client)["expression_effects"] is None


def test_locus_detail_carries_expression_effects(tmp_path, monkeypatch):
    lead = {
        "chrom": "chr4", "start": 3074945, "stop": 3074968, "motif": "CCG", "relation": "adjacent",
        "distance_bp": 12, "same_motif": False,
        "targets": [{"gene_id": "ENSG1", "gene": "HTT", "effects": [
            {"cell_type": "CD8_TEM", "coeff": -0.048, "se": 0.008, "pval": 1.8e-9, "qval": 1e-7}]}],
    }
    path = _annotation_file(tmp_path, {"HD_HTT": {"polymorphic": True, "leads": [lead]}})
    monkeypatch.setattr(main, "SC_ETR_ANNOTATIONS_JSON", str(path))
    with TestClient(app) as client:
        ee = _locus(client)["expression_effects"]
        other = _locus(client, "FTDALS1_C9orf72")["expression_effects"]
    assert ee["available"] is True and ee["polymorphic"] is True and ee["leads"] == [lead]
    assert ee["source"] == {"paper": {"doi": "x"}}
    assert other is None  # the file has no entry for this locus


def test_a_rebuilt_file_is_picked_up_without_a_restart(tmp_path, monkeypatch):
    path = _annotation_file(tmp_path, {"HD_HTT": {"polymorphic": False, "leads": []}})
    monkeypatch.setattr(main, "SC_ETR_ANNOTATIONS_JSON", str(path))
    with TestClient(app) as client:
        assert _locus(client)["expression_effects"]["polymorphic"] is False
        path.write_text(json.dumps({
            "schema_version": SCHEMA_VERSION, "source": {},
            "loci": {"HD_HTT": {"polymorphic": True, "leads": []}},
        }))
        assert _locus(client)["expression_effects"]["polymorphic"] is True


def test_an_unreadable_file_is_reported_without_breaking_the_locus_page(tmp_path, monkeypatch):
    path = tmp_path / "ann.json"
    path.write_text('{"schema_version": 99, "loci": {}, "source": {}}')
    monkeypatch.setattr(main, "SC_ETR_ANNOTATIONS_JSON", str(path))
    with TestClient(app) as client:
        locus = _locus(client)
    assert locus["expression_effects"]["available"] is False
    assert "schema_version" in locus["expression_effects"]["error"]
    assert locus["gene"] == "HTT"
