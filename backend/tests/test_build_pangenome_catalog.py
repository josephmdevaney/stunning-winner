import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "real_data" / "build_pangenome_catalog.py"
spec = importlib.util.spec_from_file_location("build_pangenome_catalog", SCRIPT)
build = importlib.util.module_from_spec(spec)
spec.loader.exec_module(build)

PORTAL_ID = build.PORTAL_ID

CATALOG = "\n".join(
    [
        "chr1\t100\t120\tID=1-100-120-CAG;MOTIFS=CAG;STRUC=<TR:1-100-120-CAG>",
        # a variation cluster: its ID is its constituents' IDs, comma-joined
        "chr1\t200\t260\tID=1-200-230-GA,1-235-260-GA;MOTIFS=GA;STRUC=<VC:1-200-260>",
        "chr1\t400\t430\tID=1-400-430-T;MOTIFS=T;STRUC=<TR:1-400-430-T>",
        "chr2\t500\t540\tID=2-500-540-CCG;MOTIFS=CCG;STRUC=<TR:2-500-540-CCG>",
    ]
) + "\n"


def paper_row(table, locus, chrom, start, end, portal_id):
    return {
        "Source table": table, "Paper locus_id": locus, "Chrom": chrom,
        "Start (0-based)": str(start), "End": str(end), PORTAL_ID: portal_id,
    }


def resolve(tmp_path, rows):
    bed = tmp_path / "catalog.bed"
    bed.write_text(CATALOG)
    return build.resolve(rows, str(bed))


def test_a_portal_id_selects_the_locus_it_names(tmp_path):
    by_id, left_out = resolve(tmp_path, [paper_row("3a", "a", "chr1", 100, 120, "1-100-120-CAG")])
    assert list(by_id) == ["1-100-120-CAG"] and left_out == []


def test_a_portal_id_inside_a_cluster_selects_the_cluster_once(tmp_path):
    rows = [
        paper_row("3a", "a", "chr1", 200, 230, "1-200-230-GA"),
        paper_row("3e", "b", "chr1", 235, 260, "1-235-260-GA"),
    ]
    by_id, left_out = resolve(tmp_path, rows)
    assert list(by_id) == ["1-200-230-GA,1-235-260-GA"] and left_out == []


def test_a_portal_id_the_catalog_lacks_falls_back_to_the_locus_overlapping_the_paper_locus_most(tmp_path):
    # 250-420 overlaps the cluster by 10bp and the 1-400-430-T locus by 20bp
    by_id, left_out = resolve(tmp_path, [paper_row("3f", "c", "chr1", 250, 420, "1-999-999-A")])
    assert list(by_id) == ["1-400-430-T"] and left_out == []
    # no portal match at all: 190-265 overlaps the cluster by 60bp and nothing else
    by_id, _ = resolve(tmp_path, [paper_row("3f", "c", "chr1", 190, 265, "")])
    assert list(by_id) == ["1-200-230-GA,1-235-260-GA"]


def test_a_paper_locus_nothing_overlaps_is_left_out_not_invented(tmp_path):
    row = paper_row("3e", "d", "chr2", 900, 950, "2-900-950-A")
    by_id, left_out = resolve(tmp_path, [row])
    assert by_id == {} and left_out == [row]


def test_rows_without_coordinates_are_split_off(tmp_path):
    csv_path = tmp_path / "found.csv"
    csv_path.write_text(
        "Source table,Paper locus_id,Chrom,Start (0-based),End,TRExplorer portal match: locus ID\n"
        "3a,MUC6,,,,\n"
        "3a,a,chr1,100,120,1-100-120-CAG\n"
    )
    with_coords, without = build.read_found_loci(str(csv_path))
    assert [r["Paper locus_id"] for r in with_coords] == ["a"]
    assert [r["Paper locus_id"] for r in without] == ["MUC6"]
