"""Tests for coverage-based sex inference from mosdepth summary files."""

from roh import sex_inference


def _write_summary(tmp_path, rows):
    path = tmp_path / "sample.mosdepth.summary.txt"
    lines = ["chrom\tlength\tbases\tmean\tmin\tmax"]
    for chrom, length, mean in rows:
        bases = int(length * mean)
        lines.append(f"{chrom}\t{length}\t{bases}\t{mean}\t0\t100")
        lines.append(f"{chrom}_region\t{length}\t{bases}\t{mean}\t0\t100")
    path.write_text("\n".join(lines) + "\n")
    return path


def _autosome_rows(mean=16.0):
    return [(f"chr{i}", 100_000_000, mean) for i in range(1, 23)]


def test_infers_male_from_half_chrx_and_present_chry(tmp_path):
    rows = _autosome_rows() + [("chrX", 100_000_000, 8.0), ("chrY", 50_000_000, 5.0)]
    path = _write_summary(tmp_path, rows)
    result = sex_inference.infer_sex_from_mosdepth(path)
    assert result["sex"] == "M"
    assert result["chrX_autosomal_ratio"] == 0.5
    assert result["chrY_autosomal_ratio"] > sex_inference.CHRY_PRESENT_MIN_RATIO


def test_infers_female_from_full_chrx_and_absent_chry(tmp_path):
    rows = _autosome_rows() + [("chrX", 100_000_000, 16.0), ("chrY", 50_000_000, 0.05)]
    path = _write_summary(tmp_path, rows)
    result = sex_inference.infer_sex_from_mosdepth(path)
    assert result["sex"] == "F"
    assert result["chrX_autosomal_ratio"] == 1.0


def test_ambiguous_coverage_does_not_guess(tmp_path):
    # chrX ratio between the two clusters, chrY present but weak -- not a clean
    # signal either way, so this must not be forced into M or F.
    rows = _autosome_rows() + [("chrX", 100_000_000, 12.0), ("chrY", 50_000_000, 1.0)]
    path = _write_summary(tmp_path, rows)
    result = sex_inference.infer_sex_from_mosdepth(path)
    assert result["sex"] is None


def test_missing_chrx_returns_none(tmp_path):
    path = _write_summary(tmp_path, _autosome_rows())
    result = sex_inference.infer_sex_from_mosdepth(path)
    assert result["sex"] is None
    assert "reason" in result
