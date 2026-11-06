"""Kivvi suite: viewing PacBio Kivvi (HiFi large-unit VNTR genotyper) results for
LPA KIV2 and D4Z4 inside the stunning-winner app.

This package only *reads* Kivvi's outputs (`<prefix>.kivvi.<target>.json/.vcf/.svg`).
Kivvi itself is PacBio-proprietary and is neither vendored nor redistributed here;
run it separately against a WGS BAM and point KIVVI_RESULTS_ROOT at the outputs.
"""
