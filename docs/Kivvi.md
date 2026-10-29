# Kivvi suite

The Kivvi suite shows the output of PacBio's [Kivvi](https://github.com/PacificBiosciences/kivvi), which genotypes the D4Z4 repeat (FSHD) and the LPA KIV2 repeat from HiFi whole-genome data. It only reads Kivvi's files; it does not run Kivvi.

Backend: `backend/kivvi/` (routes under `/api/kivvi`). Frontend: `web/src/kivvi/`. Tests: `backend/tests/test_kivvi.py`.

## Data

Samples are found under `$KIVVI_RESULTS_ROOT` (default `~/kivvi_tool/runs`), laid out as `<root>/<sample>/<d4z4|kiv2>/<sample>.kivvi.<target>.{json,vcf,svg}` (the VCF can also be `.vcf.gz`, read as is; if both exist the plain `.vcf` wins). Adding a sample is a matter of putting a Kivvi run there. The root lives outside the repository, and real per-participant outputs are never committed.

## What is and is not applied

Kivvi defines one clinical rule, which is applied: contraction of D4Z4 to 1–10 units on a chr4 qA haplotype with an intact PolyA signal is the FSHD1 range, except on the 4A166 (`Group2.1`) haplotype, and a contracted qA allele on another chromosome is flagged for review. Kivvi defines no methylation cut-off and no KIV2 short-allele cut-off, so none is applied.

## FSHD2 evidence

The **FSHD2 evidence** screen lays out what Kivvi reports that bears on FSHD2. It makes no call. The Kivvi preprint (Chen et al., bioRxiv 2026, doi 10.64898/2026.04.10.717730) gives no FSHD2 rule or cut-off, and validated Kivvi on FSHD1 and control samples, not on confirmed FSHD2 cases.

The screen shows:

- **Sample methylation**, Kivvi's `methylation_all_sites`: methylated reads over all reads across every D4Z4 unit. Kivvi's documentation notes it is low in FSHD2 samples.
- **Each allele**: size, chromosome, distal haplotype and methylation (Kivvi's per-allele value, which covers the allele's last three repeat units), and whether it is *permissive*: chr4, `qAIntactPolyA`, and not the 4A166 haplotype, the kind that can express DUX4.
- **Where each size sits**: 1–10 units (the FSHD1 range), 11–20 (the range the preprint calls commonly indicated for FSHD2), above 20, or not resolved for a partly assembled allele whose lower bound does not settle it. These labels describe a size; nothing depends on them.
- **What a finding rests on**: a permissive allele (shown), hypomethylation across D4Z4 including non-permissive alleles (shown, no cut-off), and a pathogenic variant in SMCHD1, DNMT3B or LRIF1, which is **not assessed**: Kivvi does not report those genes, and the screen does not look for them.

Cautions from the preprint are on the screen: the six candidate SMCHD1 samples had unknown FSHD status; Kivvi's accuracy in the 11–20 range was not measured; about 10% of FSHD2 cases carry a permissive allele of 8–10 units; a contracted permissive allele with low methylation can mean FSHD1 and FSHD2 together; the preprint left cell lines out of its methylation analysis because methylation in lymphoblastoid cell lines can be unreliable.

The preprint summarises an allele's methylation as the median over its last five repeat units, and the sample as the median over all CpG sites in all units; the installed Kivvi (v1.1.0) documents the per-allele value over the last three units, and that is what the screen shows.

## Limits

- The only sample on hand is PacBio's public demo (HG03453), a control with all alleles above 10 units. The FSHD1-range and review paths, and everything on the FSHD2 screen other than its no-contraction, high-methylation case, are covered by synthetic tests only, not by real data.
- A sample with Kivvi output but no small-variant calls has no way to show SMCHD1, DNMT3B or LRIF1 variants here.
