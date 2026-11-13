# Kivvi suite

The Kivvi suite shows the output of PacBio's [Kivvi](https://github.com/PacificBiosciences/kivvi), which genotypes the D4Z4 repeat (FSHD) and the LPA KIV2 repeat from HiFi whole-genome data. It only reads Kivvi's files; it does not run Kivvi.

Backend: `backend/kivvi/` (routes under `/api/kivvi`). Frontend: `web/src/kivvi/`. Tests: `backend/tests/test_kivvi.py`.

## Data

Samples are found under `$KIVVI_RESULTS_ROOT` (default `~/kivvi_tool/runs`), laid out as `<root>/<sample>/<d4z4|kiv2>/<sample>.kivvi.<target>.{json,vcf,svg}` (the VCF can also be `.vcf.gz`, read as is; if both exist the plain `.vcf` wins). Adding a sample is a matter of putting a Kivvi run there. The root lives outside the repository, and real per-participant outputs are never committed.

## What is and is not applied

Kivvi defines one clinical rule, which is applied: contraction of D4Z4 to 1–10 units on a chr4 qA haplotype with an intact PolyA signal is the FSHD1 range, except on the 4A166 (`Group2.1`) haplotype, and a contracted qA allele on another chromosome is flagged for review. Kivvi defines no methylation cut-off and no KIV2 short-allele cut-off, so no KIV2 cut-off is applied. The one cut-off in the app is a hypomethylation cut-off on the FSHD2 evidence screen (below), set by the project rather than taken from Kivvi.

## FSHD2 evidence

The **FSHD2 evidence** screen lays out what Kivvi reports that bears on FSHD2. It makes no FSHD2 call. The Kivvi preprint (Chen et al., bioRxiv 2026, doi 10.64898/2026.04.10.717730) gives no FSHD2 rule or methylation cut-off, and validated Kivvi on FSHD1 and control samples, not on confirmed FSHD2 cases.

The screen shows:

- **Sample methylation**, Kivvi's `methylation_all_sites`: methylated reads over all reads across every D4Z4 unit. Kivvi's documentation notes it is low in FSHD2 samples.
- **Hypomethylation**: a sample is **positive** when its sample methylation is below **30%** (negative at 30% or above, no data when Kivvi reports none). The 30% cut-off was set for this tool by the project team; it is not from Kivvi or the preprint, and it is applied to the sample-level value, not to individual alleles.
- **Each allele**: size, chromosome, distal haplotype and methylation (Kivvi's per-allele value, which covers the allele's last three repeat units), and whether it is *permissive*: chr4, `qAIntactPolyA`, and not the 4A166 haplotype, the kind that can express DUX4.
- **Where each size sits**: 1–10 units (the FSHD1 range), 11–20 (the range the preprint calls commonly indicated for FSHD2), above 20, or not resolved for a partly assembled allele whose lower bound does not settle it. These labels describe a size; nothing depends on them.
- **What a finding rests on**: a permissive allele, hypomethylation across D4Z4 (scored against the 30% cut-off), and a variant in SMCHD1, DNMT3B or LRIF1 (supplied by you, see below). The evidence tier counts how many of the three a sample shows.

### Evidence tiers

The tiers are a research view of how much of an FSHD2 picture a sample shows, not a diagnosis. FSHD2 needs a DUX4-permissive allele plus a second hit that lowers D4Z4 methylation, so the permissive allele is a gate:

| Tier | Meaning |
| --- | --- |
| **1 · all three findings** | permissive allele, hypomethylation (below 30%) and a counting variant |
| **2 · permissive allele and one other finding** | permissive allele plus either hypomethylation or a counting variant |
| **3 · permissive allele only** | a permissive allele with neither |
| **Not supported · no permissive allele found** | no tier, however low the methylation or damaging the variant |

A variant **counts** when it is in SMCHD1, DNMT3B or LRIF1 and is either loss-of-function (frameshift, nonsense or canonical ±1/2 splice, read from its HGVS notation) or classified pathogenic / likely pathogenic in ClinVar. A missense variant of uncertain significance is shown but does not count. A sample missing from the variants file is **not assessed**, which is not the same as none found, so its tier may still rise (the screen says so). "No permissive allele found" means none among the alleles Kivvi assembled: with fewer than four resolved, or an allele Kivvi could not place on a chromosome, a permissive allele may be missing from the data, and the screen notes it.

### Supplying the variants

Kivvi reports D4Z4 only, so the variants come from you, as a JSON file at `$KIVVI_FSHD2_VARIANTS` or, by default, `fshd2_variants.json` in the results root (outside the repo; it holds participant data and is never committed):

```json
{
  "samples": {
    "SAMPLE_A": {"variants": [
      {"gene": "SMCHD1", "hgvs_c": "c.181dup", "hgvs_p": "p.(C61Lfs*9)",
       "chrom": "chr18", "pos": 2656255, "build": "hg38", "genotype": "heterozygous",
       "clinvar_url": "https://www.ncbi.nlm.nih.gov/clinvar/variation/<id>/",
       "clinvar_classification": "Pathogenic", "note": "optional free text"}
    ]},
    "SAMPLE_B": {"variants": []}
  }
}
```

`gene`, `hgvs_c` and (when there is one) `hgvs_p` are what the tiers read; the rest is shown. An empty list means the sample was checked and none was found; leave a sample out if it was not checked. The file is reread when it changes.

Cautions from the preprint are on the screen: the six candidate SMCHD1 samples had unknown FSHD status; Kivvi's accuracy in the 11–20 range was not measured; about 10% of FSHD2 cases carry a permissive allele of 8–10 units; a contracted permissive allele with low methylation can mean FSHD1 and FSHD2 together; the preprint left cell lines out of its methylation analysis because methylation in lymphoblastoid cell lines can be unreliable.

The preprint summarises an allele's methylation as the median over its last five repeat units, and the sample as the median over all CpG sites in all units; the installed Kivvi (v1.1.0) documents the per-allele value over the last three units, and that is what the screen shows.

## Limits

- The only sample on hand is PacBio's public demo (HG03453), a control with all alleles above 10 units. The FSHD1-range and review paths, and everything on the FSHD2 screen other than its no-contraction, high-methylation case, are covered by synthetic tests only, not by real data.
- A sample with Kivvi output but no small-variant calls has no way to show SMCHD1, DNMT3B or LRIF1 variants here.
