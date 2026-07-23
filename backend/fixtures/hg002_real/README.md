# HG002 real-data fixture

`HG002.strchive_loci.trgt.vcf.gz` is a genuine TRGT v5.1.0 genotype call on real HG002
PacBio HiFi reads at all 81 STRchive disease loci (80 original + FRAXG_BCLAF3, added
2026-07-23) -- not simulated, not hand-built.

## Provenance

- **Reads**: `HG002_GRCh38.haplotag.10x.bam`, PacBio SequelII CCS (HiFi) reads for
  HG002/NA24385 (Ashkenazim son), aligned to GRCh38, from NIST/GIAB's public dataset
  (`s3://giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_SequelII_CCS_11kb/HG002_GRCh38/`).
  Only reads overlapping the 80 STRchive loci (+/-2kb padding) were fetched, via
  `samtools view` remote range requests against the BAM's `.bai` index -- the full BAM
  is 69GB and was never downloaded in full.
- **Reference**: `GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz`, the same
  reference GIAB aligned these reads to (from the same S3 bucket, under
  `release/references/GRCh38/`). Per-chromosome sequences were rebuilt at their real,
  full lengths (so trgt's contig-length checks against the BAM header pass) but only
  actually populated with real bases in a window around each locus -- the rest is
  N-filled placeholder, since trgt never reads there.
- **Repeat catalog**: `strchive_trgt.bed`, built by combining STRchive's 80 disease
  loci (coordinates, IDs) with PacBio's own official
  `PacificBiosciences/trgt` repo's `repeats/pathogenic_repeats.hg38.bed` catalog for
  MOTIFS/STRUC where the two matched by gene + coordinates (53/80), falling back to a
  simple single-motif STRUC built from STRchive's own motif field otherwise.
- **Genotyper**: the actual `trgt` v5.1.0 Linux binary, downloaded from PacBio's own
  GitHub release (`PacificBiosciences/trgt` release `v5.1.0`), run as:

  ```
  trgt genotype --genome padded_reference.fasta --reads hg002_strchive_loci.bam \
      --repeats strchive_trgt.bed --output-prefix hg002_real --sample-name HG002 \
      --karyotype XX --threads 4
  ```

  Note the `--karyotype XX`: HG002/NA24385 is real male donor, but this run was
  genotyped with a simplifying diploid-XX assumption on chrX rather than proper
  hemizygous calling (unlike the real male-hemizygous carrier-panel run in
  `fixtures/carrier_panel`, which uses `trgt`'s actual chrX/chrY hemizygous calling
  -- see that fixture's own README). `sample_sex.json`'s "M" for HG002 is the real
  donor's sex, sourced from GIAB's own sample metadata, independent of how this
  particular VCF happened to be genotyped.

## Why this matters

This is real HG002 HiFi data run through the real PacBio genotyper -- the actual
research-use benchmark sample NIST/GIAB and the broader genomics community use as a
truth set. It's licensed "For Research Use Only. Not for use in diagnostic
procedures." per PacBio's GIAB data terms, which is exactly this project's use.

Running the real pipeline against real data (rather than only the hand-built
`synthetic_trgt.vcf` fixture) surfaced two genuine correctness bugs in the annotation
code that the synthetic fixture never exercised:

1. Compound-motif loci (e.g. HTT's `(CAG)nCAACAG(CCG)n`) report one MC sub-count per
   STRUC segment -- summing them (instead of picking out the disease-relevant motif's
   own segment) would misclassify benign CAG counts as intermediate/pathogenic due to
   an unrelated flanking CCG tract.
2. TRGT's FORMAT arrays (AL/ALLR/SD/MC/AP/AM) are indexed by rank among the two
   *called* haplotypes' original allele indices, not sized to cover every declared
   REF+ALT allele -- a GT like `1/2` (neither haplotype calling REF) gets a 2-element
   array for indices {1,2}, not a 3-element REF+ALT1+ALT2 array.

Both are fixed in `repeatscope/trgt.py` and covered by regression tests in
`tests/test_trgt.py` using these exact real values.

HG002 is a well-characterized reference sample with no known pathogenic repeat
expansions, so a mostly-Normal panel with one real intermediate finding (ATXN2, a
known real population polymorphism) is the expected, correct result -- see
`HMNR7_VWA1`'s note in `strchive.py`'s `_is_contraction_locus` docstring for one
remaining known STRchive data-modeling edge case this pipeline can't fully resolve.

## FRAXG_BCLAF3 (added 2026-07-23)

Sourced from a *different*, genuinely genome-wide TRGT run for this same donor
(`HG002_rep1_analysis_v3_0_2/HG002.GRCh38.trgt.sorted.phased.vcf.gz`, a full
PacBio HiFi-Solve-style workflow output, not the padded/loci-only BAM described
above) -- the record at chrX:19,990,892-19,991,003 (`TRID=chrX_19990892_19991003_GCC`)
was extracted as-is (unmodified REF/ALT/FORMAT) and appended to this file.

That TRID and region are wider than FRAXG_BCLAF3's own catalog entry
(chrX:19,990,922-19,990,973, motif CCG) because the genome-wide run's catalog defines
this repeat's boundaries more broadly and reads the same repeat unit in a different
phase (GCC is a cyclic rotation of CCG -- confirmed via `motifs.is_cyclic_rotation`).
The record's real TRID/coordinates were kept unmodified rather than relabeled to match
this fixture's other rows' disease-specific IDs, preserving traceability to its actual
source; `annotate_strchive_panel` resolves it to `FRAXG_BCLAF3` via
`TrgtIndex.find_overlap`'s coordinate-overlap fallback, verified end-to-end.

The record is a genuine hemizygous TRGT call (`GT=1`, single-value FORMAT arrays --
HG002 is male, and this genome-wide run used proper chrX/chrY hemizygous calling,
unlike this fixture's other 80 loci, genotyped with `--karyotype XX`; see above).
`repeatscope/trgt.py`'s `_parse_gt`/`TrgtCall.hemizygous` already handle this VCF
convention correctly, and the record is otherwise untouched -- no fields fabricated
or adjusted. Real value: 21 GCC repeats (depth 15), within FRAXG_BCLAF3's benign
range (<=23) -- Normal, consistent with HG002 having no known pathogenic repeat
expansions.
