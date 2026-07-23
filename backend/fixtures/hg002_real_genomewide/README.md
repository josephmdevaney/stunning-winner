# HG002 genome-wide TRGT fixture (Git LFS)

`HG002.GRCh38.trgt.genomewide.vcf.gz` is the real, unmodified output of a genuinely
genome-wide TRGT run for HG002/NA24385 -- 937,195 repeat regions with real genotype
calls, not the padded/loci-only run described in `../hg002_real/README.md`.

## Provenance

- **Source**: `HG002_rep1_analysis_v3_0_2/HG002.GRCh38.trgt.sorted.phased.vcf.gz`, one
  output file of a full PacBio HiFi-Solve-style workflow run for this donor (the same
  directory also has the BAM, structural variants, small variants, methylation pileups,
  etc. for the same run -- this fixture only takes the TRGT VCF).
- **Copied unmodified**: same TRID/coordinates/REF/ALT/FORMAT values as the source file,
  no records added, removed, or edited.
- **Tracked via Git LFS** (see `.gitattributes`): at 119MB gzip-compressed, this exceeds
  GitHub's 100MB single-file limit for a regular git blob.

## What it's used for

- `project/uploads/repeat_catalog.hg38.bed` (the genome-wide "Repeat Catalog" screen's
  region definitions) was built by extracting this same VCF's own TRID/END/MOTIFS/STRUC
  INFO fields -- see the catalog-replacement commit.
- This file is the actual **sample genotype source** for that screen: `catalog_vcf_dir`
  in `config/modes.json` points genome-wide mode's catalog lookups here specifically,
  separately from `vcf_dir` (which still points at the smaller, 81-locus
  `../hg002_real/HG002.strchive_loci.trgt.vcf.gz` used by the Overview/Locus Detail
  disease panel). Keeping them separate matters: auto-discovering both VCFs from the
  same directory would silently let one overwrite the other as "the" path for sample
  HG002, since both files report the same sample name.
