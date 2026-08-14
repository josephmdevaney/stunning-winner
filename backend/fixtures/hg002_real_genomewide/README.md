# HG002 genome-wide TRGT fixture (Git LFS)

`HG002.GRCh38.trgt.genomewide.vcf.gz` is the real, unmodified output of a genuinely
genome-wide TRGT run for HG002/NA24385, genotyped against the Broad Institute's
TRExplorer v2.0 catalog -- see `project/uploads/README.md` for the catalog's own
provenance.

## Provenance (replaced 2026-08-14)

- **Genotyper**: the real `trgt` v5.1.0 binary (same version already used
  elsewhere in this project -- see `../hg002_real/README.md`), run against real
  HG002 PacBio HiFi reads by the user directly (own local TRGT build/run, not
  downloaded from a public bucket like the loci-only fixture).
- **Catalog**: TRExplorer v2.0's `variation_clusters_and_isolated_TRs`
  catalog -- 5,339,786 regions, hg38, the same file this project's own
  `project/uploads/repeat_catalog.hg38.bed` now uses. TRID-for-TRID identical
  to the catalog: every single one of the catalog's 5,339,786 regions has a
  real matching genotype call in this VCF (verified directly --
  `TrgtIndex.by_trid` hits 100% of catalog region IDs, no coordinate-overlap
  fallback needed at all for this pairing).
- **Copied unmodified**: same TRID/coordinates/REF/ALT/FORMAT values as the
  source file, no records added, removed, or edited.
- **Tracked via Git LFS** (see `.gitattributes`): at 278MB gzip-compressed,
  this exceeds GitHub's 100MB single-file limit for a regular git blob.

## What it's used for

- This file is the **sample genotype source** for the genome-wide "Repeat Catalog"
  screen: `catalog_vcf_dir` in `config/modes.json` points genome-wide mode's catalog
  lookups here specifically, separately from `vcf_dir` (which still points at the
  smaller, 81-locus `../hg002_real/HG002.strchive_loci.trgt.vcf.gz` used by the
  Overview/Locus Detail disease panel). Keeping them separate matters:
  auto-discovering both VCFs from the same directory would silently let one
  overwrite the other as "the" path for sample HG002, since both files report the
  same sample name.
- **Performance note**: parsing this VCF's ~5.34M records into `TrgtIndex` takes
  ~2 minutes on first access per (mode, sample) -- a one-time cost, cached for
  the process's lifetime afterward (see `repeatscope/api/main.py`'s
  `_get_trgt_index`). Not yet optimized further; `iter_trgt_records`'s
  per-record FORMAT-field parsing is the dominant cost, same bottleneck
  identified (but not addressed) when benchmarking the old, smaller VCF.

## History

- **2026-07-23 (first added)**: a 937,195-region genome-wide run, itself the
  original source of the (now-replaced) 937,195-region repeat catalog.
- **2026-08-14 (replaced)**: swapped for this 5,339,786-region run against
  TRExplorer v2.0, alongside the matching catalog swap in
  `project/uploads/repeat_catalog.hg38.bed`.
