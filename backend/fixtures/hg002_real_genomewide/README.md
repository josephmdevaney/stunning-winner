# HG002 genome-wide TRGT fixture

`HG002.GRCh38.trgt.genomewide.vcf.gz` holds real, unmodified genotype calls
from a genuinely genome-wide TRGT run for HG002/NA24385, currently rescoped
to just the 24,956 regions RExPRT predicts pathogenic -- see
`project/uploads/README.md` for the catalog's own provenance and why.

## Provenance (rescoped 2026-08-21)

- **Rescoped from**: the real, unmodified full genome-wide TRGT VCF described
  under "History" below, via real `bcftools view -R` region extraction
  against the 24,956 RExPRT-positive regions (see
  `backend/scripts/real_data/build_rexprt_pathogenic_subset.py`) -- not
  resimulated, re-genotyped, or approximated. Same TRID/coordinates/REF/ALT/
  FORMAT values as the source records for every region kept.
- **Coverage**: every one of the catalog's 24,956 regions has a real matching
  genotype call in this VCF (TRID-exact match verified directly --
  `TrgtIndex.by_trid` hits 100% of catalog region IDs).
- **No longer needs Git LFS**: at ~1.5MB gzip-compressed (down from 278MB
  before rescoping), this is well under GitHub's 100MB single-file limit for
  a regular git blob -- see `.gitattributes`.

## What it's used for

- This file is the **sample genotype source** for the genome-wide "Repeat Catalog"
  screen: `catalog_vcf_dir` in `config/modes.json` points genome-wide mode's catalog
  lookups here specifically, separately from `vcf_dir` (which still points at the
  smaller, 81-locus `../hg002_real/HG002.strchive_loci.trgt.vcf.gz` used by the
  Overview/Locus Detail disease panel). Keeping them separate matters:
  auto-discovering both VCFs from the same directory would silently let one
  overwrite the other as "the" path for sample HG002, since both files report the
  same sample name.

## History

- **2026-07-23 (first added)**: a 937,195-region genome-wide run, itself the
  original source of the repeat catalog at the time.
- **2026-08-14 (replaced)**: swapped for a full 5,339,786-region run genotyped
  directly against the full TRExplorer v2.0 catalog. Parsing that VCF's ~5.34M
  records into `TrgtIndex` measured ~2 minutes / several GB RAM on first
  access per (mode, sample) under normal conditions -- and, on 2026-08-21,
  was observed taking 15+ minutes without finishing on a machine under real
  memory pressure (~60MB free system RAM at the time) -- a real, measured
  scaling cost, not a theoretical one.
- **2026-08-21 (rescoped)**: narrowed to just the 24,956 RExPRT-predicted
  -pathogenic regions in response to the cost above, once it was clear the
  product's actual use of this screen was "browse candidate loci RExPRT
  flagged," not "browse all 5.3M reference definitions." The full
  5,339,786-record 2026-08-14 VCF remains recoverable from git history (it
  was a real Git LFS object) if the full genome-wide run is ever needed
  again -- e.g. to browse regions RExPRT scored benign or never evaluated.
