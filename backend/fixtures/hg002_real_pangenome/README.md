# HG002 TRGT calls for the Pangenome Repeat Catalog

`HG002.GRCh38.trgt.pangenome.vcf.gz` holds real, unmodified TRGT v5.1.0 calls for HG002/NA24385 at
the 5,709 loci of `project/uploads/pangenome_repeat_catalog.hg38.bed`, extracted from the same
genome-wide run as `../hg002_real_genomewide/` (real `bcftools view -R` extraction, then a TRID
filter, by `backend/scripts/real_data/build_pangenome_catalog.py`; see `project/uploads/README.md`).
Every locus has a record; TRGT called all but one.

It is its own directory because it is a second VCF of the same sample: `catalog_vcf_dirs` in
`config/modes.json` names it for the `pangenome` catalog, and keeping it apart from
`../hg002_real_genomewide/` keeps each directory to one VCF per sample.

A mode with no `catalog_vcf_dirs` entry for the catalog falls back to its `catalog_vcf_dir`
genome-wide VCF, and a locus is genotyped there only if a TRGT record carries its exact ID.
