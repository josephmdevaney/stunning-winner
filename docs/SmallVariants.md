# Small variants (SNV / indel)

The **Small variants** suite is a genome-wide review of a PacBio HiFi small-variant call set (DeepVariant output, or anything shaped like it). It does not start from a gene list: every call in the VCF is in the table, and the filters narrow it. Each call is annotated with its gene body, ClinVar's classification, and population frequency from gnomAD and HPRC.

It shows what the evidence says and where there is none. It does not classify variants: the ClinVar assertion is ClinVar's, and an ACMG/AMP criterion with no data behind it is shown as *not evaluated*, never as met or ruled out.

Backend: `backend/smallvar/` (routes under `/api/smallvar`). Frontend: `web/src/smallvar/`. Tests: `backend/tests/test_smallvar.py`.

## Data

Samples are single-sample VCFs under `$SMALLVAR_RESULTS_ROOT` (default `~/smallvar_data/samples`), laid out as `<root>/<sample>/<name>.vcf[.gz]`, with an optional `sample.json`:

```json
{ "label": "…", "sex": "male", "sex_source": "…", "caller": "…" }
```

A small public HG002 slice (the GIAB v4.2.1 truth set, chr22) ships in `backend/fixtures/hg002_smallvar/` so the suite works on a fresh checkout. Real samples stay outside the repo.

Each VCF is read once into a SQLite index under `$SMALLVAR_CACHE` (default `~/smallvar_data/cache`) on a background thread, and the screens poll `/status` until it is ready. The index is rebuilt when the VCF or any annotation resource changes. On the whole HG002 genome (4,095,610 calls) the build takes about 105 seconds; filtered queries then return in 0.02-0.5 seconds. A multi-allelic record becomes one call per ALT allele the genotype carries.

## Annotation resources

All optional. A missing resource leaves its columns empty, the API says so, and a filter that depends on it is left out rather than applied (otherwise "ClinVar pathogenic" would match nothing and read as a negative result).

| Resource | Where it comes from | Notes |
|---|---|---|
| ClinVar | `scripts/real_data/build_smallvar_clinvar.py` from NCBI `variant_summary.txt.gz` → `~/smallvar_data/clinvar_plp.grch38.tsv.gz` (`SMALLVAR_CLINVAR`) | GRCh38 small variants classified pathogenic, likely pathogenic or conflicting (October 2026: 529,569 variants). A variant outside this subset is not shown to be absent from ClinVar. Matched by exact allele on trimmed alleles. |
| gnomAD v3.0 genomes | `gnomad.hg38.v3.custom.v1.zip` in PacBio's HiFi WDL resource bundle v1.0.2 (`static_resources/GRCh38/slivar_gnotate/`), `SMALLVAR_GNOMAD` | PASS sites only; overall allele frequency, allele count and homozygote count (no allele number, no ancestry-specific or popmax frequency). Built 2020-10-08 from gnomAD release 3.0 (about 71,700 genomes). Absent means not seen at a PASS site in that release. |
| HPRC | `hprc.deepvariant.glnexus.hg38.v1.zip` in the same bundle, `SMALLVAR_HPRC` | About 100 long-read samples. Useful for spotting recurrent long-read calls, not for ruling variants out. The cohort includes GIAB HG001, HG003 and HG004, so for related samples such as HG002 it is not independent. |
| Gene bodies | `project/uploads/genes.hg38.bed` (`SMALLVAR_GENES_BED`) | RefSeq gene loci. Overlap with a gene body only, not coding consequence. |

### Reading the PacBio frequency files

The two frequency resources are slivar "gnotate" zips, which only slivar documents how to read. `backend/smallvar/gnotate.py` reads them directly in Python, without slivar. The layout was worked out from the data and checked against HG002: 98.4% of chr22 SNVs, 98.2% of deletions and 97.3% of insertions are found by their predicted code, and a control with shuffled alleles finds none.

Per chromosome, a sorted array of uint64 codes plus parallel float32 arrays (the allele count and homozygote count are float32 too):

- bits 36-63: position, **0-based** (VCF POS - 1)
- bits 32-35: `len(ref)`; bits 28-31: `len(alt)`
- bits 2-27: the ref and alt bases, 2 bits each (A, C, G, T = 0..3), first base most significant; bits 0-1 reserved
- alleles with `len(ref) + len(alt) > 13` do not fit and are listed, with their values, in `long-alleles.txt` (0-based position again)

Matching uses the alleles as written in the VCF, as gnomAD writes them. The arrays are deflate-compressed, so a chromosome is decompressed whole; chromosomes are processed one at a time.

## What the screens show

- **Overview**: call-set metrics counted from the file; the funnel from all calls to the review list (each step is a real filter on the Variants screen, so the count can be reproduced there); ClinVar matches and the genes carrying them; and which ACMG/AMP criteria have evidence here and which do not.
- **Variants**: every call, filtered by ClinVar class and review stars, gnomAD and HPRC frequency (with a choice for calls absent from gnomAD), zygosity, type, chromosome, a gene body, a gene list (case-insensitive), and GQ; sortable and paged. A filter combination that returns nothing shows how many calls survive each filter in turn. A row opens to the ClinVar record, frequencies, genotype and quality, and flags a ClinVar P/LP call whose gnomAD frequency is above the 5% stand-alone benign (BA1) threshold. The ACMG SF v3.3 gene list is available as a one-click gene-list fill; the suite does not apply the ACMG SF reporting rules.

Zygosity on chrX outside the pseudo-autosomal regions is shown as hemizygous for a male sample. Sex comes from `sample.json` when given, with its source; otherwise it is inferred from the heterozygous fraction of chrX non-PAR calls against the autosomes, and left *undetermined* when there are too few chrX calls.

## Reference data

`backend/smallvar/reference/acmg_sf_v3_3.tsv` is the ACMG SF v3.3 gene table (84 genes, with the inheritance, reportable genotype, variant restriction and rule IDs from the implementation spec). A test checks it against the counts the spec states. It feeds the gene-list fill today and is ready for a secondary-findings screen (rules R1-R9).

## Limits

- No coding consequence, in-silico score (CADD, REVEL, SpliceAI, AlphaMissense), phenotype (HPO) or parental genotype, so PVS1/PM4/BP3, PP3/BP4/BP7, PP4 and PM3/PS2/PM6 are not evaluated. The PacBio bundle also carries Ensembl gene models and the reference FASTA, which would allow consequence annotation.
- gnomAD is v3.0, not the current release. Replacing it means building a newer gnotate file; the reader needs no change.
- The only sample tested is the HG002 GIAB truth set, which is a consensus of several callers, not DeepVariant output: it is unphased, and its GQ and DP are summed across datasets, so they are not caller-style values. Phasing, GQ and VAF handling are covered by synthetic tests, not yet by a real DeepVariant callset.
- Compound heterozygosity is not called: the VCF's phase sets are stored per call, but no screen uses them yet.
