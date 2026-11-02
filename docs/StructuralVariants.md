# Structural variants

The **Structural variants** suite is a genome-wide review of the SV calls in a PacBio HiFi call set: sawfish, pbsv or sniffles, as a single sample or a joint call across a family. Every call is in the table and the filters narrow it. Each call is annotated with the genes it touches (at exon level, on the MANE Select transcript), population matches, genome context, and for a trio, an inheritance class.

It shows what the evidence says and where there is none. It does not classify variants, and an evidence source that is not installed is shown as *not evaluated*, never as met or ruled out.

Backend: `backend/structvar/` (routes under `/api/structvar`). Frontend: `web/src/structvar/`. Tests: `backend/tests/test_structvar.py`.

## Data

Call sets are SV VCFs under `$STRUCTVAR_RESULTS_ROOT` (default `~/sv_data/samples`), laid out as `<root>/<callset>/<name>.vcf[.gz]`, with an optional `sample.json`:

```json
{ "label": "…", "caller": "…", "support": "reads", "proband": "SAMPLE_A", "parents": ["SAMPLE_B", "SAMPLE_C"] }
```

`support` names what the genotype's `AD` counts, shown beside each call and on the support filter; it is "reads" unless the file says otherwise. The truth set's `AD` is not read depth (a `1|1` is `0,2`, a `1|0` is `1,1`), so its sample.json says "haplotype calls".

The suite ships one public call set, the **GIAB HG002 v5.0q SV benchmark, chr22 only** (`backend/fixtures/hg002_structvar/`: 2,276 deletions and insertions, the GRCh38 `stvar` VCF cut to chr22 with INFO reduced to `SVTYPE` and `SVLEN`). It is a benchmark derived from the HG002 diploid assembly, not a read-based caller, and the Genome in a Bottle consortium makes it public. Everything else is outside the repo: SV calls from real families are participants' data. With no annotation resources installed, the fixture shows the calls and their genotypes and reports the gene, population and track columns as *not evaluated*.

Each VCF is read once into a SQLite index under `$STRUCTVAR_CACHE` (default `~/sv_data/cache`) on a background thread, and the screens poll `/status` until it is ready. A joint call of a trio indexes in a few seconds; the index is rebuilt when the VCF, an annotation resource, or the declared roles change.

### Reading the VCF

`structvar/vcf.py` normalises how each caller writes things: sequence-resolved insertions and deletions (a full inserted sequence, or a long REF with a one-base ALT) versus symbolic `<DEL>`, `<DUP>`, `<DUP:TANDEM>`, `<INV>`, `<CNV>`; breakends with their mate; imprecise depth-based copy-number events whose genotypes carry `CN`; and half-called genotypes such as `./1`, which still carry the variant. Support comes from whichever fields the caller wrote (`AD`, `DV`/`DR`, `DP`, `GQ`, or `CNQ`).

A genotype is read for how many alleles it carries and how many it writes. A haploid call, such as a bare `1` on chrX or chrY outside the PARs in a male (the truth set writes these), has ploidy 1, so its carrier is **hemizygous** and is not called heterozygous. The zygosity filter keeps this distinction: *heterozygous* needs two written alleles, and *homozygous or hemizygous* takes both a `1/1` and a haploid `1`. A half-called genotype such as `1|.` carries the variant; `0|.` is a reference call on one haplotype and says nothing about the other.

A record whose REF or ALT is not valid VCF is skipped and counted, and the overview says so. In one pbsv file a single insertion at a contig end has a NUL byte for its REF. `bcftools` stops at it and silently drops every record after, so a tool that trusts it loses data. A truncated or corrupt gzip is reported as such by the build (it fails; it is not read partway).

## Who is the child

For a call set with exactly three samples, the proband is chosen from the genotypes: each sample is tried as the child of the other two, counting Mendelian violations on PASS calls made in all three (a child carrying a variant neither parent has, a child homozygous when a parent lacks it, and so on). A sample is accepted only if its violation rate is below half of the runner-up's, and a tie is never accepted. The basis is stored and shown ("X is the child: violation rate 2.3% against 18.5% and 17.5%"). A `sample.json` naming the proband and parents overrides this, and with any other number of samples no family structure is assumed.

An SV the proband carries then has an inheritance class: **inherited** if either parent carries it, **de novo candidate** if both parents have a reference call, **unknown** if a parent has no call and neither carries it. "De novo" here means only that neither parent has a call. On real trios most such calls are a parent's missed call (hundreds per genome against an expected rate of well under one), so the screen shows the parents' read depth beside it, the filters for parental depth and supporting reads, and the funnel's de novo step comes after the filters that remove repeats and segmental duplications.

## Annotation resources

All optional. A missing resource leaves its columns empty, the API reports it, and a filter that depends on it is left out of the query instead of being applied.

| Resource | Source | Notes |
|---|---|---|
| Gene models | `scripts/real_data/build_sv_genes.py`: UCSC's MANE track (`hgdownload.soe.ucsc.edu/gbdb/hg38/mane/mane.bb`, from the NCBI / EMBL-EBI MANE project) plus Ensembl 101 gene models and gene-level files from PacBio's HiFi WDL resource bundle (v1.0.2) | One transcript per gene: the **MANE Select** transcript where a gene has one (about 19,400 genes), so exon numbers match clinical reports and carry the RefSeq `NM_` accession; otherwise the Ensembl 101 transcript with the longest coding sequence. The gene span is the Ensembl locus. |
| Population SVs | `scripts/real_data/build_sv_catalogs.py` from the bundle's `sv_pop_vcfs/` | gnomAD-SV via dbVar (short-read; frequency, per ancestry), HPRC pbsv (long-read; frequency, about 100 samples), DECODE ONT (presence only), EEE (number of discovery genomes). |
| Genome tracks | PacBio's bundle, `GRCh38/annotation/` (UCSC tables) | Segmental duplications, RepeatMasker and tandem repeats (SVs up to 100 kb, queried with `tabix`), and gap / telomere / heterochromatin / PAR / V(D)J regions. These carry UCSC's data-use terms, which are worth checking before commercial use. |
| Gene context | the bundle's HPO, inheritance, gnomAD v2.1.1 LoF and ClinVar gene-condition files | pLI and observed/expected LoF, inheritance, HPO phenotype terms, and ClinVar condition names (as of 2022-12), shown beside a gene hit. |

The MANE file is a bigBed; `structvar/bigbed.py` reads it directly (header, chromosome tree, R-tree of zlib blocks), with no external tool.

### Gene impact

An SV is compared with each gene it touches and described by what it does to that transcript's exons: **whole gene** (a deletion or duplication covering it), **coding exon**, **non-coding exon / UTR**, or **intronic**, with the exon numbers hit (numbered along the transcript's strand, so a deletion of DMD exons 45-50 reports exactly those). Insertions and breakends are points; an inversion is a hit where a breakpoint falls, and a gene wholly inside the inverted segment is not. This is overlap with the gene model, not a predicted consequence: whether a deletion keeps the reading frame is not worked out.

### Population matches

Matching is by breakpoints and size, as svpack and Truvari do, not by exact allele: DEL / DUP / INV / CNV need a reciprocal overlap of at least 50%; INS must be within 500 bp with the shorter insertion at least half the length of the longer. Insertion sequences are not compared, so two different insertions at one locus can match. A DEL or DUP can also match a catalog CNV (whose direction the source does not state). When several records match, the best similarity is reported and the highest frequency among them. A call that matches no catalog is labelled *novel*, which is a statement about these four catalogs only.

## What the screens show

- **Overview**: counts read from the file (by type, by size, the FILTER values, breakends, imprecise calls); the samples and how the family structure was decided; the funnel from all calls to the review list, where every step is a real filter on the SV screen so its count can be reproduced there; the protein-coding genes hit by rare calls; and which evidence sources are installed.
- **Structural variants**: the table, with filters for type, size, genotype, inheritance and parental depth, supporting reads, GQ, gene impact, a gene list (case-insensitive), population frequency, "in no catalog", segmental duplication, repeats and special regions, and chromosome. By default it shows the proband's PASS calls that are not breakends, because the VCF holds the whole family's calls. A filter combination that returns nothing shows how many calls survive each filter in turn. A row opens to the call, the family's genotypes and read support, the gene hits with exon numbers and the MANE accession, the four population matches, the genome context, gene-level context, and a locus view drawing the SV over the genes' exons.

## Checked on the GIAB HG002 truth set

The full GRCh38 v5.0q set (89,995 deletions and insertions after dropping other types; run locally, not committed) indexed in under 6 seconds and reconciled with `bcftools` count for count: 39,081 deletions, 50,914 insertions, 73,694 with FILTER `PASS` or unset, none invalid or skipped. What it showed:

- **FILTER means something else here.** `HET1`, `HET2`, `GAP1`, `GAP2`, `DIPX` and `DIPY` are callability flags (heterozygous or uncalled in a haplotype of the assembly), not caller failures, and a call can carry two. The overview lists each value separately, so one call with `HET1;HET2` counts under both.
- **Hemizygous calls.** 1,964 calls carry a haploid genotype (1,653 on chrX, 311 on chrY). Before ploidy was read they were labelled heterozygous.
- **Population matches rise with size.** Of the PASS calls, 78% of deletions and 72% of insertions match one of the four catalogs, almost all through HPRC (long reads, so it sees the small events). gnomAD-SV matches about 9%, because short reads miss most events under a few hundred bp, so *novel against gnomAD-SV alone* is not meaningful for small SVs.
- **Most SVs are in repeats.** 94% of PASS deletions (30,595 of 32,551) and 90% of PASS insertions (36,993 of 41,143) are mostly repeat by the definition the filter uses, and the largest class is a tandem repeat (81% of deletions, 73% of insertions). Most structural variants do sit in repeats, and the label agrees with the truth set's own annotation where it can be compared: 94% of the calls the truth VCF flags `TRF` (49,260 of 52,185) are labelled tandem repeat here, and for deletions the Alu / L1 / LTR class agrees with its RepeatMasker class 90% / 76% / 59% of the time. The label is broader than the truth set's flag, though (76% of PASS calls labelled tandem repeat here, against 58% of all calls flagged `TRF`), because the UCSC simple-repeat and low-complexity tracks are included. So **Not mostly repeat** drops almost every call on a truth-set-like callset. It is a filter for screening false de novo calls, where repeats are the usual source, and the funnel puts it after the others (the checkbox says so). For an insertion the label is the repeat at the insertion point, since the inserted sequence is not compared.

## Limits

- No ClinGen dosage-sensitivity curation, so the ACMG/ClinGen CNV scoring is not applied; no CADD-SV; no needLR frequencies (all shown as not evaluated).
- Breakends are listed, not interpreted: a translocation or inversion breakpoint is not tied to its mate beyond the mate position, and most breakends are artifacts.
- The three callers are not merged, so the same event can appear in each caller's call set.
- Population catalogs are thin in segmental duplications, which is also where many calls are least reliable, so "novel" is common there. Use the segmental-duplication and repeat filters together with it.
- The MANE Select transcript is one transcript per gene; an event that matters in another isoform is not seen.
- The committed call set is the chr22 slice of an assembly-derived benchmark: it shows how the screens read a callset, not how sawfish, pbsv or sniffles behave. The suite has been run on real family call sets from those three callers, which are not in the repo.
