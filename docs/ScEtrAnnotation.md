# Expression effects on the disease panel (TenK10K sc-eTR)

The locus detail screen for a disease-panel locus can show an **Expression effects in immune cells** card. It reports, for the STRchive disease repeat, whether the repeat was polymorphic in the TenK10K cohort and which repeats at or near it were the lead expression association for a gene, with the effect in each immune cell type.

It is a research annotation about *common* repeat-length variation and gene expression. It says nothing about pathogenic expansions, and it never feeds the locus classification.

Backend: `backend/repeatscope/sc_etr.py` (loader), `GET /api/locus/{id}` (`expression_effects`). Build script: `backend/scripts/real_data/build_sc_etr_disease_annotation.py`. Frontend: `web/src/components/ExpressionEffectsCard.tsx`. Tests: `backend/tests/test_sc_etr.py` (synthetic data only).

## Source

Tanudisastro et al., *Tandem repeat variation shapes immune cell type-specific gene expression* (bioRxiv 10.1101/2024.11.02.621562, v4, 2026-03-03): 1,925 people (Tasmanian Ophthalmic Biobank and BioHEART), ExpansionHunter genotypes at about 1.7 million repeats, associaTR mapping against pseudobulk expression in 28 immune cell types. The released tables are on Zenodo, 10.5281/zenodo.15009519.

**Licence: CC BY-NC-ND 4.0** (the paper and the Zenodo record). The repo therefore holds no data from it. The built annotation file is read from outside the repo and the card is absent when the file is.

## Building the annotation

Unzip two files from the Zenodo record, `qval_egenes.zip` (54 MB) and `polymorphic_catalog.zip` (97 MB, a 2.2 GB JSON), then:

```bash
python3 backend/scripts/real_data/build_sc_etr_disease_annotation.py \
    --qval-dir /path/to/qval_egenes \
    --catalog-json /path/to/filtered_polymorphic_catalog.json \
    --strchive project/uploads/STRchive-loci.json \
    --out ~/sc_etr_annotations/disease_panel.hg38.json
```

It takes about 25 seconds. The API reads `~/sc_etr_annotations/disease_panel.hg38.json` by default, or the path in `REPEATSCOPE_SC_ETR_ANNOTATIONS`. A rebuilt file is picked up without a restart.

### How the join works

- Each `qval_egenes` table is one cell type. A row is one tested gene, with list-valued columns holding its lead repeat (or tied leads) and a **gene-level** q-value. Rows with q < 0.05 are kept. `CD4_TCM_permuted_qval.tsv` is the permutation null and is skipped.
- A lead is given only as `chr`, `pos` and `motif`. `pos` is the **0-based start**, the same as the catalog's `ReferenceRegion` start and as this project's BED files. Each lead is matched to the catalog repeat that starts there and has its motif, which gives its end. Multi-repeat catalog loci are split into their component repeats first, because each component is its own tested variant.
- A lead belongs to a STRchive locus if it overlaps it (`overlaps`) or lies within 1 kb (`adjacent`, a different repeat next door). `same_motif` compares the lead's motif with the locus's motifs up to rotation and strand, since STRchive writes CAG where the catalog writes CTG.
- `polymorphic` means a repeat of the paper's catalog overlaps the locus.

### Result of the real build

All 60,571 lead repeats place on the genome, and the 15,889 genes with a lead match the paper's count. For the 82 panel loci: 78 are polymorphic in the cohort, 12 have a lead repeat at the disease repeat, and 5 more have one only nearby.

Five of the paper's seven fine-mapped disease loci come out as a lead at their own repeat, regulating their own gene: C9orf72, CSTB, DIP2B, GLS and PPP2R2B. CACNA1A's lead is 981 bp from the disease repeat and regulates a neighbouring gene (STX10); CBL has no lead at its repeat here.

## What the card does and does not say

- **Lead only.** The released tables hold the top repeat for each gene in each cell type. A repeat that is not listed is not shown to have no effect. (The 15,889 genes match the paper, but the tables give 60,571 distinct lead repeats against the paper's 69,210 sc-eTRs; the paper probably counts non-lead associations too.)
- **Lead is not causal.** Fine-mapping results are not in the Zenodo files and are not shown.
- **Effect.** The fixed-effects meta-analysis coefficient per summed diploid repeat length, on normalized expression. The q-value is gene-level, not per repeat. Direction is shown as ▲/▼ and in blue/purple; the red/green of the status palette is not used, because an expression change is neither pathogenic nor normal.
- **Neighbouring repeats are labelled.** HTT's lead, for example, is the CCG repeat 12 bp from the CAG disease repeat, and is shown as a different repeat. A lead whose motif differs from the disease motif but overlaps the tract (CSTB's pentamer inside its dodecamer repeat) says so.
- **Common range only.** These are repeat lengths seen in European-ancestry cohorts (950 people with no sign of ocular disease and 975 with diagnosed or suspected cardiac disorders), not expansions and not people with the disease.

## Limits

- Per-repeat results for non-lead repeats are in the per-cell-type meta-analysis archives on Zenodo (4.6-7.9 GB each) and are not used.
- Supplementary Tables 1-11 on bioRxiv (fine-mapping, disease loci, PanelApp) were not used.
- Only the 82 STRchive loci are annotated. About 1.3% of the Repeat Catalog (326 of 24,956 loci) and 1.6% of the Pangenome Repeat Catalog (94 of 5,709) have a lead at the same locus or with the same motif, so a catalog-wide column would be sparse.
