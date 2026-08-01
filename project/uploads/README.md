# project/uploads fixtures

## rexprt_scores.hg38.tsv (added 2026-08-01)

Real, published tandem-repeat pathogenicity predictions from **RExPRT**
(Danzi et al. 2024, "Deep learning-based approach for the classification of
pathogenic and benign tandem repeat expansions,"
[PMC10832122](https://pmc.ncbi.nlm.nih.gov/articles/PMC10832122/); code+data at
[github.com/ZuchnerLab/RExPRT](https://github.com/ZuchnerLab/RExPRT)), joined onto
`repeat_catalog.hg38.bed`'s 937,195 genome-wide regions. Not a re-implementation of
their model -- these are the lab's own precomputed ensemble SVM+XGBoost scores,
carried through unmodified.

### Provenance

- **Source**: `Reference TR scores/reference_TRsAnnotated_RExPRTscores.txt` from the
  RExPRT repo (Git LFS, 333MB, MIT license), which scores ~836,546 reference TRs
  (motif length 3-8bp) on **hg19** coordinates using real genomic-context features:
  GERP conservation, pLI/LOEUF gene constraint, exon/UTR/promoter location,
  transcription-factor-binding-site overlap (RAD21/SMC3/CTCF), TAD boundaries, eSTR
  (expression-associated) status, nervous-system expression, GC content, and Alu
  proximity.
- **Join**: `backend/scripts/real_data/build_rexprt_catalog_scores.py` lifts each
  scored TR's hg19 coordinates to hg38 (UCSC `hg19ToHg38` chain, via `pyliftover`)
  and matches by coordinate overlap against this project's own non-overlapping
  genome-wide catalog, keeping the best-overlap match per catalog region (8,933
  ties resolved out of 225,303 final matches -- the catalog and RExPRT's reference
  set weren't built from the identical repeat-finder run, so a handful of
  near-duplicate boundaries exist).
- **Coverage**: 225,303 of the catalog's 937,195 regions matched (832,547 of
  836,545 source rows lifted successfully; the remainder either failed liftover or
  had no overlapping catalog region). Only ~211k catalog regions have a 3-8bp
  motif at all -- RExPRT's trained range -- so the bulk of the catalog's
  regions (longer-motif VNTR-like repeats) were never scored by RExPRT and are
  absent from this file. **Absence means "not evaluated," never "benign."**
  `repeatscope/rexprt.py` and the API surface this explicitly rather than
  defaulting an unmatched region to any score.
- **Columns kept**: `motif`, `gene`, `gene_type`, `region` (exon/intron/intergenic),
  `location` (First/Middle/Last), `gerp`, `loeuf`, `pli`, `promoter`, `utr_5`,
  `utr_3`, `tad`, `estr`, `open_reg_region`, `rad21`, `smc3`, `ctcf`, `alu`,
  `alu_distance`, `gc_content`, `tissue`, `svm`, `xgb`, `ensemble_score`,
  `ensemble_max`, `ensemble_binary`, plus the original `hg19_chrom/start/end` for
  traceability back to RExPRT's own source row. RExPRT's raw S2SNet topological
  graph indices (`eSh*`/`eTr*`/`eH`/`eW`/`eS6`/`eS`/`J`/`eX*`) are real but not
  clinically interpretable on their own and are intentionally left out of this
  join -- the full source file is the reference if anyone needs them.

### Known limitation -- verified against a real example before use

RExPRT scores a repeat's **reference-length genomic context**, not any individual
sample's expanded allele -- and the paper documents that intronic pathogenic
repeats are prone to false negatives (RFC1, FXN, CNBP, ATXN10 misclassified as
benign in training, attributed to near-zero GERP scores typical of introns).

Spot-checking this join against **C9orf72** (the well-known FTD/ALS GGGGCC
expansion, STRchive `FTDALS1_C9orf72`, chr9:27,573,484-27,573,546) confirms this
is a live, not just theoretical, limitation: RExPRT's real score for that exact
locus is `ensemble_max=0.421`, `ensemble_binary=0` -- **classified benign**,
despite being one of the best-established pathogenic repeat expansions known.
The product surfaces RExPRT results on the genome-wide catalog as a **research
signal for candidate prioritization**, explicitly labeled as such, and never as a
replacement for the STRchive disease panel's own curated, length-threshold-based
clinical classification of its 81 known loci.
