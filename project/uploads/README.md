# project/uploads fixtures

## repeat_catalog.hg38.bed (replaced 2026-08-14: TRExplorer v2.0)

The genome-wide repeat catalog is the **`variation_clusters_and_isolated_TRs`**
catalog from the Broad Institute's TRExplorer v2.0 release (Weisburd,
Dolzhenko et al. 2026,
[biorxiv.org/content/10.1101/2024.10.04.615514](https://www.biorxiv.org/content/10.1101/2024.10.04.615514v5);
data at [github.com/broadinstitute/trexplorer-catalog](https://github.com/broadinstitute/trexplorer-catalog),
release [v2.0](https://github.com/broadinstitute/trexplorer-catalog/releases/tag/v2.0),
asset `TRExplorer.variation_clusters_and_isolated_TRs_v2.hg38.TRGT.bed.gz`,
MIT license) -- **5,339,786 real TR loci**, hg38, TRGT `--repeats` BED format
(same `ID=;MOTIFS=;STRUC=` convention as before). Replaces the previous
937,195-region catalog (itself extracted from a real HG002 genome-wide TRGT
run -- see `backend/fixtures/hg002_real_genomewide/README.md`), a ~5.7x
expansion.

The v2.0 release's headline "5,599,658 TR loci" figure (4,863,041 carried over
from v1.0 plus 736,617 new definitions in v2.0 -- known disease loci,
polymorphic TRs from T2T assemblies, HipSTR/Adotto/VamosV3 catalogs, ClinVar
indels that are TRs, and several published TR studies) describes the
release's separate **main** catalog file instead. This project uses the
**variation-clusters** variant -- the one the Broad's own release notes
recommend for actual TRGT genotyping, since it merges some nearby/overlapping
simple-repeat calls into combined "variation cluster" regions -- because a
real genome-wide TRGT run for HG002 exists against this exact catalog (see
below), giving **100% real TRID-matched genotype coverage** rather than
partial coordinate-overlap fallback. Real, verified data end-to-end took
priority over matching the release's headline row count exactly.

**Real, verified overlap in the data**: this catalog is *not* guaranteed
non-overlapping (unlike the old one) -- direct measurement against the
release's main file found 922,396 rows (~16.5%) overlapping a preceding
region on the same chromosome, with runs of over 1,000 consecutive
overlapping entries in places; the variation-clusters variant used here
exists specifically to consolidate much of that, but real overlapping entries
still occur. This surfaced a real correctness bug in `TrgtIndex.find_overlap`
(a fixed 5-element backward scan that assumed sparse, non-overlapping data)
and the equivalent logic in `build_rexprt_catalog_scores.py` -- both fixed to
use the same max-span-bounded window technique `GeneIndex` already used
correctly for overlapping gene bodies. See the git history for
`repeatscope/trgt.py`.

## HG002.GRCh38.trgt.genomewide.vcf.gz (replaced 2026-08-14)

A real, genuinely genome-wide TRGT v5.1.0 run for HG002, genotyped directly
against the catalog above -- see
`backend/fixtures/hg002_real_genomewide/README.md` for full provenance. Every
one of the catalog's 5,339,786 regions has a real matching genotype call in
this VCF (TRID-exact match verified directly, no coordinate-overlap fallback
needed for this pairing).

## rexprt_scores.hg38.tsv (added 2026-08-01; rebuilt 2026-08-14 for the v2.0 catalog)

Real, published tandem-repeat pathogenicity predictions from **RExPRT**
(Danzi et al. 2024, "Deep learning-based approach for the classification of
pathogenic and benign tandem repeat expansions,"
[PMC10832122](https://pmc.ncbi.nlm.nih.gov/articles/PMC10832122/); code+data at
[github.com/ZuchnerLab/RExPRT](https://github.com/ZuchnerLab/RExPRT)), joined onto
`repeat_catalog.hg38.bed`'s genome-wide regions. Not a re-implementation of
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
  and matches by coordinate overlap against the genome-wide catalog, keeping the
  best-overlap match per catalog region (8,602 ties resolved out of 727,004 final
  matches against the v2.0 catalog -- the catalog and RExPRT's reference set
  weren't built from the identical repeat-finder run, so near-duplicate boundaries
  exist, more so now given the v2.0 catalog's own real overlapping regions -- see
  above).
- **Coverage**: 727,004 of the v2.0 catalog's 5,339,786 regions matched (832,547 of
  836,545 source rows lifted successfully; the remainder either failed liftover or
  had no overlapping catalog region) -- up from 225,303 of 937,195 against the old,
  smaller catalog, since a denser catalog gives RExPRT's fixed reference set more
  candidate regions to land on. The bulk of the catalog's regions (motifs outside
  RExPRT's 3-8bp trained range, or simply no matching source row) were never scored
  by RExPRT and are absent from this file. **Absence means "not evaluated," never
  "benign."** `repeatscope/rexprt.py` and the API surface this explicitly rather
  than defaulting an unmatched region to any score.
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
is a live, not just theoretical, limitation: RExPRT's real score for the
overlapping v2.0 catalog region chr9:27,573,523-27,573,546 (this catalog's
variation-cluster consolidation combines that same repeat with adjacent
sub-repeats -- see above -- this region carries the identical underlying
RExPRT source row as every other coordinate framing of this locus checked) is
`ensemble_max=0.421`, `ensemble_binary=0` -- **classified benign**, despite
being one of the best-established pathogenic repeat expansions known.
The product surfaces RExPRT results on the genome-wide catalog as a **research
signal for candidate prioritization**, explicitly labeled as such, and never as a
replacement for the STRchive disease panel's own curated, length-threshold-based
clinical classification of its 81 known loci.
