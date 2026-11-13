# project/uploads fixtures

## repeat_catalog.hg38.bed (rescoped 2026-08-21: RExPRT-predicted-pathogenic only)

**Current scope**: 24,956 regions -- every genome-wide catalog region RExPRT
predicts pathogenic (`ensemble_binary=1` in `rexprt_scores.hg38.tsv`), and
nothing else. This is a deliberate, real reduction in scope, not a bug: the
full TRExplorer v2.0 catalog (5,339,786 regions, see below) is real, correct
data, but building the backend's in-memory genotype index from a matching
5.3M-record VCF took minutes and multiple GB of RAM under normal conditions,
and measurably longer (15+ minutes, still not finished) on a machine under
real memory pressure at the time. Since the actual use of the Repeat Catalog
screen was "browse candidate loci RExPRT flagged," not "browse all 5.3M
reference definitions," the *data* itself was narrowed to that candidate set
rather than filtering it at request time -- filtering after the fact still
requires building the full-size index first, so it wouldn't have helped.

Built by `backend/scripts/real_data/build_rexprt_pathogenic_subset.py` from
the full v2.0 catalog BED, RExPRT's own scores, and the matching full
genome-wide VCF (real `bcftools view -R` region extraction, not
reimplemented) -- verified to reproduce byte-identical output against the
manually-run extraction it replaced. Coverage stays 100% real TRID-matched
genotype (24,956 of 24,956 regions; the very small gap between this and the
24,941 RExPRT-positive coordinate count is a handful of real duplicate
-coordinate catalog entries with different IDs at the same position). The
`n_rexprt_evaluated`/`n_rexprt_predicted_pathogenic` API stats still report
against the *full* 727,004-row RExPRT reference table (unaffected by this
catalog rescoping -- see below), so the "predicted pathogenic only" UI filter
is now close to a no-op against this catalog specifically (nearly every row
already qualifies by construction).

The full-scope 5,339,786-region catalog and its matching real genome-wide
VCF are not kept in this repo after rescoping (they were ~536MB/278MB) but
remain fully recoverable: they're real Git LFS objects in the git history of
the commit that first added the TRExplorer v2.0 catalog (`git log -- 
project/uploads/repeat_catalog.hg38.bed`), or can be rebuilt from scratch per
the provenance below.

### Provenance of the full v2.0 catalog (git history; no longer deployed)

The genome-wide repeat catalog was the **`variation_clusters_and_isolated_TRs`**
catalog from the Broad Institute's TRExplorer v2.0 release (Weisburd,
Dolzhenko et al. 2026,
[biorxiv.org/content/10.1101/2024.10.04.615514](https://www.biorxiv.org/content/10.1101/2024.10.04.615514v5);
data at [github.com/broadinstitute/trexplorer-catalog](https://github.com/broadinstitute/trexplorer-catalog),
release [v2.0](https://github.com/broadinstitute/trexplorer-catalog/releases/tag/v2.0),
asset `TRExplorer.variation_clusters_and_isolated_TRs_v2.hg38.TRGT.bed.gz`,
MIT license) -- 5,339,786 real TR loci, hg38, TRGT `--repeats` BED format
(same `ID=;MOTIFS=;STRUC=` convention as the current, rescoped file).

The v2.0 release's headline "5,599,658 TR loci" figure (4,863,041 carried over
from v1.0 plus 736,617 new definitions in v2.0 -- known disease loci,
polymorphic TRs from T2T assemblies, HipSTR/Adotto/VamosV3 catalogs, ClinVar
indels that are TRs, and several published TR studies) describes the
release's separate **main** catalog file instead. This project used the
**variation-clusters** variant -- the one the Broad's own release notes
recommend for actual TRGT genotyping, since it merges some nearby/overlapping
simple-repeat calls into combined "variation cluster" regions -- because a
real genome-wide TRGT run for HG002 existed against this exact catalog,
giving 100% real TRID-matched genotype coverage rather than partial
coordinate-overlap fallback.

**Real, verified overlap in the data**: this catalog is *not* guaranteed
non-overlapping (unlike the pre-v2.0 catalog) -- direct measurement against
the release's main file found 922,396 rows (~16.5%) overlapping a preceding
region on the same chromosome, with runs of over 1,000 consecutive
overlapping entries in places; the variation-clusters variant consolidates
much of that, but real overlapping entries still occur. This surfaced a real
correctness bug in `TrgtIndex.find_overlap` (a fixed 5-element backward scan
that assumed sparse, non-overlapping data) and the equivalent logic in
`build_rexprt_catalog_scores.py` -- both fixed to use the same max-span
-bounded window technique `GeneIndex` already used correctly for overlapping
gene bodies. See the git history for `repeatscope/trgt.py`.

## HG002.GRCh38.trgt.genomewide.vcf.gz (rescoped 2026-08-21, alongside the catalog above)

A real, genuinely genome-wide TRGT v5.1.0 run for HG002, rescoped down to just
the records for the 24,956 regions above via real `bcftools view -R` region
extraction (see `build_rexprt_pathogenic_subset.py`) -- not resimulated or
re-genotyped. Every one of the 24,956 catalog regions has a real matching
genotype call (TRID-exact match verified directly). Before rescoping, this
was a full 5,339,786-record genome-wide run genotyped directly against the
full v2.0 catalog -- see git history for that provenance if the full run is
ever needed again.

## pangenome_repeat_catalog.hg38.bed (added 2026-10-26: the second catalog, "Pangenome Repeat Catalog")

5,709 regions -- the TRExplorer v2.0 variation-cluster loci that the pangenome SV paper's
Supplementary Table 3 (Lin et al., medRxiv
[10.64898/2026.08.21.26361050](https://www.medrxiv.org/content/10.64898/2026.08.21.26361050v1.full),
CC BY; tables 3a, 3b, 3d, 3e and 3f) *found* in TRExplorer. Same TRGT `--repeats` BED format and
ID convention as `repeat_catalog.hg38.bed`, and picked on the Repeat Catalog screen alongside it
(`GET /api/catalogs`, `GET /api/catalog?catalog=pangenome`).

Built by `backend/scripts/real_data/build_pangenome_catalog.py` from the found rows of the
comparison sheet `TRExplorer_v2_vs_SuppTable3_comparison` (the 8,048 found rows of the "In
Repeat Catalog" and "Not in Repeat Catalog" tabs) and the full v2.0 variation-cluster catalog.
Each found paper locus maps to the one variation-cluster locus that TRGT would genotype it as:

- the variation-cluster locus containing the locus the TRExplorer portal matched to the paper
  locus, by ID (a cluster's ID is the comma-joined IDs of the TRs merged into it);
- otherwise the variation-cluster locus overlapping the paper locus by the most bases (71 of the
  8,041 paper rows with coordinates: 67 had no portal locus at all, 4 a portal locus that is not in
  the variation-cluster catalog);
- otherwise left out: 2 of those 71 have no variation-cluster locus overlapping them (table 3e
  `chr17_41265399_41265765`, table 3f `chr9_137048431_137048637`), and 7 rows (MUC6, ACAN, ATXN8OS,
  HOXA13-I/II/III, in tables 3a/3d/3e) give a gene name and no coordinates. No definition was
  invented for any of them.

8,048 found rows -> 5,709 distinct loci, since many paper loci share a cluster. One locus is
the portal's match for a paper locus that abuts rather than overlaps it
(`11-10808754-10808761-GTC` for `chr11_10808761_10808828`, in tables 3a and 3f). 961 of the 5,709 are also in the
Repeat Catalog (those are its RExPRT-positive loci; the other 4,748 are new to this app). C9orf72's
GGCCCC repeat (`9-27573484-27573546-GGCCCC`) is in this catalog, which the Repeat Catalog leaves out
because RExPRT scores it benign.

Genotypes come from `backend/fixtures/hg002_real_pangenome/` (HG002, real `bcftools view -R` extraction
from the same full genome-wide TRGT run as the Repeat Catalog's VCF; every locus has a record and
TRGT called all but one, chr20:43,189,969-43,189,980). Genotypes are matched by TRID only: a locus
no TRGT record names is shown ungenotyped, never given the call of a neighbouring locus.
RExPRT scores and gene annotation are the same tables the Repeat Catalog uses. The full
variation-cluster catalog and genome-wide VCF the script reads are not kept in the repo; they are
the Git LFS objects of commit `5783ff5` (see above).

## rexprt_scores.hg38.tsv (added 2026-08-01; rebuilt 2026-08-14 for the v2.0 catalog)

Real, published tandem-repeat pathogenicity predictions from **RExPRT**
(Danzi et al. 2024, "Deep learning-based approach for the classification of
pathogenic and benign tandem repeat expansions,"
[PMC10832122](https://pmc.ncbi.nlm.nih.gov/articles/PMC10832122/); code+data at
[github.com/ZuchnerLab/RExPRT](https://github.com/ZuchnerLab/RExPRT)). Not a
re-implementation of their model -- these are the lab's own precomputed
ensemble SVM+XGBoost scores, carried through unmodified. **Unaffected by the
2026-08-21 catalog rescoping above** -- this remains the full reference table
(both pathogenic- and benign-predicted regions) matched against the full
5,339,786-region v2.0 catalog, independent of how many of those regions the
deployed `repeat_catalog.hg38.bed` currently includes.

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
  and matches by coordinate overlap against the full genome-wide catalog, keeping
  the best-overlap match per catalog region (8,602 ties resolved out of 727,004
  final matches against the v2.0 catalog -- the catalog and RExPRT's reference set
  weren't built from the identical repeat-finder run, so near-duplicate boundaries
  exist, more so given the v2.0 catalog's own real overlapping regions -- see
  above).
- **Coverage**: 727,004 of the full v2.0 catalog's 5,339,786 regions matched
  (832,547 of 836,545 source rows lifted successfully; the remainder either failed
  liftover or had no overlapping catalog region) -- up from 225,303 of 937,195
  against the old, smaller catalog, since a denser catalog gives RExPRT's fixed
  reference set more candidate regions to land on. The bulk of the catalog's
  regions (motifs outside RExPRT's 3-8bp trained range, or simply no matching
  source row) were never scored by RExPRT and are absent from this file.
  **Absence means "not evaluated," never "benign."** `repeatscope/rexprt.py` and
  the API surface this explicitly rather than defaulting an unmatched region to
  any score.
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
expansion, STRchive `FTDALS1_C9orf72`, chr9:27,573,484-27,573,546) confirmed
this is a live, not just theoretical, limitation: RExPRT's real score for the
overlapping v2.0 catalog region chr9:27,573,523-27,573,546 was
`ensemble_max=0.421`, `ensemble_binary=0` -- **classified benign**, despite
being one of the best-established pathogenic repeat expansions known. **As a
direct consequence, C9orf72 is absent from the current (RExPRT-positive-only)
`repeat_catalog.hg38.bed`** -- a real, visible cost of narrowing the catalog
to RExPRT's own predictions, not a hidden one. The full `rexprt_scores.hg38.tsv`
reference table still carries this real score for anyone who needs it. The
product surfaces RExPRT results on the genome-wide catalog as a **research
signal for candidate prioritization**, explicitly labeled as such, and never as a
replacement for the STRchive disease panel's own curated, length-threshold-based
clinical classification of its 81 known loci (which still includes C9orf72,
unaffected by any of this).
