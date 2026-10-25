# Epigenetics suite

The Epigenetics suite reads CpG methylation from a genome's PacBio pileup and reports what it shows at the human imprinted loci. It works from the genome alone: no parent samples, no array, no reference cohort.

| Part | What it shows | Status |
|---|---|---|
| **Imprinting** | 48 imprinted DMRs, which imprinting disorders their pattern fits, and a multi-locus check | Built |
| **Region Detail** | CpG-by-CpG methylation for one DMR, split by haplotype | Built |
| **Methods** | Thresholds, analysis window, data provenance, how each call is made | Built |
| **Episignatures** | Genome-wide DMR scoring with EpiSignaLR's classifiers | Not built (see [Next](#next-episignatures-from-episignalr)) |

Backend: `backend/epigenetics/` (routes under `/api/epigenetics`). Frontend: `web/src/epigenetics/`. Tests: `backend/tests/test_epigenetics.py`, `test_epigenetics_api.py`.

**Research use only.** A methylation pattern is not a diagnosis. The thresholds are defaults chosen for this panel and have not been clinically validated.

## Input

The suite reads the [pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools) bed files of a sample, with their tabix indexes:

- `<prefix>.combined.bed.gz` (required)
- `<prefix>.hap1.bed.gz` and `<prefix>.hap2.bed.gz` (optional; they come from a haplotagged BAM, for example after HiPhase)

Each row is one CpG: its methylation percent, coverage and the modified and unmodified counts. Region queries go through the `tabix` binary (htslib), so a DMR costs milliseconds and the bed files are never loaded whole.

Samples are registered in `backend/config/epigenetics.local.json`, which is gitignored (`epigenetics.local.json.example` shows the format). Either give a `dir` and `prefix`, or an explicit path per track. A sample with no `.tbi` index is not registered. Set `EPIGENETICS_CONFIG` to use a config elsewhere. Genomes stay where they are; nothing is copied into the repository.

The reader records what made the file: the pb-CpG-tools version, the pileup mode (`model` or `count`), the site mode and the minimum coverage, all read from the bed header and shown on Methods.

## The imprinted-DMR panel

`backend/epigenetics/data/imprinted_dmrs.json` is generated, not hand-written, by `backend/scripts/reference_data/build_imprinted_dmrs.py`. It fetches Table 1 of the 2024 interim joint statement on multi-locus imprinting disturbance ([Clin Epigenetics, PMC11295890](https://pmc.ncbi.nlm.nih.gov/articles/PMC11295890/), CC BY 4.0) and parses it. Re-run it to rebuild the panel.

The table has 48 DMRs in hg38, in three categories:

| Category | Count | Meaning |
|---|---|---|
| **Clinical** (CA) | 15 | Clinically associated: used to diagnose imprinting disorders |
| **Research** (NC) | 27 | Imprinted, but not used in diagnosis |
| **Polymorphic** (X) | 6 | Polymorphic or tissue-variable. Shown, never interpreted |

Each entry carries its coordinates, the published CpG count, the imprint origin (oocyte or sperm gDMR) and, where known, which parental allele is methylated (44 of 48). The seven loci the clinical rules lean on also carry a tighter **core window**, such as the SNRPN promoter CpG island, because pooling a whole 10 kb DMR dilutes the signal.

The MEG3/DLK1:IG-DMR coordinates differ between the 2022 and 2024 consensus papers. The entry keeps the 2024 ones and a `coordinate_notes` field says so; the UI marks the region as disputed.

## How a region is called

For each DMR the suite pools the CpGs in its analysis window as Σ modified ÷ Σ (modified + unmodified), so deeper sites count for more. This reproduces the per-DMR percentages in the reference spreadsheet exactly.

A call needs at least 5 CpGs and 8× mean depth (4× per haplotype). Otherwise the region is **low coverage** and nothing is claimed.

With haplotypes (phased):

| Pattern | Call |
|---|---|
| The two haplotypes differ by at least 40 points | **Imprint intact** (allele-specific methylation) |
| Both under 20% | **Loss of methylation** |
| Both over 80% | **Gain of methylation** |
| Anything else | **Borderline** |

Without haplotypes, the combined value is used: under 20% is a loss, over 80% a gain, 30–70% an intact imprint, anything else borderline. The three thresholds (loss, gain, allelic difference) can be changed on Methods; every response echoes the ones it used. A threshold change re-classifies from cached CpGs, so it is instant.

A polymorphic region is reported as "Not interpreted" whatever it shows. A research region's finding is outlined, not filled, so it is not read as diagnostic.

## Disorders and the multi-locus check

`disorders.py` holds one rule per disorder, over the clinical DMRs only:

| Disorder | Pattern |
|---|---|
| Prader-Willi | SNURF:TSS-DMR gain |
| Angelman | SNURF:TSS-DMR loss |
| Beckwith-Wiedemann | KCNQ1OT1:TSS-DMR loss, or H19/IGF2:IG-DMR gain |
| Silver-Russell | H19/IGF2:IG-DMR loss (MEST and GRB10 gain are supportive) |
| Transient neonatal diabetes | PLAGL1:alt-TSS-DMR loss |
| Pseudohypoparathyroidism 1B | GNAS A/B loss (XL and AS1 loss, NESP gain are supportive) |
| Temple (14q32) | MEG3:TSS-DMR loss (IG-DMR loss is supportive) |
| Kagami-Ogata (14q32) | MEG3:TSS-DMR and IG-DMR gain |

Each gets one of **Pattern present**, **Supportive only** (only the supporting loci fit), **No pattern** or **Insufficient data** (a needed locus could not be called). A multi-locus disturbance is reported when clinical DMRs on **two or more chromosomes** are aberrant: it counts chromosomes, not DMRs, because several DMRs in one cluster are one event.

What a methylation pattern cannot do is on every screen: it cannot tell a deletion, uniparental disomy and an imprinting defect apart, tissue and mosaicism change what blood shows, and a "pattern present" is a lead for a validated test. For the UPD question, use the ROH suite alongside this one.

## API

| Route | Returns |
|---|---|
| `GET /api/epigenetics/samples` | Registered samples and which tracks each has |
| `GET /api/epigenetics/imprinting` | Regions, disorders, multi-locus check, provenance, thresholds. Params: `sample`, `window` (`core` or `dmr`), `loss_below`, `gain_above`, `allelic_min_diff` |
| `GET /api/epigenetics/imprinting/region` | CpG-level points for one DMR and its flank. Params: `sample`, `id`, `flank` (0–50000), `window` |

A bad threshold or window is a 422, an unknown sample or DMR a 404, a missing `tabix` a 503 with the reason.

## Validation

- **Reference spreadsheet.** All 24 per-DMR values in the sheet that came with the EpiSignaLR paper were reproduced exactly from HG002's beds.
- **HG002** (public reference genome, Revio, pb-CpG-tools 3.0.0): 41 of 48 regions show an intact imprint, 13 of the 15 clinical ones. All eight disorders read "No pattern", and the multi-locus check is negative (14 of 15 clinical DMRs assessable). The exceptions are all explained: the MEG3/DLK1 IG-DMR is borderline and its coordinates are disputed; IGF2R:Int2 and the SNRPN alternative TSS are borderline (the latter is an outlier in the published data too); MKRN3 is polymorphic; PPIEL:Ex1 and GPR1-AS read as both-alleles-methylated, but they are a research and a polymorphic region; IGF2:Ex9 has too little coverage.
- **Not validated.** No genome with a known imprinting disorder was available, so the positive case (a real PWS, BWS or Temple genome) has been tested only on synthetic data. The thresholds have not been calibrated against a cohort.

## Limits

- The panel comes from one consensus paper; a locus it omits is not assessed.
- Whole-genome pileups made in `model` mode are not the same scale as the array-derived or `count`-mode data other tools were trained on.
- Mosaicism below the thresholds' sensitivity, and tissue, are not modelled.

## Next: episignatures from EpiSignaLR

Not built yet. The plan, from the paper and repository:

- EpiSignaLR scores genome-wide DMR methylation against 34 episignature classifiers (elastic net and linear SVM, with probabilities; 0.5 positive, 0.4–0.6 indeterminate).
- Its repository has no license, so its code and models are **fetched at setup into a gitignored directory and never committed**.
- The models are pickled `joblib` files. Loading a pickle runs code, so scoring would run in a separate subprocess from a pinned checkout, never in the API process.
- PacBio controls sit shifted toward the cases in the paper, and the models were trained on count-mode data. Every episignature result would carry both warnings.
