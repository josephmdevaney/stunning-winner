# Epigenetics suite

The Epigenetics suite reads CpG methylation from a genome's PacBio pileup and reports what it shows at the human imprinted loci. It works from the genome alone: no parent samples, no array, no reference cohort.

| Part | What it shows | Status |
|---|---|---|
| **Imprinting** | 48 imprinted DMRs, which imprinting disorders their pattern fits, and a multi-locus check | Built |
| **Region Detail** | CpG-by-CpG methylation for one DMR, split by haplotype | Built |
| **Methods** | Thresholds, analysis window, data provenance, how each call is made | Built |
| **Episignatures / Signature Detail** | Scores for three conditions from EpiSignaLR's sequencing-derived models, and the DMRs behind each score | Built, needs a one-time setup (see [Episignatures](#episignatures-from-episignalr)) |
| **Array-derived classifiers** (on Episignatures) | EpiSignaLR's 34 array-derived classifiers, ranked, with each condition's CpGs against the reference profiles | Built, same setup |

Backend: `backend/epigenetics/` (routes under `/api/epigenetics`; episignatures in `backend/epigenetics/episignature/`). Frontend: `web/src/epigenetics/`. Tests: `backend/tests/test_epigenetics.py`, `test_epigenetics_api.py`, `test_episignature.py`.

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

## Episignatures from EpiSignaLR

[EpiSignaLR](https://github.com/millerlaboratory/EpiSignaLR) (Miller laboratory; a medRxiv preprint, October 2026) scores a long-read genome against methylation signatures of three conditions: **SETD1B-related disorder**, **Kleefstra syndrome (EHMT1)** and **Kabuki syndrome (KMT2D)**, and against 34 array-derived classifiers. The three signatures are each an ElasticNet-selected set of DMRs and a linear SVM, trained on nanopore genomes of blood DNA from 9 or 10 cases and 10 controls. This suite scores a registered genome with both.

### Setup

```
python backend/scripts/setup_episignalr.py
```

The EpiSignaLR repository publishes no licence, so none of it is committed. The script fetches each file from one pinned commit into `backend/epigenetics/episignature/vendor/` (gitignored), checks every file against a SHA-256 pinned in the script, and records what it fetched in `vendor/VERSION.json`. Until it has been run, the Episignatures screen says so and shows the command.

The models and the array reference tables are pickles, and a pickle can run code when it is loaded, so they are never loaded in the API process. `episignature/convert.py` reads each model in a separate `python -I` process that allows only numpy's array reconstruction, stands in for every scikit-learn class with an empty object that keeps the state it is given, and refuses anything else (`os.system` and the like). It needs only numpy and joblib; if the Python running the setup lacks them, the script makes a small virtualenv for the two inside the vendor folder. The output is plain JSON: the scaler, the selected features, the SVM's weights and the Platt sigmoid's two numbers. The converted weights are checked against the ones the authors publish in each model's metadata.

### How a genome is scored

1. **DMR values.** For each DMR in the model's panel, the pooled fraction Σ`est_mod_count` ÷ Σ`cov` over the DMR's listed CpGs that have at least 2 reads, from the combined pileup. This is EpiSignaLR's own PacBio calculation (`--weighted --min-coverage 2`).
2. **Score.** The selected DMRs are standardised with the training mean and spread, and the decision value is the SVM intercept plus one weight × value term per DMR. The case probability is scikit-learn's: a Platt sigmoid, refined by libsvm's two-class step, which is ported as written so the numbers match.
3. **Call.** 0.5 or more is a match; 0.4 to 0.6 is indeterminate; below is no signature. A DMR with no covered CpG is replaced by the training mean, as in EpiSignaLR; if that applies to more than a quarter of the selected DMRs, no call is made.

Each model's screen shows the decision value split into its per-DMR terms, and every selected DMR beside the training controls and cases, so a score can be checked against the data behind it.

### Validation

- **The probability** was compared with the authors' scikit-learn 1.8.0 pipeline on about 370 profiles per model (the training rows, random profiles, and profiles with missing DMRs): decision values and probabilities agree to 10⁻¹⁵. Tests keep a few of those numbers.
- **End to end**, the authors' own scripts (`region_average_methylation.py`, `run_enet_svm.py`) were run on HG002's pileup: every DMR value is identical and the scores match: SETD1B 0.0065, Kleefstra 0.5232, Kabuki 0.3602.

### Reading the scores

- **Platform.** The models are nanopore-trained. In the paper, PacBio controls sit at the edge of the control group nearest the cases, and PacBio libraries made 9 of the 12 indeterminate or false-positive calls among other groups' samples for the Kabuki and Kleefstra models. A match or an indeterminate score is a lead; a negative score is not clearance.
- **Pileup mode.** EpiSignaLR's PacBio pileups are made with `--pileup-mode count --modsites-mode reference`. A pileup in another mode (HG002's is `model`, `denovo`) gives different methylation values, so its scores are not comparable with the published ones. The screen says which mode the pileup was made in and what to re-run.
- **Sample type.** The signatures come from blood DNA (the paper notes dried-blood-spot DNA scoring differently). HG002's reference DNA is from a lymphoblastoid cell line, so it is not a like-for-like control. It scores Kleefstra 0.52 (indeterminate); with a cell line, a pileup in another mode and 20 training genomes, that is a reason to doubt the score, not a finding.
- **Training size.** The cohorts are 19 to 20 people. Kabuki's leave-one-out specificity is 50%: half of its controls were called cases in the authors' own test.

### The 34 array-derived classifiers

EpiSignaLR also scores a genome with 34 classifiers adapted from NSBEpi (Geysens et al.), built from published median array methylation (Aref-Eshghi et al.): one profile per condition and one for controls, over each condition's 69 to 153 signature CpGs (3,642 distinct CpGs). The conditions and genes are in the preprint's Supplemental Table 4 and in `episignature/data/array_conditions.json`.

- **CpG values.** A genome's methylation at each signature CpG is the percent methylated of the pileup row at that position (no coverage filter, as in EpiSignaLR); a CpG with no row is missing. All 3,642 positions are read in one `tabix -R` call.
- **The SVM is retrained for every genome.** EpiSignaLR adds the genome as an extra row to the condition's reference table, fills every gap (in the reference or the genome) with that CpG's mean over all rows, and fits a linear SVM with the condition's own profile as the single positive example (weighted 10 to 1). A genome's values reach the fit through any reference gap, so the fit cannot be stored once. `episignature/array_svm.py` repeats it and returns each condition's decision value; a score of 0.35 or more is a match, and more than one can match.
- **Where it runs.** Training needs scikit-learn, so the scorer is a separate process in a Python made by setup (the one running it if it has numpy, joblib and scikit-learn, otherwise a small virtualenv in the vendor folder; its path is recorded in `VERSION.json`). It reads only the converted reference JSON and a genome's values, never a pickle. Setup also converts the reference pickle (a dict of pandas tables) with `episignature/convert_array.py`, which never imports pandas and checks that the tables list their CpGs in the order of EpiSignaLR's coordinate file, the order that pairs a genome's values with the reference rows.
- **Cut-off.** EpiSignaLR's script and methods use 0.35; the preprint's results text says 0.3. 0.35 is used and echoed in every response.

Validation: on HG002's pileup, the authors' own `run_SVM.py` (scikit-learn 1.8.0) and this scorer (scikit-learn 1.6.1, the newest that supports Python 3.9) give the same score for all 34 conditions, within 10⁻¹³. A synthetic genome's scores from the authors' script are kept as golden numbers in the tests.

How they behaved in the preprint's genomes: the Kleefstra classifier called all 28 genomes of people with Kleefstra syndrome; the Kabuki classifier missed the two individuals with KDM6A variants; the SETD1B classifier called 1 of 3 PacBio genomes. A few controls scored above the cut-off for other conditions. With 34 classifiers an occasional match is expected.

HG002 reaches 1.02 for ICF2_3_4 (immunodeficiency-centromeric instability-facial anomalies syndrome), the only condition above the cut-off. That is not a finding. HG002's reference DNA is a lymphoblastoid cell line, and the genome sits near the ICF reference profile at many of its CpGs, but this was not tested further.

## Next

- **A pileup in count mode**, so scores can be compared with the paper's, if a genome's BAM is at hand and `aligned_bam_to_cpg_scores` is installed.
- **Reference ranges** from a control cohort, so a DMR value can be placed against the same platform's controls rather than the nanopore training set.
