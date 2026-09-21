# RepeatScope — Runs of Homozygosity (ROH): Parts & Methods

This document describes the ROH module of RepeatScope: its architecture (backend, API, frontend, tests) and its analytical methodology (parsing, region classification, gene overlap, the inbreeding-coefficient estimate, and the UPD screening heuristic), including the caveats around each.

> Reference build: `backend/repeatscope/roh.py`, `roh_reference.py`, `api/main.py`, `tests/test_roh.py`; frontend `web/src/api.ts`, `web/src/screens/ROHScreen.tsx`.

## Table of Contents

1. [Overview](#1-overview)
2. [Inputs & Data Sources](#2-inputs--data-sources)
3. [Parts (Architecture)](#3-parts-architecture)
   - 3.1 [Backend module — `backend/repeatscope/roh.py`](#31-backend-module--backendrepeatscoperohpy)
   - 3.2 [Reference data — `backend/repeatscope/roh_reference.py`](#32-reference-data--backendrepeatscoperoh_referencepy)
   - 3.3 [API layer — `backend/repeatscope/api/main.py`](#33-api-layer--backendrepeatscopeapimainpy)
   - 3.4 [Frontend client — `web/src/api.ts`](#34-frontend-client--websrcapits)
   - 3.5 [UI — `web/src/screens/ROHScreen.tsx`](#35-ui--websrcscreensrohscreentsx)
   - 3.6 [Tests — `backend/tests/test_roh.py`](#36-tests--backendteststest_rohpy)
4. [Methods (Methodology)](#4-methods-methodology)
   - 4.1 [Parsing](#41-parsing)
   - 4.2 [Region position classification](#42-region-position-classification)
   - 4.3 [Length banding](#43-length-banding)
   - 4.4 [Gene overlap (autozygosity mapping)](#44-gene-overlap-autozygosity-mapping)
   - 4.5 [Autosomal inbreeding coefficient (F) estimate](#45-autosomal-inbreeding-coefficient-f-estimate)
   - 4.6 [UPD screening heuristic](#46-upd-screening-heuristic)
   - 4.7 [AR candidate gene surfacing](#47-ar-candidate-gene-surfacing)
   - 4.8 [Consanguinity verdict](#48-consanguinity-verdict)
   - 4.9 [Relationship F-value reference table](#49-relationship-f-value-reference-table)
   - 4.10 [Non-autosomal handling (chrX, chrY, chrM)](#410-non-autosomal-handling-chrx-chry-chrm)
5. [API Reference](#5-api-reference)
6. [Response Schema — Selected Fields](#6-response-schema--selected-fields)
7. [Limitations & Caveats](#7-limitations--caveats)
8. [Test Coverage](#8-test-coverage)

## 1. Overview

ROH (Runs of Homozygosity) is a self-contained analysis module within RepeatScope that surfaces autozygosity, consanguinity signal, and possible uniparental disomy (UPD) from real `bcftools roh` output. Unlike the rest of RepeatScope — which is organized around a repeat-panel "mode" and sample — ROH is keyed only by sample: homozygosity mapping is a whole-genome SNP-based analysis and does not depend on a tandem-repeat catalog/panel at all.

Currently one real sample has an ROH run: HG002, from real `bcftools roh` RG-region output. There is no synthetic or demo sample — by design, the module surfaces real autozygosity/consanguinity/UPD calls only, not placeholder data.

## 2. Inputs & Data Sources

- **Primary input:** `bcftools roh` RG (region-summary) lines — tab-separated `sample, chrom, start, end, length_bp, n_markers, quality`. ST (per-marker) lines and header/comment lines are ignored.
- **Fixture:** `backend/fixtures/roh/HG002.GRCh38.roh.txt` (sample HG002, genome build GRCh38).
- **Static reference data** (`roh_reference.py`): GRCh38 chromosome lengths and centromere midpoints (UCSC/NCBI), a curated OMIM disease-gene model, UPD-screening length thresholds, and theoretical relationship-based F values.

## 3. Parts (Architecture)

The module spans four layers: a backend analysis module, a small API surface, a typed frontend client, and a dedicated screen with several chart components.

### 3.1 Backend module — `backend/repeatscope/roh.py`

Pure-Python module (293 lines) that parses real `bcftools roh` output and derives all downstream analysis from those real regions only — no simulated values are introduced at any step.

| Function | Purpose |
|---|---|
| `parse_roh_file(path)` | Parses RG (region-summary) lines from a real bcftools roh output file into `RohRegion` records; ignores ST lines and comments; drops a known stray empty tab-separated field rather than failing. |
| `classify_position(chrom, start, end)` | Classifies a region as telomeric / pericentromeric / interstitial using real GRCh38 chromosome lengths and centromere midpoints. Returns `None` for chromosomes outside that reference (chrY, chrM). |
| `length_band(length_bp)` | Buckets a region's length into `<1 Mb` / `1–3 Mb` / `3–5 Mb` / `≥5 Mb` for display/coloring. |
| `gene_overlaps(chrom, start, end)` | Returns curated OMIM disease-gene model entries whose coordinates overlap the region. |
| `upd_flag(chrom, position, length_bp)` | Screening-only UPD flag using the SNP-microarray-style length/position thresholds in `roh_reference.py`. Autosomes only. |
| `_region_dict(r)` | Assembles the full per-region dict served by the API: adds position class, gene overlaps, length band, imprinted flag, UPD flag, and F-value eligibility to a parsed `RohRegion`. |
| `build_regions(sample, path)` / `regions_for_sample(sample)` | Builds (and, for the cached variant, memoizes) the full list of region dicts for a sample. |
| `build_overview(sample, path)` / `overview_for_sample(sample)` | Aggregates per-sample summary statistics: F-value, UPD flags, AR-candidate genes, per-chromosome rollup, and the relationship F-value reference table. |
| `sample_keys()` | Returns the list of sample IDs that have a real ROH run (currently just HG002). |
| `query_regions(sample, chrom, min_len_mb, page, page_size)` | Filters and paginates a sample's regions for the UI's regions table. |

### 3.2 Reference data — `backend/repeatscope/roh_reference.py`

Static constants module (144 lines): assembly geometry, screening thresholds, and the curated gene model, documented in its own module docstring as to what is a public reference fact versus a user-curated / heuristic value.

| Constant | Description |
|---|---|
| `CHROM_LENGTH` | GRCh38 chr1–22 + chrX + chrY base-pair lengths (standard assembly values). |
| `CENTROMERE_MID` | Approximate GRCh38 centromere midpoints (from UCSC cytoBand ranges); used only to detect centromere-spanning regions, not exact arm boundaries. |
| `AUTOSOMES` / `CHROM_ORDER` | chr1–22 in karyotype order, then chrX / chrY / chrM for reporting order. |
| `IMPRINTED_CHROMS` | `{chr6, chr7, chr11, chr14, chr15, chr20}` — autosomes carrying a well-established imprinted locus, used by the UPD heuristic. |
| `AUTOSOMAL_GENOME_LENGTH` | Sum of `CHROM_LENGTH` over the autosomes; denominator for the F-value estimate. |
| `TELOMERE_PROXIMITY_BP` | 1,000,000 — distance from a chromosome terminus counted as telomeric. |
| `F_VALUE_MIN_LENGTH_BP` | 3,000,000 — minimum ROH length counted toward the autosomal F-value (IBD-informative floor). |
| `AR_CANDIDATE_MIN_LENGTH_BP` | 1,000,000 — separate, lower floor for surfacing an AR gene as an autozygosity-mapping candidate. |
| `UPD_TELOMERIC_IMPRINTED_BP` | 5,000,000 — telomeric UPD-screening threshold. |
| `UPD_INTERSTITIAL_IMPRINTED_BP` | 10,000,000 — interstitial threshold on an imprinted chromosome. |
| `UPD_INTERSTITIAL_NONIMPRINTED_BP` | 15,000,000 — interstitial threshold elsewhere. |
| `RELATIONSHIP_F_VALUES` | Theoretical F = (1/2)^n by pedigree relationship, for the consanguinity ladder reference points. |
| `GENE_MODEL` | The 12-entry curated OMIM disease-gene table (see Section 4.4). |

### 3.3 API layer — `backend/repeatscope/api/main.py`

Three routes, added as a self-contained section of the FastAPI app (not scoped by the `mode`/panel query parameters the rest of RepeatScope's API uses). A helper, `_resolve_roh_sample`, defaults to the first available sample when none is given and returns 404 for an unrecognized sample.

| Route | Params | Returns |
|---|---|---|
| `GET /api/roh/samples` | — | `{ samples: string[], default: string }` |
| `GET /api/roh/overview` | `sample?` (defaults to first available) | `RohOverviewResponse` — `by_chromosome`, `stats`, `upd_flags`, `ar_candidate_genes`, `relationship_f_values`, methodology notes |
| `GET /api/roh/regions` | `sample?, chrom (default "all"), min_len_mb, page, page_size` | `RohRegionsResponse` — filtered/paginated region rows |

An unrecognized sample value returns HTTP 404 on both `/api/roh/overview` and `/api/roh/regions`.

### 3.4 Frontend client — `web/src/api.ts`

TypeScript interfaces mirror the backend response shapes exactly: `RohGeneHit`, `RohRegion`, `RohByChromosome`, `RohArCandidateGene`, `RohRelationshipFValue`, `RohOverviewResponse`, `RohRegionsResponse`, `RohSamplesResponse`. The `rohApi` client exposes `samples()`, `overview(sample)`, and `regions(params)`, each a thin fetch wrapper returning the typed response.

### 3.5 UI — `web/src/screens/ROHScreen.tsx`

A 785-line screen assembled from the following components:

- **ROHScreen** — top-level screen; loads the overview and full region list on sample change and lays out the stat row, UPD banner, ideogram, consanguinity card, AR-candidate list, per-chromosome track, and regions table.
- **StatGrid / StatCard** — shared stat-tile components: autosomal F (%), autosomal ROH total, UPD flag count, AR candidate gene count.
- **IdeogramCard / Ideogram** — genome-wide autozygosity map: all 24 chromosomes drawn to genomic scale in two columns, ROH regions painted by length band, with centromere, imprinted-chromosome, and UPD-flag legend markers; filterable by length band.
- **PerChromosomeTrackCard / ChromosomeTrack** — a zoomed single-chromosome track: ROH segments to scale, centromere marker, an Mb ruler, and OMIM disease-gene tick marks colored by inheritance (AR/AD).
- **ConsanguinityCard / ConsanguinityLadder / consanguinityVerdict** — F at the ≥1 Mb and ≥3 Mb floors, a horizontal ladder plotting the sample's F against theoretical relationship-degree F values, and an ACMG-threshold-based verdict banner (indicated / possible / not indicated).
- **RegionsTable / Field / PageButton** — a filterable, paginated table of individual ROH regions: chromosome, coordinates, length, position class, marker count, quality, and overlapping genes/UPD flag.

### 3.6 Tests — `backend/tests/test_roh.py`

184 lines, 20 tests, all run against the real HG002 fixture (no mocked ROH data). See [Section 8](#8-test-coverage) for the full breakdown.

## 4. Methods (Methodology)

### 4.1 Parsing

Only RG (region) lines are used; ST (per-SNP state) lines are ignored entirely. A known paste artifact in the HG002 fixture (a stray empty tab-separated field on one chr6 row) is defended against by dropping empty fields rather than failing the parse — a genuine RG row never contains one, and the recovered `length_bp` is verified to equal `end − start + 1` like every other row.

### 4.2 Region position classification

Each region is classified against real GRCh38 chromosome lengths and approximate centromere midpoints (from UCSC cytoBand ranges):

- **Telomeric** — region start or end within 1 Mb of a chromosome terminus.
- **Pericentromeric** — region spans the chromosome's centromere midpoint.
- **Interstitial** — neither of the above.
- **Unclassified (null)** — chromosomes without reference length data in the model (chrY, chrM).

### 4.3 Length banding

Each region is bucketed into `<1 Mb`, `1–3 Mb`, `3–5 Mb`, or `≥5 Mb`, used for ideogram/track color-coding and the regions-table length filter.

### 4.4 Gene overlap (autozygosity mapping)

Each region is checked for overlap against a curated, user-verified OMIM disease-gene model (coordinates cross-checked against NCBI Gene / OMIM), tagged autosomal recessive (AR) or autosomal dominant (AD):

| Symbol | Chrom | Start | End | Inh. | Disease | OMIM |
|---|---|---|---|---|---|---|
| ZAP70 | chr2 | 97,713,576 | 97,756,364 | AR | Immunodeficiency 48 (SCID) | 176947 |
| CNGA3 | chr2 | 98,346,188 | 98,399,797 | AR | Achromatopsia 2 | 600053 |
| HFE | chr6 | 26,087,429 | 26,098,343 | AR | Hemochromatosis, type 1 | 613609 |
| SHH | chr7 | 155,799,980 | 155,812,463 | AD | Holoprosencephaly 3 | 600725 |
| MNX1 | chr7 | 156,994,915 | 157,004,509 | AD | Currarino syndrome | 142994 |
| TJP2 | chr9 | 69,121,264 | 69,255,208 | AR | Cholestasis, progressive familial intrahepatic 4 | 607709 |
| TMC1 | chr9 | 72,521,608 | 72,838,297 | AR | Deafness, DFNB7/11 (AR); DFNA36 (AD) | 606706 |
| RB1 | chr13 | 48,303,751 | 48,481,890 | AD | Retinoblastoma | 614041 |
| OTOA | chr16 | 21,689,000 | 21,819,000 | AR | Deafness, AR 22 | 607038 |
| COG7 | chr16 | 23,368,000 | 23,432,000 | AR | Congenital disorder of glycosylation IIe | 606978 |
| CLN3 | chr16 | 28,466,000 | 28,490,000 | AR | Neuronal ceroid lipofuscinosis 3 (Batten) | 607042 |
| TUFM | chr16 | 28,850,000 | 28,857,000 | AR | Combined oxidative phosphorylation deficiency 4 | 602389 |

### 4.5 Autosomal inbreeding coefficient (F) estimate

**F = (total length of autosomal ROH ≥ 3 Mb) ÷ (total autosomal genome length)**

The 3 Mb floor (`F_VALUE_MIN_LENGTH_BP`) excludes short background ROH arising from linkage disequilibrium or genetic drift, which are common even in outbred individuals and are not informative of identity-by-descent (IBD). This is a conventional analytical floor, not a value pinned to a specific external citation. chrX, chrY, and chrM are excluded from this calculation entirely (see [4.10](#410-non-autosomal-handling-chrx-chry-chrm)).

### 4.6 UPD screening heuristic

A region is flagged as a possible UPD candidate under the following autosome-only rules:

| Condition | Result |
|---|---|
| Telomeric, any autosome, ≥ 5 Mb | Flagged (imprinted-chromosome status noted if applicable) |
| Interstitial, imprinted chromosome (chr6/7/11/14/15/20), ≥ 10 Mb | Flagged |
| Interstitial, non-imprinted autosome, ≥ 15 Mb | Flagged |
| chrX / chrY / chrM, any length | Never flagged (heuristic is autosome-only) |

> These thresholds follow a commonly cited SNP-microarray ROH/UPD screening convention, often attributed to Hoppman-Chaney et al. **That attribution has NOT been independently verified against a primary source in this codebase.** Flags are explicitly surfaced (in both the API response and the UI) as screening candidates for orthogonal confirmation — methylation-specific analysis, microsatellite or trio testing — not as diagnostic calls.

### 4.7 AR candidate gene surfacing

Autosomal ROH ≥ 1 Mb (`AR_CANDIDATE_MIN_LENGTH_BP`) that overlap an AR-inheritance gene in the model are surfaced as autozygosity-mapping candidates. This floor is deliberately lower than — and kept separate from — the F-value's 3 Mb IBD-informative floor: a real AR disease-causing variant can sit within a smaller genuinely-autozygous ROH, and applying the F-value's stricter floor here would silently drop real candidates. Where more than one qualifying region overlaps the same gene, the largest overlapping region is kept.

### 4.8 Consanguinity verdict

The UI derives an interpretive verdict from F(≥3 Mb), framed against ACMG 2021 technical-standard consanguinity-reporting conventions:

- **F(≥3 Mb) < 2%** — not indicated.
- **F(≥3 Mb) 2–<10%** — possible (within the ACMG 2–3% consanguinity reporting range; consider discussing close-relation history).
- **F(≥3 Mb) ≥ 10%** — indicated (meets or exceeds first-cousin-level autozygosity; consider discussing history and an expanded AR candidate-gene review).

A ladder visualization plots the sample's F(≥1 Mb) and F(≥3 Mb) alongside theoretical F values for named pedigree relationships, for visual context.

### 4.9 Relationship F-value reference table

Theoretical inbreeding coefficient F = (1/2)^n, for n meioses separating the parents' common ancestor — standard textbook population-genetics values, not a disease-specific guideline:

| Relationship | F | F % |
|---|---|---|
| Parent-child / full siblings | 0.25 | 25% |
| Uncle-niece / half-siblings / double first cousins | 0.125 | 12.5% |
| First cousins | 0.0625 | 6.25% |
| First cousins once removed | 0.03125 | 3.125% |
| Second cousins | 0.015625 | 1.5625% |
| Second cousins once removed | 0.0078125 | 0.78125% |
| Third cousins | 0.00390625 | 0.390625% |

### 4.10 Non-autosomal handling (chrX, chrY, chrM)

chrX, chrY, and chrM regions are real `bcftools roh` calls and are reported in full, but excluded from the autosomal F-value: chrX/chrY homozygosity in a male sample reflects hemizygosity, not autozygosity, and mitochondrial DNA has no diploid homozygosity concept to begin with.

## 5. API Reference

All three routes are unscoped by RepeatScope's mode/panel query parameters (see [3.3](#33-api-layer--backendrepeatscopeapimainpy)); the only relevant key is `sample`.

| Route | Params | Returns |
|---|---|---|
| `GET /api/roh/samples` | — | `{ samples: string[], default: string }` |
| `GET /api/roh/overview` | `sample?` | `RohOverviewResponse` |
| `GET /api/roh/regions` | `sample?, chrom, min_len_mb, page, page_size` | `RohRegionsResponse` |

An unrecognized `sample` returns HTTP 404 on `/api/roh/overview` and `/api/roh/regions`.

## 6. Response Schema — Selected Fields

`RohOverviewResponse.stats` (the aggregate summary block returned by `/api/roh/overview`):

| Field | Type | Description |
|---|---|---|
| `n_regions_total` / `n_regions_autosomal` / `n_regions_chrX/Y/M` | int | Region counts overall and by chromosome class. |
| `n_regions_f_eligible` | int | Autosomal regions ≥ `F_VALUE_MIN_LENGTH_BP`. |
| `total_autosomal_roh_bp` / `total_f_eligible_roh_bp` | int | Summed region length in bp, all autosomal vs. F-eligible only. |
| `f_value` / `f_value_pct` | float | Autosomal inbreeding coefficient estimate, as a fraction and as a percentage. |
| `n_upd_flags` | int | Count of regions flagged by the UPD screening heuristic. |
| `n_ar_candidate_genes` | int | Count of distinct AR genes surfaced as autozygosity-mapping candidates. |
| `f_value_min_length_mb` / `ar_candidate_min_length_mb` | int | The two length floors (3 Mb / 1 Mb), echoed back for UI display. |

`RohRegion` (one row of `/api/roh/regions`, and each entry in `upd_flags`):

- `sample, chrom, start, end, length_bp, length_mb, n_markers, quality` — as parsed/derived from the bcftools roh RG line.
- `position` — telomeric / pericentromeric / interstitial / null.
- `length_band, imprinted_chrom, autosomal, f_value_eligible` — derived classification flags.
- `genes` — array of overlapping OMIM gene-model entries (symbol, inheritance, disease, OMIM ID, coordinates).
- `upd_flag` — the screening message string, or null.

## 7. Limitations & Caveats

- Single real sample supported (HG002) — by design, there is no synthetic/demo ROH sample; the tool surfaces real calls only.
- UPD screening thresholds: the attribution to Hoppman-Chaney-style SNP-microarray convention has not been independently verified against a primary source in this codebase. Treat flags as screening candidates, not diagnostic calls.
- Centromere positions are approximate midpoints (from UCSC cytoBand ranges), used only to detect centromere-spanning regions — not exact acrocentric/metacentric arm boundaries.
- The OMIM disease-gene model is a curated 12-gene subset (chr2, 6, 7, 9, 13, 16), not a comprehensive AR/AD gene panel.
- The F-value and AR-candidate floors (3 Mb / 1 Mb) are conventional analytical thresholds, not values pinned to a single external citation.

## 8. Test Coverage

20 tests in `backend/tests/test_roh.py`, all exercised against the real HG002 fixture (8,966 parsed RG regions):

| Area | What's covered |
|---|---|
| Parsing | Real RG-line parsing (8,966 regions for HG002), stray-empty-field recovery, comment-line exclusion. |
| Classification | Telomeric/pericentromeric/interstitial boundaries; length-band boundaries. |
| Gene overlap | Real HFE-region overlap on chr6; empty result off-locus. |
| UPD eligibility | Autosome-only gating; imprinted vs. non-imprinted thresholds at and below the boundary. |
| Overview / F-value | Internal consistency of region-count totals; F-value in the expected low (non-consanguineous) range for HG002; denominator matches `AUTOSOMAL_GENOME_LENGTH`. |
| Real findings | The one real chr9 UPD flag in the HG002 fixture is detected with the correct coordinates and classification. |
| AR-candidate floor | TJP2 (overlapped by a real but ~49 kb ROH) is correctly excluded by the 1 Mb floor; the 1 Mb / 3 Mb floors are confirmed distinct. |
| Relationship table | F values present, ordered by closeness of relationship, parent-child/sibling value = 0.25. |
| Query / pagination | Chromosome and length filtering; stable pagination and totals across pages. |
| API routes | `/api/roh/samples`, `/api/roh/overview`, `/api/roh/regions` on the real HG002 sample; 404 on an unknown sample; other RepeatScope routes remain unaffected. |
