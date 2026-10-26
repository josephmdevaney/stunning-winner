// Typed client for the Epigenetics backend (backend/epigenetics/router.py), served under
// /api/epigenetics on the same FastAPI app as the other suites.

export type Status = 'imprinted' | 'loss' | 'gain' | 'borderline' | 'low_coverage' | 'no_data'
export type Category = 'CA' | 'NC' | 'X'
export type DisorderStatus = 'pattern_present' | 'supportive_only' | 'not_observed' | 'insufficient_data'
export type Window = 'core' | 'dmr'

export interface TrackStats {
  n_cpgs: number
  mean_coverage: number
  methylation: number // pooled percent
}

export interface Interval {
  chrom: string
  start: number
  end: number
}

export interface RegionResult {
  id: string
  name: string
  chrom: string
  category: Category // CA clinically associated | NC non-clinical | X polymorphic (not interpreted)
  primary_diagnostic: boolean
  interpretable: boolean
  methylated_allele: 'maternal' | 'paternal' | null
  disorders: string[]
  disorder_notes: string[]
  window: Window
  combined: TrackStats | null
  hap1: TrackStats | null
  hap2: TrackStats | null
  status: Status
  phasing: 'phased' | 'unphased' | 'none'
  pattern: string | null
  allelic_difference: number | null
  notes: string[]
  interval: Interval
  dmr: Interval
  n_cpgs_published: number
  imprint_origin: string
  coordinate_notes: string | null
}

export interface DisorderTrigger {
  dmr: string
  name: string
  expected: 'loss' | 'gain'
  primary: boolean
  observed: Status | 'not_assessed'
  hit: boolean
  methylation: number | null
}

export interface DisorderSummary {
  id: string
  name: string
  locus: string
  status: DisorderStatus
  meaning: string
  mechanisms: string
  triggers: DisorderTrigger[]
}

export interface Mlid {
  scope: string
  assessed: number
  considered: number
  aberrant: { id: string; name: string; chrom: string; status: Status }[]
  chromosomes: string[]
  possible: boolean
}

export interface Thresholds {
  min_cpgs: number
  min_mean_coverage: number
  min_hap_mean_coverage: number
  loss_below: number
  gain_above: number
  imprint_low: number
  imprint_high: number
  allelic_min_diff: number
}

export interface Provenance {
  pb_cpg_tools_version: string | null
  pileup_mode: string | null
  modsites_mode: string | null
  min_coverage: number | null
  tracks: string[]
  genome_build: string
}

export interface ImprintingOverview {
  sample: string
  label: string
  note: string | null
  window: Window
  thresholds: Thresholds
  provenance: Provenance
  panel: { title: string; journal: string; pmcid: string; license: string }
  summary: { regions: number; called: number; by_status: Record<string, number> }
  regions: RegionResult[]
  disorders: DisorderSummary[]
  mlid: Mlid
  caveats: string
}

export interface TrackPoint {
  pos: number
  methylation: number
  coverage: number
}

export interface RegionTrack {
  id: string
  name: string
  chrom: string
  dmr: { start: number; end: number }
  analysis_window: { start: number; end: number; which: Window }
  view: { start: number; end: number }
  tracks: Record<'combined' | 'hap1' | 'hap2', TrackPoint[]>
}

export interface SampleInfo {
  key: string
  label: string
  note: string | null
  tracks: string[]
}

export interface SamplesResponse {
  samples: SampleInfo[]
  default: string
}

export interface ThresholdOverrides {
  loss_below?: number
  gain_above?: number
  allelic_min_diff?: number
}

async function getJSON<T>(path: string): Promise<T> {
  const res = await fetch(path)
  if (!res.ok) {
    let detail = ''
    try {
      detail = (await res.json()).detail ?? ''
    } catch {
      /* not JSON */
    }
    throw new Error(detail || `${path} -> ${res.status}`)
  }
  return res.json() as Promise<T>
}

function query(params: Record<string, string | number | undefined>): string {
  const q = Object.entries(params)
    .filter(([, v]) => v !== undefined && v !== '')
    .map(([k, v]) => `${k}=${encodeURIComponent(String(v))}`)
    .join('&')
  return q ? `?${q}` : ''
}

export const epigeneticsApi = {
  samples: () => getJSON<SamplesResponse>('/api/epigenetics/samples'),
  imprinting: (sample: string, window: Window, overrides: ThresholdOverrides = {}) =>
    getJSON<ImprintingOverview>(`/api/epigenetics/imprinting${query({ sample, window, ...overrides })}`),
  region: (sample: string, id: string, window: Window, flank = 3000) =>
    getJSON<RegionTrack>(`/api/epigenetics/imprinting/region${query({ sample, id, window, flank })}`),
  episignatures: (sample: string) => getJSON<EpisignatureResponse>(`/api/epigenetics/episignatures${query({ sample })}`),
  arraySignatures: (sample: string) => getJSON<ArrayResponse>(`/api/epigenetics/episignatures/array${query({ sample })}`),
}

// ---- Episignatures (EpiSignaLR) -------------------------------------------------------------------------------

export type EpiCall = 'positive' | 'indeterminate' | 'negative' | 'insufficient_data'

export interface EpiWarning {
  id: string
  severity: 'caution' | 'note'
  text: string
}

export interface EpiFeature {
  dmr: string
  region: Interval | null
  weight: number
  z: number
  contribution: number // this DMR's share of the SVM decision value; positive pushes toward "case"
  value: number | null // the genome's methylation fraction (0-1), null when the DMR had no covered CpG
  imputed: boolean
  n_cpgs: number
  n_cpgs_listed: number
  mean_coverage: number
  training_mean: number
  controls: number[]
  cases: number[]
}

export interface EpiModel {
  id: string
  name: string
  gene: string
  probability: number | null
  decision: number
  intercept: number
  call: EpiCall
  threshold: number
  indeterminate_zone: [number, number]
  features: EpiFeature[]
  dmrs: { panel: number; selected: number; with_data: number; imputed_selected: number }
  training: { cases: number; controls: number }
  performance: {
    design: string
    auc: number
    sensitivity: number
    specificity: number
    balanced_accuracy: number
    permutation_p: number
    confusion_matrix: number[][]
  }
  loo: { group: 0 | 1; score: number }[]
  warnings: EpiWarning[]
}

export interface EpisignatureOverview {
  installed: true
  sample: string
  label: string
  source: { repo: string; commit: string; fetched_at: string; source: string; licence: string }
  provenance: Provenance
  models: EpiModel[]
}

export interface EpisignatureNotInstalled {
  installed: false
  setup_command: string
  repo: string
  available: { id: string; name: string; gene: string }[]
  sample: string
  label: string
  models: []
}

export type EpisignatureResponse = EpisignatureOverview | EpisignatureNotInstalled

// ---- Array-derived episignatures (EpiSignaLR's 34 classifiers) ----------------------------------------------

export interface ArrayCondition {
  id: string
  name: string
  genes: string
  score: number // SVM decision value; at or above the threshold is a match
  matches: boolean
  n_cpgs: number
  n_found: number
  cpgs: {
    coords: [string, number][] // chromosome, 1-based position
    control: (number | null)[] // reference median methylation (0-1) of controls
    condition: (number | null)[] // reference median methylation of this condition
    sample: (number | null)[] // this genome, null where the pileup has no row
  }
}

export interface ArrayOverview {
  installed: true
  sample: string
  label: string
  threshold: number
  call: 'control' | 'one' | 'several'
  matches: string[]
  conditions: ArrayCondition[] // highest score first
  cpgs: { total: number; found: number }
  scorer: { scikit_learn: string; reference_scikit_learn: string }
  warnings: EpiWarning[]
}

export interface ArrayNotInstalled {
  installed: false
  setup_command: string
  sample: string
  label: string
  conditions: []
}

export type ArrayResponse = ArrayOverview | ArrayNotInstalled
