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
}
