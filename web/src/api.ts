// Typed client for the RepeatScope FastAPI backend (backend/repeatscope/api/main.py).

export interface Meta {
  sample: string
  source_vcf: string
  caller: string
  sex: string | null
}

export interface LocusSummary {
  id: string
  gene: string
  disease: string
  region: string
  motif: string
  a1: number | null
  a2: number | null
  zygosity: string | null
  depth: number | null
  status: string
  status_color: string
  genotyped: boolean
  note: string
  clickable: boolean
}

export interface OverviewResponse {
  sample: string
  caller: string
  stats: {
    loci_analyzed: number
    loci_genotyped: number
    pathogenic_findings: number
    mean_read_depth: number | null
    qc_status: string
  }
  loci: LocusSummary[]
}

export interface LocusDetail {
  id: string
  gene: string
  disease: string
  chrom: string
  start: number
  stop: number
  motif: string
  inheritance: string
  evidence: string
  benign_min: number | null
  benign_max: number | null
  intermediate_min: number | null
  intermediate_max: number | null
  pathogenic_min: number | null
  pathogenic_max: number | null
  genotyped: boolean
  trid: string | null
  phased: boolean
  a1: number | null
  a2: number | null
  a1_bp: number | null
  a2_bp: number | null
  a1_range_bp: [number, number] | null
  a2_range_bp: [number, number] | null
  a1_depth: number | null
  a2_depth: number | null
  depth: number | null
  a1_purity: number | null
  a2_purity: number | null
  a1_methylation: number | null
  a2_methylation: number | null
  a1_motif: string | null
  a2_motif: string | null
  zygosity: string | null
  status: string
  status_color: string
  note: string
  region: string
  pileup_reads: { value: number; cluster: 1 | 2 }[]
  pileup_source: string
  methylation_source: string
  caller: string
}

export interface UnavailableMetric {
  label: string
  reason: string
}

export interface SampleQCResponse {
  sample: string
  overall_status: string
  loci_passing: number
  loci_total: number
  mean_locus_depth: number | null
  rows: { gene: string; id: string; depth: number | null; min_required: number; flag: string }[]
  unavailable_metrics: UnavailableMetric[]
}

export interface RunQCResponse {
  sample: string
  repeat_regions_total: number
  repeat_regions_genotyped: number
  repeat_region_mean_depth: number | null
  mean_depth_by_chromosome: {
    chrom: string
    mean_repeat_region_depth: number | null
    n_regions: number
    n_genotyped: number
  }[]
  note: string
  unavailable_metrics: UnavailableMetric[]
}

export interface CatalogRow {
  chrom: string
  start: number
  stop: number
  length: number
  motifs: string
  struc: string
  genotyped: boolean
  a1: number | null
  a2: number | null
  depth: number | null
}

export interface CatalogResponse {
  total_all: number
  total_matching: number
  unique_motifs: number
  chromosomes: string[]
  page: number
  page_size: number
  rows: CatalogRow[]
}

async function getJSON<T>(path: string): Promise<T> {
  const res = await fetch(path)
  if (!res.ok) throw new Error(`${path} -> ${res.status}`)
  return res.json() as Promise<T>
}

interface ModeSample {
  mode?: string
  sample?: string
}

function withModeSample(path: string, ms: ModeSample, extra?: URLSearchParams): string {
  const q = extra ?? new URLSearchParams()
  if (ms.mode) q.set('mode', ms.mode)
  if (ms.sample) q.set('sample', ms.sample)
  const qs = q.toString()
  return qs ? `${path}?${qs}` : path
}

export interface ModesResponse {
  modes: { key: string; label: string; sample_count: number }[]
  default: string
}

export interface SamplesResponse {
  samples: string[]
  default: string
  mode: string
}

export const api = {
  modes: () => getJSON<ModesResponse>('/api/modes'),
  samples: (mode?: string) => getJSON<SamplesResponse>(withModeSample('/api/samples', { mode })),
  meta: (mode?: string, sample?: string) => getJSON<Meta>(withModeSample('/api/meta', { mode, sample })),
  overview: (mode?: string, sample?: string) =>
    getJSON<OverviewResponse>(withModeSample('/api/overview', { mode, sample })),
  locus: (id: string, mode?: string, sample?: string) =>
    getJSON<LocusDetail>(withModeSample(`/api/locus/${encodeURIComponent(id)}`, { mode, sample })),
  sampleQC: (mode?: string, sample?: string) =>
    getJSON<SampleQCResponse>(withModeSample('/api/qc/sample', { mode, sample })),
  runQC: (mode?: string, sample?: string) =>
    getJSON<RunQCResponse>(withModeSample('/api/qc/run', { mode, sample })),
  catalog: (params: {
    mode?: string; sample?: string; chrom?: string; motif?: string; minLen?: number; page?: number; pageSize?: number
  }) => {
    const q = new URLSearchParams()
    if (params.chrom) q.set('chrom', params.chrom)
    if (params.motif) q.set('motif', params.motif)
    if (params.minLen) q.set('min_len', String(params.minLen))
    q.set('page', String(params.page ?? 0))
    q.set('page_size', String(params.pageSize ?? 50))
    return getJSON<CatalogResponse>(withModeSample('/api/catalog', { mode: params.mode, sample: params.sample }, q))
  },
}

