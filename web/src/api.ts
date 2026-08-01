// Typed client for the RepeatScope FastAPI backend (backend/repeatscope/api/main.py).

// LocusSummary.id prefix for gene findings (non-repeat variants, e.g. Paraphase/
// sawfish calls) -- see backend/repeatscope/api/main.py. Distinguishes them from
// repeat-locus ids so callers know to route to GeneFindingScreen instead of
// LocusDetailScreen (see OverviewScreen and TopBar's search).
export const NON_REPEAT_ID_PREFIX = '__non_repeat__'

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

export interface GeneFindingDetail {
  gene: string
  disease: string
  region: string
  tool: string
  mechanism: string
  status: string
  status_color: string
  zygosity: string | null
  note: string
  variants: string[]
  raw_fields: Record<string, string>
  genotyped: boolean
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
  scoped_to_panel: boolean
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
  // Real gene(s) this region falls inside (rare but real for overlapping/
  // antisense transcripts), or -- if intergenic -- the nearest gene and its
  // distance/side. genes is empty and nearest_gene is null only if no gene
  // reference data was available server-side, not a fabricated "no gene" claim.
  genes: string[]
  nearest_gene: string | null
  nearest_distance_bp: number | null
  nearest_direction: 'upstream' | 'downstream' | null
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

export interface RohGeneHit {
  symbol: string
  inheritance: string
  disease: string
  omim: string
  start: number
  end: number
}

export interface RohRegion {
  sample: string
  chrom: string
  start: number
  end: number
  length_bp: number
  length_mb: number
  n_markers: number
  quality: number
  position: string | null
  length_band: string
  imprinted_chrom: boolean
  genes: RohGeneHit[]
  upd_flag: string | null
  autosomal: boolean
  f_value_eligible: boolean
}

export interface RohByChromosome {
  chrom: string
  chrom_length_bp: number | null
  centromere_mid_bp: number | null
  n_regions: number
  total_length_bp: number
  imprinted: boolean
}

export interface RohArCandidateGene {
  symbol: string
  disease: string
  omim: string
  chrom: string
  roh_length_mb: number
}

export interface RohRelationshipFValue {
  relationship: string
  f_value: number
  f_value_pct: number
}

export interface RohOverviewResponse {
  sample: string
  genome_build: string
  caller: string
  by_chromosome: RohByChromosome[]
  stats: {
    n_regions_total: number
    n_regions_autosomal: number
    n_regions_chrX: number
    n_regions_chrY: number
    n_regions_chrM: number
    n_regions_f_eligible: number
    total_autosomal_roh_bp: number
    total_f_eligible_roh_bp: number
    f_value_min_length_mb: number
    ar_candidate_min_length_mb: number
    autosomal_genome_length_bp: number
    f_value: number
    f_value_pct: number
    n_upd_flags: number
    n_ar_candidate_genes: number
  }
  upd_flags: RohRegion[]
  ar_candidate_genes: RohArCandidateGene[]
  relationship_f_values: RohRelationshipFValue[]
  non_autosomal_note: string
  upd_methodology_note: string
}

export interface RohRegionsResponse {
  total_all: number
  total_matching: number
  chromosomes: string[]
  page: number
  page_size: number
  rows: RohRegion[]
}

export interface RohSamplesResponse {
  samples: string[]
  default: string
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
  geneFinding: (gene: string, mode?: string, sample?: string) => {
    const q = new URLSearchParams()
    q.set('gene', gene)
    return getJSON<GeneFindingDetail>(withModeSample('/api/gene-finding', { mode, sample }, q))
  },
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

// ROH is self-contained: keyed by `sample` only, independent of RepeatScope's
// assay `mode`/panel system (bcftools roh doesn't depend on a repeat panel).
export const rohApi = {
  samples: () => getJSON<RohSamplesResponse>('/api/roh/samples'),
  overview: (sample?: string) =>
    getJSON<RohOverviewResponse>(sample ? `/api/roh/overview?sample=${encodeURIComponent(sample)}` : '/api/roh/overview'),
  regions: (params: { sample?: string; chrom?: string; minLenMb?: number; page?: number; pageSize?: number }) => {
    const q = new URLSearchParams()
    if (params.sample) q.set('sample', params.sample)
    if (params.chrom) q.set('chrom', params.chrom)
    if (params.minLenMb) q.set('min_len_mb', String(params.minLenMb))
    q.set('page', String(params.page ?? 0))
    q.set('page_size', String(params.pageSize ?? 100))
    return getJSON<RohRegionsResponse>(`/api/roh/regions?${q.toString()}`)
  },
}

