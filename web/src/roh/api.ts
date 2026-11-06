// Typed client for the ROH backend routes (/api/roh/*, see backend/repeatscope/api/main.py).
// ROH is self-contained: keyed by `sample` only, independent of RepeatScope's assay
// `mode`/panel system (bcftools roh doesn't depend on a repeat panel).

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
  hemizygous: boolean
}

export interface RohByChromosome {
  chrom: string
  chrom_length_bp: number | null
  centromere_mid_bp: number | null
  n_regions: number
  total_length_bp: number
  imprinted: boolean
  hemizygous: boolean
}

export interface RohSexInference {
  sex: 'M' | 'F' | null
  method: string | null
  reason?: string
  autosomal_mean_depth?: number
  chrX_mean_depth?: number
  chrY_mean_depth?: number
  chrX_autosomal_ratio?: number
  chrY_autosomal_ratio?: number
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
  sex_inference: RohSexInference
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
