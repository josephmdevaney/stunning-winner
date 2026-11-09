// Typed client for the small-variant backend (backend/smallvar/router.py), served under
// /api/smallvar on the same FastAPI app as RepeatScope.

export interface SampleInfo {
  key: string
  label: string
}

export interface SamplesResponse {
  samples: SampleInfo[]
  default: string
}

export type BuildState = 'ready' | 'building' | 'stale' | 'error'

export interface BuildStatus {
  state: BuildState
  progress?: { phase: string; records: number; seconds: number }
  error?: string
}

export interface SexCall {
  sex: string | null
  source: string | null
}

export interface Meta {
  sample: string
  label: string
  reference: string
  caller: string | null
  vcf: string
  status: BuildStatus
  built?: string | null
  n_variants?: number | null
  sex?: SexCall | null
}

export interface Qc {
  records: number
  calls_passing: number
  calls_filtered: number
  snv: number
  mnv: number
  insertions: number
  deletions: number
  complex: number
  het: number
  hom: number
  het_hom_ratio: number | null
  ti_tv: number | null
  phased_het_fraction: number | null
  indel_sizes: { bin: string; count: number }[]
  per_chrom: { chrom: string; count: number }[]
  symbolic_or_star_alleles: number
  sample_name: string | null
  sample_columns: number
}

export interface Resource {
  key: string
  label: string
  installed: boolean
  detail: string
}

export interface EvidenceRow {
  criteria: string
  evidence: string
  available: boolean
  note: string
}

export interface Overview {
  sample: string
  qc: Qc
  funnel: { label: string; count: number }[]
  clinvar: { plp: number; conflicting: number } | null
  clinvar_plp_rare: number | null
  top_plp_genes: { gene: string; count: number }[]
  gnomad_coverage: { found: number; of: number; common_over_1pct: number } | null
  resources: Resource[]
  evidence: EvidenceRow[]
  sex: SexCall
}

export interface ClinvarInfo {
  classification: string
  group: 'plp' | 'conflicting'
  stars: number
  review: string
  variation_id: string
  name: string
  conditions: string
  last_evaluated: string
  n_submitters: number
  genes: string
}

export interface FreqInfo {
  af: number
  ac: number
  hom: number
}

export interface Variant {
  id: number
  chrom: string
  pos: number
  ref: string
  alt: string
  kind: 'snv' | 'mnv' | 'ins' | 'del' | 'complex'
  genes: string[]
  genotype: { gt: string; copies: number; zygosity: 'het' | 'hom' | 'hemizygous'; phased: boolean; phase_set: string | null; haplotypes: string }
  quality: { qual: number | null; filter: string; passed: boolean; dp: number | null; gq: number | null; vaf: number | null }
  clinvar: ClinvarInfo | null
  gnomad: FreqInfo | null
  hprc: FreqInfo | null
}

export interface VariantsResponse {
  total_all: number
  total_matching: number
  page: number
  page_size: number
  rows: Variant[]
  funnel: { label: string; count: number }[] | null
  ignored_filters: string[]
  gnomad_installed: boolean
  hprc_installed: boolean
  clinvar_installed: boolean
}

export interface GenePreset {
  key: string
  label: string
  genes: string[]
  note: string
}

export interface FilterState {
  passOnly: boolean
  kinds: string[]
  zygosity: 'any' | 'het' | 'hom'
  chrom: string
  inGene: boolean
  genes: string
  clinvar: 'any' | 'plp' | 'plp_conflicting' | 'has_match' | 'no_match'
  minStars: number
  maxAf: string
  afMissing: 'include' | 'exclude'
  maxHprcAf: string
  minGq: string
  sort: 'position' | 'af' | 'stars' | 'gq'
  order: 'asc' | 'desc'
}

export const DEFAULT_FILTERS: FilterState = {
  passOnly: true, kinds: [], zygosity: 'any', chrom: '', inGene: false, genes: '', clinvar: 'any', minStars: 0, maxAf: '',
  afMissing: 'include', maxHprcAf: '', minGq: '', sort: 'position', order: 'asc',
}

async function getJSON<T>(path: string): Promise<T> {
  const res = await fetch(path)
  if (!res.ok) throw new Error(`${path} -> ${res.status}`)
  return res.json() as Promise<T>
}

function qs(params: Record<string, string | number | boolean | undefined>): string {
  const parts = Object.entries(params)
    .filter(([, v]) => v !== undefined && v !== '')
    .map(([k, v]) => `${k}=${encodeURIComponent(String(v))}`)
  return parts.length ? `?${parts.join('&')}` : ''
}

export const smallvarApi = {
  samples: () => getJSON<SamplesResponse>('/api/smallvar/samples'),
  // Asking for status starts the index build the first time a sample is requested.
  status: (sample: string) => getJSON<BuildStatus>(`/api/smallvar/status${qs({ sample })}`),
  meta: (sample: string) => getJSON<Meta>(`/api/smallvar/meta${qs({ sample })}`),
  overview: (sample: string) => getJSON<Overview>(`/api/smallvar/overview${qs({ sample })}`),
  genePresets: () => getJSON<{ presets: GenePreset[] }>('/api/smallvar/gene-presets'),
  variants: (sample: string, f: FilterState, page: number, pageSize: number, explain = false) =>
    getJSON<VariantsResponse>(
      `/api/smallvar/variants${qs({
        sample,
        pass_only: f.passOnly,
        kinds: f.kinds.join(','),
        zygosity: f.zygosity,
        chrom: f.chrom,
        in_gene: f.inGene || undefined,
        genes: f.genes.replace(/\s+/g, ','),
        clinvar: f.clinvar,
        min_stars: f.minStars || undefined,
        max_af: f.maxAf === '' ? undefined : Number(f.maxAf),
        af_missing: f.afMissing,
        max_hprc_af: f.maxHprcAf === '' ? undefined : Number(f.maxHprcAf),
        min_gq: f.minGq === '' ? undefined : Number(f.minGq),
        sort: f.sort,
        order: f.order,
        page,
        page_size: pageSize,
        explain: explain || undefined,
      })}`,
    ),
}
