// Typed client for the structural-variant backend (backend/structvar/router.py), served under
// /api/structvar on the same FastAPI app as RepeatScope.

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

export interface Roles {
  proband: string | null
  parents: string[]
  source: 'inferred' | 'sample.json' | 'none'
  basis: string
}

export interface Meta {
  sample: string
  label: string
  reference: string
  caller: string | null
  vcf: string
  status: BuildStatus
  built?: string
  n_svs?: number
  roles?: Roles
  source?: string | null
  samples?: string[]
}

export type SvType = 'DEL' | 'DUP' | 'INS' | 'INV' | 'CNV' | 'BND'

export interface Qc {
  records: number
  skipped: number
  invalid_records: number
  source: string | null
  samples: string[]
  by_type_all: Record<SvType, number>
  by_type_pass: Record<SvType, number>
  filters: { filter: string; count: number }[]
  pass_non_bnd: number
  imprecise: number
  size_bins: ({ bin: string } & Partial<Record<SvType, number>>)[]
  per_sample: { sample: string; role: string; carrier: number; ref: number; nocall: number }[]
  roles: Roles
  de_novo_candidates: number | null
  with_gene: number
  with_population_match: number
}

export interface Resource {
  key: string
  label: string
  installed: boolean
  detail: string
}

export interface Overview {
  sample: string
  qc: Qc
  funnel: { label: string; count: number }[]
  top_genes: { gene: string; count: number }[]
  resources: Resource[]
  trio: boolean
}

export interface GeneHit {
  symbol: string
  gene_id: string
  biotype: string
  relation: 'whole_gene' | 'coding_exon' | 'noncoding_exon' | 'intronic'
  exons: number[]
  n_exons: number
  strand: string
  transcript: string
  model_source: string
  rank: number
}

export interface PopMatch {
  svtype: string
  similarity: number
  n: number
  af: number | null
  ac: number | null
  an: number | null
  popmax: number | null
}

export interface FamilyMember {
  sample: string
  role: 'proband' | 'parent' | 'member'
  gt: string
  copies: number
  partial: boolean
  gq: number | null
  alt: number | null
  ref: number | null
  dp: number | null
  cn: number | null
}

export interface Sv {
  id: number
  vid: string
  chrom: string
  pos: number
  end: number
  svtype: SvType
  size: number | null
  filter: string
  passed: boolean
  qual: number | null
  imprecise: boolean
  symbolic: boolean
  mate: { chrom: string; pos: number } | null
  event: string | null
  focus: { gt: string | null; copies: number | null; gq: number | null; alt: number | null; dp: number | null; cn: number | null }
  inheritance: 'de_novo' | 'inherited' | 'unknown' | null
  parent_min_dp: number | null
  family: FamilyMember[]
  n_genes: number
  genes: GeneHit[]
  impact: string
  population: Record<string, PopMatch>
  pop_af_max: number | null
  pop_n: number
  segdup_frac: number | null
  odd: string[]
  repeat: string | null
  repeat_frac: number | null
}

export interface SvsResponse {
  total_all: number
  total_matching: number
  page: number
  page_size: number
  rows: Sv[]
  funnel: { label: string; count: number }[] | null
  ignored_filters: string[]
  population_installed: boolean
  genes_installed: boolean
  tracks_installed: boolean
  trio: boolean
}

export interface RegionGene {
  symbol: string
  gene_id: string
  biotype: string
  start: number
  end: number
  strand: string
  transcript: string
  model_source: string
  exons: [number, number][]
  cds: [number, number] | null
}

export interface RegionResponse {
  chrom: string
  start: number
  end: number
  genes: RegionGene[]
  n_genes: number
  svs: { id: number; pos: number; end: number; svtype: SvType; size: number | null; copies: number; impact: string; imprecise: boolean }[]
}

export interface GeneInfo {
  symbol: string
  gene_id: string
  biotype: string
  transcript: string
  ensembl_transcript: string
  model_source: string
  n_exons: number
  pli: number | null
  oe_lof: number | null
  inheritance: string[]
  hpo_n: number
  hpo_terms: string[]
  diseases: string[]
}

export interface FilterState {
  passOnly: boolean
  includeBnd: boolean
  types: string[]
  minSize: string
  maxSize: string
  carried: boolean
  zygosity: 'any' | 'het' | 'hom'
  inheritance: 'any' | 'de_novo' | 'inherited' | 'unknown'
  minParentDp: string
  minGq: string
  minAlt: string
  impact: 'any' | 'gene' | 'exon' | 'coding' | 'whole_gene'
  genes: string
  chrom: string
  maxPopAf: string
  novel: boolean
  maxSegdup: string
  noOdd: boolean
  noRepeat: boolean
  preciseOnly: boolean
  sort: 'position' | 'size' | 'impact' | 'af' | 'gq'
  order: 'asc' | 'desc'
}

export const DEFAULT_FILTERS: FilterState = {
  passOnly: true, includeBnd: false, types: [], minSize: '', maxSize: '', carried: true, zygosity: 'any', inheritance: 'any', minParentDp: '',
  minGq: '', minAlt: '', impact: 'any', genes: '', chrom: '', maxPopAf: '', novel: false, maxSegdup: '', noOdd: false, noRepeat: false,
  preciseOnly: false, sort: 'position', order: 'asc',
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

const num = (s: string): number | undefined => (s === '' || Number.isNaN(Number(s)) ? undefined : Number(s))

export const structvarApi = {
  samples: () => getJSON<SamplesResponse>('/api/structvar/samples'),
  // Asking for status starts the index build the first time a call set is requested.
  status: (sample: string) => getJSON<BuildStatus>(`/api/structvar/status${qs({ sample })}`),
  meta: (sample: string) => getJSON<Meta>(`/api/structvar/meta${qs({ sample })}`),
  overview: (sample: string) => getJSON<Overview>(`/api/structvar/overview${qs({ sample })}`),
  gene: (symbol: string) => getJSON<GeneInfo>(`/api/structvar/genes/${encodeURIComponent(symbol)}`),
  region: (sample: string, chrom: string, start: number, end: number) =>
    getJSON<RegionResponse>(`/api/structvar/region${qs({ sample, chrom, start, end })}`),
  svs: (sample: string, f: FilterState, page: number, pageSize: number) =>
    getJSON<SvsResponse>(
      `/api/structvar/svs${qs({
        sample,
        pass_only: f.passOnly,
        include_bnd: f.includeBnd || undefined,
        types: f.types.join(','),
        min_size: num(f.minSize),
        max_size: num(f.maxSize),
        carried: f.carried,
        zygosity: f.zygosity,
        inheritance: f.inheritance,
        min_parent_dp: num(f.minParentDp),
        min_gq: num(f.minGq),
        min_alt: num(f.minAlt),
        impact: f.impact,
        genes: f.genes.replace(/\s+/g, ','),
        chrom: f.chrom,
        max_pop_af: num(f.maxPopAf),
        novel: f.novel || undefined,
        max_segdup: f.maxSegdup === '' ? undefined : Number(f.maxSegdup) / 100,
        no_odd: f.noOdd || undefined,
        no_repeat: f.noRepeat || undefined,
        precise_only: f.preciseOnly || undefined,
        sort: f.sort,
        order: f.order,
        page,
        page_size: pageSize,
      })}`,
    ),
}
