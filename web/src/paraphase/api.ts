// Typed client for the Paraphase Annotation backend (backend/paraphase_annotate/router.py),
// served under /api/paraphase on the same FastAPI app as RepeatScope.

export type RunState = 'not_run' | 'blocked' | 'running' | 'done' | 'failed'

export interface StatusResponse {
  state: RunState
  reason?: string
  log_tail?: string
}

export interface RegionSummary {
  region: string
  region_median_depth: number | null
  gene: string | null
  total_haplotypes: number | null
  n_haplotypes_with_pathogenic: number | null
  n_haplotypes_without_pathogenic: number | null
  // PacBio's own bar for clinical significance in a recessive-disease model:
  // zero haplotypes free of a pathogenic variant is the long-read analog of
  // biallelic mutations (see the Paraphase annotation module's own README).
  zero_functional_haplotypes: boolean
  info: string | null
}

export interface VariantRecord {
  variant: string
  transcript: string | null
  consequence_terms: string[]
  impact: string | null
  hgvsc: string | null
  hgvsp: string | null
  clinvar_id: string | null
  // Real, freshly-verified ClinVar data (see backend/paraphase_annotate/enrich.py) --
  // classification is ClinVar's own germline_classification description
  // (e.g. "Pathogenic", "Conflicting classifications of pathogenicity"), not
  // assumed from the fact that this variant appears in this list at all.
  clinvar_classification: string | null
  clinvar_review_status: string | null
  clinvar_condition: string | null
  // True only for ACMG's own definitive terms (Pathogenic / Likely pathogenic
  // / Pathogenic/Likely pathogenic) when ClinVar-backed, or a VEP HIGH-impact
  // call when there's no ClinVar entry at all -- never true for Conflicting,
  // Uncertain significance, or any other non-definitive classification.
  is_definitive_pathogenic: boolean
}

export interface HaplotypeDetail {
  region: string
  gene: string
  transcript: string
  transcript_start: number | null
  transcript_end: number | null
  haplotype: string
  haplotype_boundary: string
  // fully_contained | partial | truncated | no_overlap -- a truncated
  // haplotype is a pseudogene copy, not a functional one (see
  // annotation/README.md's own definitions).
  boundary_check: string
  known_nonfunctional: boolean
  pathogenic_variants: VariantRecord[]
}

// How an annotation run was set up (backend/paraphase_annotate/runner.py's
// run_meta.json): which interpret_paraphase.py release, whether the JSON had to
// be adapted to the layout that release reads, and which Paraphase made the data.
export interface RunMeta {
  sample: string
  started_at: string
  finished_at?: string
  script_tag: string | null
  script_fetched_at: string | null
  python: string | null
  input_schema: string | null
  script_schema: string
  adapted_json: boolean
  data_paraphase_version: string | null
  regions_requested: string[]
  regions_missing_from_json: string[]
}

// What setup_paraphase_annotate.py vendored (paraphase_annotate/vendor/VERSION.json).
export interface ScriptInfo {
  paraphase_tag: string | null
  fetched_at: string | null
  source: string | null
}

export interface OverviewResponse extends StatusResponse {
  sample: string
  regions: RegionSummary[]
  run?: RunMeta | null
  script?: ScriptInfo | null
}

// --- Paraphase's own JSON (backend/paraphase_annotate/regions.py) -------------

export type FlagKind = 'no_data' | 'fusion' | 'sv' | 'phasing_failed' | 'cn_unresolved'

export interface JsonFlag {
  kind: FlagKind
  text: string
}

export interface JsonRegionRow {
  region: string
  genes: string[]
  medical: boolean
  total_cn: number | null
  gene_cn: number | null
  n_haplotypes: number
  two_copy_haplotypes: number
  depth_median: number | null
  depth_p80: number | null
  genome_depth: number | null
  depth_vs_genome: number | null
  headline: string | null
  flags: JsonFlag[]
  phase_region: string | null
}

export interface JsonTable {
  sample: string
  source_schema: string | null
  paraphase_version: string | null
  paraphase_command: string | null
  sample_sex: string | null
  n_regions: number
  regions: JsonRegionRow[]
}

export interface JsonHaplotype {
  name: string
  group: string | null
  kind: string | null
  boundary: number[] | null
  boundary_gene2: (number | null)[] | null
  truncated: string[]
  haplogroup: string | null
  two_copy: boolean
  unique_reads: number | null
  n_variants: number
  variants: string[]
}

export interface JsonRegionDetail {
  sample: string
  summary: JsonRegionRow
  haplotypes: JsonHaplotype[]
  alleles: string[][] | null
  fusions: { name: string; type: string | null; breakpoint: unknown }[]
  region_fields: Record<string, string | number | boolean | null>
}

// --- PacBio tool versions (backend/pacbio_tools/) ------------------------------

export type VersionGap = 'none' | 'patch' | 'minor' | 'major' | 'ahead' | null

export interface ToolObservation {
  label: string
  raw: string | null
  version: string | null
  count: number
  examples: string[]
  gap: VersionGap
}

export interface ToolRelease {
  tag: string
  version: string | null
  published: string | null
  name: string
  url: string | null
}

export type ToolStatus = 'current' | 'behind' | 'unknown' | 'not_recorded' | 'not_used'

export interface ToolVersion {
  key: string
  name: string
  repo: string
  url: string
  role: string
  used: boolean
  suite_use: string
  latest: ToolRelease | null
  observed: ToolObservation[]
  status: ToolStatus
  worst_gap: VersionGap
}

export interface ToolVersions {
  source: 'live' | 'cache' | 'snapshot'
  fetched_at: string | null
  errors: Record<string, string>
  summary: { used: number; current: number; behind: number; not_recorded: number }
  tools: ToolVersion[]
}

export interface RegionResponse extends StatusResponse {
  sample: string
  region: string | null
  haplotypes: HaplotypeDetail[]
}

export interface SamplesResponse {
  samples: string[]
}

async function getJSON<T>(path: string): Promise<T> {
  const res = await fetch(path)
  if (!res.ok) throw new Error(`${path} -> ${res.status}`)
  return res.json() as Promise<T>
}

function withSample(path: string, sample?: string): string {
  return sample ? `${path}${path.includes('?') ? '&' : '?'}sample=${encodeURIComponent(sample)}` : path
}

export const paraphaseApi = {
  jsonTable: (sample: string) => getJSON<JsonTable>(withSample('/api/paraphase/json', sample)),
  jsonRegion: (region: string, sample: string) =>
    getJSON<JsonRegionDetail>(withSample(`/api/paraphase/json/region?region=${encodeURIComponent(region)}`, sample)),
  toolVersions: (refresh = false) => getJSON<ToolVersions>(`/api/tools/versions${refresh ? '?refresh=true' : ''}`),
  samples: () => getJSON<SamplesResponse>('/api/paraphase/samples'),
  status: (sample?: string) => getJSON<StatusResponse>(withSample('/api/paraphase/status', sample)),
  run: async (sample: string) => {
    const res = await fetch(withSample('/api/paraphase/run', sample), { method: 'POST' })
    if (!res.ok) throw new Error(`run -> ${res.status}`)
    return res.json() as Promise<StatusResponse>
  },
  overview: (sample?: string) => getJSON<OverviewResponse>(withSample('/api/paraphase/overview', sample)),
  region: (region: string, sample?: string) =>
    getJSON<RegionResponse>(withSample(`/api/paraphase/region?region=${encodeURIComponent(region)}`, sample)),
}
