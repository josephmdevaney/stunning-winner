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

export interface OverviewResponse extends StatusResponse {
  sample: string
  regions: RegionSummary[]
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
