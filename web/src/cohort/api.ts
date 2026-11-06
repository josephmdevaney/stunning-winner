// Typed client for /api/cohort (see backend/repeatscope/api/main.py): one aggregated view of every
// sample in a Carrier or PureTarget V2 run.
export type CohortMode = 'carrier' | 'puretarget_v2'

export interface CohortFinding {
  gene: string
  locus_id: string | null
  disease: string
  status: string
  genotype: string | null
  zygosity: string | null
  source: string
}

export interface CohortSample {
  sample: string
  run: string
  sex: string | null
  loci_genotyped: number
  loci_total: number
  no_calls: number
  mean_depth: number | null
  low_depth_loci: number
  qc: 'Pass' | 'Low depth' | 'No data'
  reportable: CohortFinding[]
}

export interface CohortResponse {
  mode: CohortMode
  label: string
  min_locus_depth: number
  stats: {
    samples: number
    qc_pass: number
    qc_fail: number
    mean_depth: number | null
    median_depth: number | null
    samples_with_reportable: number
    reportable_findings: number
  }
  depth_histogram: { label: string; lo: number; count: number }[]
  sex: { label: string; count: number }[]
  runs: { run: string; count: number }[]
  genes: { gene: string; disease: string; samples: number }[]
  samples: CohortSample[]
}

export const cohortApi = {
  async get(mode: CohortMode): Promise<CohortResponse> {
    const r = await fetch(`/api/cohort?mode=${mode}`)
    if (!r.ok) throw new Error(`${r.status} ${r.statusText}`)
    return r.json()
  },
}
