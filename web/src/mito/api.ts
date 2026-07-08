// Typed client for the MitoLR backend (backend/mitolr/router.py), served under
// /api/mito on the same FastAPI app as RepeatScope.

export interface MitoMeta {
  sample: string
  instrument_line: string
  reference: string
  run_date: string
  caller: string
}

export interface Variant {
  id: string
  position: string
  change: string
  locus: string
  pos: number
  stop: number | null
  type: string
  consequence: string
  disease: string | null
  inheritance: string
  evidence: string
  isDeletion: boolean
  qual: number
  filter: string
  gnomadFreq: string
  origin: 'pathogenic' | 'novel'
  heteroplasmy: number
  status: string
  statusColor: string
  note: string
  depth: number
}

export interface Haplogroup {
  name: string
  confidence: number
  region: string
  markers: string[]
}

export interface OverviewResponse extends MitoMeta {
  stats: {
    variants_called: number
    pathogenic_findings: number
    haplogroup: string
    qc_status: string
  }
  variants: Variant[]
  haplogroup: Haplogroup
  consequence_options: string[]
}

export interface GeneRegion {
  name: string
  start: number
  stop: number
  type: string
}

export interface PileupRead {
  kind: 'mut' | 'wt'
  x: number
  len: number
}

export interface VariantDetailResponse extends MitoMeta {
  variant: Variant
  alt_base: string
  pileup_reads: PileupRead[]
  pileup_source: string
  gene_map: GeneRegion[]
  genome_length: number
  reference_ranges: { label: string; range: string }[]
}

export interface QCResponse extends MitoMeta {
  overall_status: string
  stats: { label: string; value: string; sub: string; dotColor: string }[]
  read_length_histogram: { kb: number; count: number }[]
  coverage_track: { pos: number; depth: number }[]
  genome_length: number
  regions: { region: string; depth: number; strandBalance: string; flag: string }[]
  min_required: number
  instrument: { label: string; value: string }[]
  coverage_y_axis_max_default: number
  unavailable_metrics: { label: string; reason: string }[]
}

async function getJSON<T>(path: string): Promise<T> {
  const res = await fetch(path)
  if (!res.ok) throw new Error(`${path} -> ${res.status}`)
  return res.json() as Promise<T>
}

export const mitoApi = {
  meta: () => getJSON<MitoMeta>('/api/mito/meta'),
  overview: () => getJSON<OverviewResponse>('/api/mito/overview'),
  variant: (id: string) => getJSON<VariantDetailResponse>(`/api/mito/variant/${encodeURIComponent(id)}`),
  qc: () => getJSON<QCResponse>('/api/mito/qc'),
}
