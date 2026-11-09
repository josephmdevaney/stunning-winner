// Typed client for the Kivvi backend (backend/kivvi/router.py), served under
// /api/kivvi on the same FastAPI app as RepeatScope.

export type Target = 'd4z4' | 'kiv2'
export type CallStatus = 'pathogenic' | 'review' | 'normal' | 'indeterminate'

export interface SampleInfo {
  key: string
  label: string
}

export interface SamplesResponse {
  samples: SampleInfo[]
  default: string
}

export interface KivviMeta {
  sample: string
  reference: string
  kivvi_version: string | null
  targets: Record<Target, boolean>
  // Set when the input was a bamlet, so genome-wide depth is unavailable.
  input_note: string | null
}

export interface AlleleCall {
  status: CallStatus
  label: string
  reason: string
}

export interface D4Z4Allele {
  name: string
  units: string[]
  size: number
  size_is_lower_bound: boolean
  size_label: string
  chromosome: string
  upstream_group: string | null
  distal_haplotype: string
  allele_type: string
  methylation: number | null
  ending_in_qal: boolean
  per_unit_methylation: number[]
  interpretation: AlleleCall
}

export interface D4Z4Summary {
  status: CallStatus
  label: string
  complete_alleles: number
  partial_alleles: number
  notes: string[]
}

// Where an allele's size sits: contracted (1-10 units), the 11-20 range the Kivvi preprint calls commonly
// indicated for FSHD2, above that, or a partial allele whose lower bound does not settle it.
export type SizeRange = 'contracted' | 'fshd2_range' | 'above' | 'unresolved'

export interface Fshd2Allele {
  index: number
  size: number
  size_label: string
  size_is_lower_bound: boolean
  size_range: SizeRange
  chromosome: string
  upstream_group: string | null
  distal_haplotype: string
  // chr4, qAIntactPolyA, and not the 4A166 haplotype: an allele that can express DUX4.
  permissive: boolean
  methylation: number | null
}

// The facts that bear on FSHD2, side by side. Kivvi defines no FSHD2 rule or cut-off, so this carries no verdict.
export interface Fshd2Evidence {
  alleles: Fshd2Allele[]
  sample_methylation: number | null
  summary: {
    permissive_alleles: number
    permissive_contracted: number
    permissive_fshd2_range: number
    alleles_with_methylation: number
    lowest_methylation: number | null
    highest_methylation: number | null
  }
}

export interface Kiv2Allele {
  index: number
  name: string
  units: string[]
  copy_number: number
  supporting_reads: number
}

export interface VariantUnit {
  unit: string
  where: string
  depth: number
  reads: number
}

export interface VariantRecord {
  chrom: string
  pos: number
  ref: string
  alt: string
  filter: string
  genotypes: string[]
  on_assembled_allele: boolean
  units: VariantUnit[]
  n_units: number
  min_vaf: number | null
}

export interface Variants {
  allele_columns: string[]
  // One per allele column, e.g. 'A1 · 28 · chr4' (D4Z4) or 'A1 · 14' (KIV2).
  allele_labels: string[]
  records: VariantRecord[]
  n_total: number
  n_on_assembled: number
}

interface TargetBase {
  sample: string
  target: Target
  label: string
  variants: Variants
  svg_available: boolean
}

export interface D4Z4Detail extends TargetBase {
  alleles: D4Z4Allele[]
  summary: D4Z4Summary
  fshd2: Fshd2Evidence
  overall_methylation: number | null
  // Kivvi's d4z4_arrays_in_cis calls; empty when no in-cis duplication was found.
  in_cis_arrays: unknown[]
  qc: { median_read_length?: number; per_allele_depth?: number }
}

export interface Kiv2Detail extends TargetBase {
  alleles: Kiv2Allele[]
  partial_alleles: { name: string; units: string[]; lower_bound: number }[]
  shortest_copy_number: number | null
  allele_cn: string
  notes: string[]
}

export interface OverviewResponse extends KivviMeta {
  d4z4: {
    summary: D4Z4Summary
    overall_methylation: number | null
    alleles: {
      size_label: string
      chromosome: string
      upstream_group: string | null
      distal_haplotype: string
      methylation: number | null
      interpretation: AlleleCall
    }[]
  } | null
  kiv2: {
    allele_cn: string
    shortest_copy_number: number | null
    alleles: { index: number; copy_number: number }[]
    notes: string[]
  } | null
}

export interface TargetQC {
  target: Target
  label: string
  repeat_depth: number | null
  genome_depth: number | null
  median_read_length: number | null
  per_allele_depth: number | null
  reads_linking_units: number
  alleles_resolved: number
  alleles_expected: number
  partial_alleles: number
  variants_total: number
  variants_on_assembled: number
}

export interface QCResponse extends KivviMeta {
  target_qc: TargetQC[]
}

async function getJSON<T>(path: string): Promise<T> {
  const res = await fetch(path)
  if (!res.ok) throw new Error(`${path} -> ${res.status}`)
  return res.json() as Promise<T>
}

function withSample(path: string, sample?: string): string {
  return sample ? `${path}${path.includes('?') ? '&' : '?'}sample=${encodeURIComponent(sample)}` : path
}

export const kivviApi = {
  samples: () => getJSON<SamplesResponse>('/api/kivvi/samples'),
  meta: (sample?: string) => getJSON<KivviMeta>(withSample('/api/kivvi/meta', sample)),
  overview: (sample?: string) => getJSON<OverviewResponse>(withSample('/api/kivvi/overview', sample)),
  d4z4: (sample?: string) => getJSON<D4Z4Detail>(withSample('/api/kivvi/target/d4z4', sample)),
  kiv2: (sample?: string) => getJSON<Kiv2Detail>(withSample('/api/kivvi/target/kiv2', sample)),
  qc: (sample?: string) => getJSON<QCResponse>(withSample('/api/kivvi/qc', sample)),
  svgUrl: (target: Target, sample: string) => withSample(`/api/kivvi/target/${target}/svg`, sample),
}
