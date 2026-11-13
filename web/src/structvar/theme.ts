// The structural-variant suite reuses RepeatScope's base tokens.
export { colors, mono } from '../theme'

export const fmtInt = (v: number | null | undefined): string => (v == null ? '—' : v.toLocaleString())
export const fmtPct = (v: number | null | undefined, digits = 0): string => (v == null ? '—' : `${(v * 100).toFixed(digits)}%`)

export function fmtSize(bp: number | null | undefined): string {
  if (bp == null) return '—'
  if (bp >= 1_000_000) return `${(bp / 1_000_000).toFixed(bp >= 10_000_000 ? 0 : 1)} Mb`
  if (bp >= 1_000) return `${(bp / 1_000).toFixed(bp >= 10_000 ? 0 : 1)} kb`
  return `${bp} bp`
}

export function fmtAf(af: number | null | undefined): string {
  if (af == null) return '—'
  if (af === 0) return '0'
  return af >= 0.01 ? af.toFixed(3) : af.toExponential(1)
}

// SV types get their own hues, separate from the status palette (red/green mean pathogenic/normal elsewhere).
export const TYPE_COLOR: Record<string, string> = {
  DEL: 'oklch(0.55 0.17 300)',
  DUP: 'oklch(0.6 0.14 220)',
  INS: 'oklch(0.6 0.12 175)',
  INV: 'oklch(0.66 0.14 75)',
  CNV: 'oklch(0.55 0.1 250)',
  BND: 'oklch(0.6 0.02 260)',
}

export const IMPACT_LABEL: Record<string, string> = {
  whole_gene: 'Whole gene',
  coding_exon: 'Coding exon',
  noncoding_exon: 'Non-coding exon / UTR',
  intronic: 'Intronic',
  none: 'No gene',
}
