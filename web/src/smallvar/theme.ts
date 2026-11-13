// The small-variant suite reuses RepeatScope's base tokens.
import { colors } from '../theme'

export { colors, mono } from '../theme'

export const fmtInt = (v: number | null | undefined): string => (v == null ? '—' : v.toLocaleString())
export const fmtPct = (v: number | null | undefined, digits = 1): string => (v == null ? '—' : `${(v * 100).toFixed(digits)}%`)

// Allele frequencies span six orders of magnitude, so small ones are shown in scientific form.
export function fmtAf(af: number | null | undefined): string {
  if (af == null) return '—'
  if (af === 0) return '0'
  if (af >= 0.01) return af.toFixed(3)
  return af.toExponential(1)
}

// ClinVar's own assertion is coloured with the status palette; a frequency flag uses the amber.
export function clinvarColor(group: 'plp' | 'conflicting'): string {
  return group === 'plp' ? colors.statusPathogenic : colors.statusUncertain
}
export function clinvarBg(group: 'plp' | 'conflicting'): string {
  return group === 'plp' ? colors.interpBg.pathogenic : colors.interpBg.uncertain
}
