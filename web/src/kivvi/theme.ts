// Kivvi reuses RepeatScope's base tokens and maps its call statuses onto the
// shared status palette.
import { colors } from '../theme'
import type { CallStatus } from './api'

export { colors, mono } from '../theme'

export function callColor(status: CallStatus): string {
  if (status === 'pathogenic') return colors.statusPathogenic
  if (status === 'review') return colors.statusUncertain
  if (status === 'normal') return colors.statusNormal
  return colors.statusUnknown
}

export function callBg(status: CallStatus): string {
  if (status === 'pathogenic') return colors.interpBg.pathogenic
  if (status === 'review') return colors.interpBg.uncertain
  if (status === 'normal') return colors.interpBg.normal
  return colors.chipBg
}

export const fmtPct = (v: number | null | undefined): string => (v == null ? '—' : `${(v * 100).toFixed(1)}%`)
export const fmtNum = (v: number | null | undefined, unit = ''): string => (v == null ? '—' : `${v.toLocaleString()}${unit}`)
