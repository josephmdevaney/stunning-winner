import { colors } from '../theme'
import type { Category, DisorderStatus, EpiCall, Status } from './api'

export const STATUS_LABEL: Record<Status, string> = {
  imprinted: 'Imprint intact',
  loss: 'Loss of methylation',
  gain: 'Gain of methylation',
  borderline: 'Borderline',
  low_coverage: 'Low coverage',
  no_data: 'No data',
}

// Green is only for an intact imprint; red for a disturbance in either direction; amber for borderline;
// gray when there isn't enough data to say.
export const STATUS_COLOR: Record<Status, string> = {
  imprinted: colors.statusNormal,
  loss: colors.statusPathogenic,
  gain: colors.statusPathogenic,
  borderline: colors.statusIntermediate,
  low_coverage: colors.statusUnknown,
  no_data: colors.statusUnknown,
}

export const DISORDER_STATUS_LABEL: Record<DisorderStatus, string> = {
  pattern_present: 'Pattern present',
  supportive_only: 'Supportive only',
  not_observed: 'No pattern',
  insufficient_data: 'Insufficient data',
}

export const DISORDER_STATUS_COLOR: Record<DisorderStatus, string> = {
  pattern_present: colors.statusPathogenic,
  supportive_only: colors.statusIntermediate,
  not_observed: colors.statusNormal,
  insufficient_data: colors.statusUnknown,
}

export const CATEGORY_LABEL: Record<Category, string> = {
  CA: 'Clinical',
  NC: 'Research',
  X: 'Polymorphic',
}

export const CATEGORY_HELP: Record<Category, string> = {
  CA: 'Clinically associated: used in imprinting-disorder diagnosis',
  NC: 'Non-clinical: imprinted, but not used in diagnosis (research only)',
  X: 'Polymorphic or tissue-variable: methylation here is not interpretable for imprinting disorders',
}

// Hap tracks are named by the phasing, not by parent of origin; two fixed, distinguishable colors.
export const TRACK_COLOR = {
  combined: 'oklch(0.3 0.02 260)',
  hap1: 'oklch(0.55 0.15 250)',
  hap2: 'oklch(0.65 0.17 55)',
} as const

export const fmtPct = (v: number | null | undefined): string => (v == null ? '—' : `${v.toFixed(1)}%`)
export const fmtPos = (n: number): string => n.toLocaleString('en-US')
export const mean = (v: number[]): number | null => (v.length ? v.reduce((a, b) => a + b, 0) / v.length : null)

// ---- Episignatures ----

export const EPI_CALL_LABEL: Record<EpiCall, string> = {
  positive: 'Matches signature',
  indeterminate: 'Indeterminate',
  negative: 'No signature',
  insufficient_data: 'Not enough data',
}

// A match is red and a clear negative green, as for the imprinting calls; the grey is for "no call".
export const EPI_CALL_COLOR: Record<EpiCall, string> = {
  positive: colors.statusPathogenic,
  indeterminate: colors.statusIntermediate,
  negative: colors.statusNormal,
  insufficient_data: colors.statusUnknown,
}

// Training controls and cases, and the genome being scored.
export const EPI_COLOR = {
  control: 'oklch(0.58 0.12 230)',
  case: 'oklch(0.58 0.19 25)',
  sample: 'oklch(0.2 0.02 260)',
} as const
