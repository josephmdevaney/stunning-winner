// MitoLR reuses RepeatScope's base tokens (colors, mono) and adds the
// mtDNA-specific status semantics and chart accents from the design bundle.
import { colors } from '../theme'

export { colors, mono } from '../theme'

export const mito = {
  // Chart accents (heteroplasmy bar, pileup, coverage, histogram).
  wildType: 'oklch(0.6 0.14 220)',
  mutant: 'oklch(0.6 0.18 25)',
  deletion: 'oklch(0.55 0.20 25)',
  threshold: 'oklch(0.4 0.02 260)',

  // Circular-map region colors, keyed by gene-map `type`.
  region: {
    control: 'oklch(0.75 0.02 90)',
    rRNA: 'oklch(0.7 0.1 300)',
    protein: 'oklch(0.7 0.1 220)',
  } as Record<string, string>,

  // Haplogroup marker chip.
  markerBg: 'oklch(0.96 0.01 220)',
  markerText: 'oklch(0.4 0.1 220)',

  // Overview origin filter chip (inactive).
  chipInactiveBg: colors.chipBg,
} as const

// Filter (VCF FILTER column) text color: PASS green, anything else amber.
export function filterColor(filter: string): string {
  return filter === 'PASS' ? 'oklch(0.5 0.14 150)' : 'oklch(0.6 0.15 75)'
}

// Detail-screen interpretation callout, tinted by status. Green is reserved for
// results that are reassuring (Common variant, Not detected); a VUS and a Rare
// variant are neither, so each gets its own tint (matching its badge) instead of
// falling through to green.
const INTERP_TINT: Record<string, { bg: string; fg: string }> = {
  Pathogenic: { bg: 'oklch(0.96 0.03 25)', fg: 'oklch(0.5 0.20 25)' },
  'Low-level heteroplasmy': { bg: 'oklch(0.96 0.03 75)', fg: 'oklch(0.5 0.15 75)' },
  VUS: { bg: colors.interpBg.uncertain, fg: colors.interpColor.uncertain },
  'Rare variant': { bg: 'oklch(0.96 0.02 230)', fg: 'oklch(0.5 0.12 230)' },
}
const REASSURING_TINT = { bg: 'oklch(0.96 0.03 150)', fg: 'oklch(0.5 0.14 150)' }

export function interpBg(status: string): string {
  return (INTERP_TINT[status] ?? REASSURING_TINT).bg
}

export function interpColor(status: string): string {
  return (INTERP_TINT[status] ?? REASSURING_TINT).fg
}
