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

// Detail-screen interpretation callout, tinted by status band.
export function interpBg(status: string): string {
  if (status === 'Pathogenic') return 'oklch(0.96 0.03 25)'
  if (status === 'Low-level heteroplasmy') return 'oklch(0.96 0.03 75)'
  return 'oklch(0.96 0.03 150)'
}

export function interpColor(status: string): string {
  if (status === 'Pathogenic') return 'oklch(0.5 0.20 25)'
  if (status === 'Low-level heteroplasmy') return 'oklch(0.5 0.15 75)'
  return 'oklch(0.5 0.14 150)'
}
