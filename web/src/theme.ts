// Color tokens lifted directly from the design bundle's .dc.html prototypes
// (project/Repeat Expansion Analyzer.dc.html and its imported screens).
export const colors = {
  pageBg: 'oklch(0.97 0.005 90)',
  textPrimary: 'oklch(0.22 0.02 260)',

  sidebarBg: 'oklch(0.22 0.04 260)',
  sidebarHover: 'oklch(0.28 0.05 260)',
  sidebarPanelBg: 'oklch(0.27 0.05 260)',
  sidebarTextMuted: 'oklch(0.6 0.02 260)',
  sidebarTextItem: 'oklch(0.85 0.02 260)',
  sidebarTextSub: 'oklch(0.68 0.03 220)',
  sidebarDotInactive: 'oklch(0.4 0.02 260)',

  accent: 'oklch(0.78 0.13 220)',
  accentButtonBg: 'oklch(0.5 0.16 220)',

  border: 'oklch(0.9 0.01 90)',
  borderSoft: 'oklch(0.92 0.01 90)',
  borderSofter: 'oklch(0.94 0.01 90)',
  borderSoftest: 'oklch(0.96 0.005 90)',
  borderInput: 'oklch(0.85 0.01 90)',
  borderInputAlt: 'oklch(0.88 0.01 90)',
  gridLine: 'oklch(0.93 0.005 90)',
  axisLine: 'oklch(0.85 0.01 90)',

  textSecondary: 'oklch(0.5 0.02 260)',
  textTertiary: 'oklch(0.55 0.02 260)',
  textQuaternary: 'oklch(0.45 0.02 260)',
  textLabel: 'oklch(0.4 0.02 260)',
  textBody: 'oklch(0.35 0.02 260)',
  chipBg: 'oklch(0.95 0.01 90)',

  statusPathogenic: 'oklch(0.55 0.20 25)',
  statusUncertain: 'oklch(0.5 0.15 300)',
  statusIntermediate: 'oklch(0.62 0.15 75)',
  statusNormal: 'oklch(0.55 0.14 150)',
  statusUnknown: 'oklch(0.6 0.02 260)',

  cluster1: 'oklch(0.6 0.14 220)',
  cluster2: 'oklch(0.6 0.18 25)',

  methylated: 'oklch(0.55 0.18 300)',
  unmethylated: 'oklch(0.8 0.01 90)',
  trackBg: 'oklch(0.93 0.01 90)',

  interpBg: {
    pathogenic: 'oklch(0.96 0.03 25)',
    uncertain: 'oklch(0.96 0.02 300)',
    intermediate: 'oklch(0.96 0.03 75)',
    normal: 'oklch(0.96 0.03 150)',
  },
  interpColor: {
    pathogenic: 'oklch(0.5 0.20 25)',
    uncertain: 'oklch(0.5 0.15 300)',
    intermediate: 'oklch(0.5 0.15 75)',
    normal: 'oklch(0.5 0.14 150)',
  },
} as const

export const mono = "'IBM Plex Mono', monospace"

export function statusColor(status: string): string {
  if (status === 'Full expansion') return colors.statusPathogenic
  if (status === 'Uncertain significance') return colors.statusUncertain
  if (status === 'Intermediate') return colors.statusIntermediate
  if (status === 'Normal') return colors.statusNormal
  return colors.statusUnknown
}

export function interpBg(status: string): string {
  if (status === 'Full expansion') return colors.interpBg.pathogenic
  if (status === 'Uncertain significance') return colors.interpBg.uncertain
  if (status === 'Intermediate') return colors.interpBg.intermediate
  return colors.interpBg.normal
}

export function interpColor(status: string): string {
  if (status === 'Full expansion') return colors.interpColor.pathogenic
  if (status === 'Uncertain significance') return colors.interpColor.uncertain
  if (status === 'Intermediate') return colors.interpColor.intermediate
  return colors.interpColor.normal
}
