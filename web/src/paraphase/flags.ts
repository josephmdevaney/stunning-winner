import { colors } from '../theme'
import type { FlagKind, ToolStatus, VersionGap } from './api'

// Flags only restate something Paraphase's own JSON says, so none of them is
// "pathogenic" red: purple for a notable structural call, amber for something
// unresolved or failed, gray for no data.
export const FLAG_COLOR: Record<FlagKind, string> = {
  fusion: colors.statusUncertain,
  sv: colors.statusUncertain,
  phasing_failed: colors.statusIntermediate,
  cn_unresolved: colors.statusIntermediate,
  no_data: colors.textTertiary,
}

export const FLAG_LABEL: Record<FlagKind, string> = {
  fusion: 'Fusion',
  sv: 'Structural variant',
  phasing_failed: 'Phasing failed',
  cn_unresolved: 'Copy number not reported',
  no_data: 'No data',
}

export const GAP_COLOR: Record<string, string> = {
  none: colors.statusNormal,
  ahead: colors.statusNormal,
  patch: colors.statusIntermediate,
  minor: colors.statusIntermediate,
  major: colors.statusPathogenic,
}

export const GAP_LABEL: Record<string, string> = {
  none: 'Current',
  ahead: 'Newer than latest release',
  patch: 'Patch behind',
  minor: 'Minor version behind',
  major: 'Major version behind',
}

export function gapColor(gap: VersionGap): string {
  return gap ? GAP_COLOR[gap] : colors.textTertiary
}

export const STATUS_LABEL: Record<ToolStatus, string> = {
  current: 'Current',
  behind: 'Behind',
  unknown: 'Version not in files',
  not_recorded: 'Not recorded',
  not_used: 'Not used here',
}

export function statusColor(status: ToolStatus, gap: VersionGap): string {
  if (status === 'current') return colors.statusNormal
  if (status === 'behind') return gapColor(gap)
  return colors.textTertiary
}
