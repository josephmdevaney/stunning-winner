import type { ReactNode } from 'react'
import { colors } from '../theme'
import type { Category, RegionResult } from './api'
import { CATEGORY_HELP, CATEGORY_LABEL, STATUS_COLOR, STATUS_LABEL, TRACK_COLOR } from './labels'

export function Chip({ color, children, title, solid }: { color: string; children: ReactNode; title?: string; solid?: boolean }) {
  return (
    <span
      title={title}
      style={{
        display: 'inline-block',
        padding: '2px 9px',
        borderRadius: 20,
        fontSize: 11,
        fontWeight: 700,
        letterSpacing: '0.02em',
        color: solid ? 'white' : color,
        border: `1px solid ${color}`,
        background: solid ? color : 'white',
        whiteSpace: 'nowrap',
      }}
    >
      {children}
    </span>
  )
}

export function CategoryChip({ category }: { category: Category }) {
  return (
    <Chip color={colors.statusUnknown} title={CATEGORY_HELP[category]}>
      {CATEGORY_LABEL[category]}
    </Chip>
  )
}

// A polymorphic (X) region's methylation isn't interpretable for imprinting disorders, so its observed pattern is
// shown as data but never as a finding.
export function RegionStatusChip({ region }: { region: RegionResult }) {
  if (!region.interpretable) {
    return (
      <Chip color={colors.statusUnknown} title={`${CATEGORY_HELP.X}. Observed: ${STATUS_LABEL[region.status]}.`}>
        Not interpreted
      </Chip>
    )
  }
  // Only a clinically associated DMR gets the filled chip; a research DMR's finding is outlined so it isn't read as diagnostic.
  const filled = region.category === 'CA' && region.status !== 'no_data' && region.status !== 'low_coverage'
  const title = [region.pattern, region.category === 'NC' ? CATEGORY_HELP.NC : null].filter(Boolean).join(' · ')
  return (
    <Chip color={STATUS_COLOR[region.status]} solid={filled} title={title || undefined}>
      {STATUS_LABEL[region.status]}
    </Chip>
  )
}

// Methylation 0-100% as a bar with its value beside it.
export function MethBar({ value, color = TRACK_COLOR.combined, width = 88 }: { value: number | null | undefined; color?: string; width?: number }) {
  if (value == null) return <span style={{ color: colors.textTertiary }}>—</span>
  return (
    <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
      <div style={{ width, height: 8, borderRadius: 4, background: colors.trackBg, overflow: 'hidden', flexShrink: 0 }}>
        <div style={{ width: `${Math.max(0, Math.min(100, value))}%`, height: '100%', background: color }} />
      </div>
      <span style={{ fontSize: 12.5, fontVariantNumeric: 'tabular-nums', minWidth: 40 }}>{value.toFixed(1)}%</span>
    </div>
  )
}

export function Message({ children, tone }: { children: ReactNode; tone?: 'error' }) {
  return (
    <div style={{ padding: 28 }}>
      <div
        style={{
          background: 'white',
          border: `1px solid ${colors.border}`,
          borderRadius: 12,
          padding: '60px 24px',
          textAlign: 'center',
          color: tone === 'error' ? colors.statusPathogenic : colors.textSecondary,
          fontSize: 14,
        }}
      >
        {children}
      </div>
    </div>
  )
}
