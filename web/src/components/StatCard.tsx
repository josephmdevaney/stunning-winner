import type { CSSProperties } from 'react'
import { colors } from '../theme'

export function StatGrid({ children }: { children: React.ReactNode }) {
  return (
    <div style={{ display: 'grid', gridTemplateColumns: 'repeat(4,1fr)', gap: 16 }}>{children}</div>
  )
}

export function StatCard({
  label,
  value,
  sub,
  valueStyle,
  dotColor,
}: {
  label: string
  value: string
  sub?: string
  valueStyle?: CSSProperties
  dotColor?: string
}) {
  return (
    <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, padding: '18px 20px' }}>
      <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
        <div style={{ fontSize: 12, fontWeight: 600, color: colors.textSecondary }}>{label}</div>
        {dotColor && <div style={{ width: 8, height: 8, borderRadius: '50%', background: dotColor }} />}
      </div>
      <div style={{ fontSize: 26, fontWeight: 700, marginTop: 8, ...valueStyle }}>{value}</div>
      {sub && <div style={{ fontSize: 12, color: colors.textTertiary, marginTop: 4 }}>{sub}</div>}
    </div>
  )
}
