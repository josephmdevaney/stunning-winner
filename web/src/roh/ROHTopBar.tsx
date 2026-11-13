import { colors } from '../theme'

export function ROHTopBar({ title, sample, subtitle }: { title: string; sample: string; subtitle: string }) {
  return (
    <div style={{ padding: '16px 28px', borderBottom: `1px solid ${colors.border}`, background: 'white' }}>
      <div style={{ fontSize: 19, fontWeight: 600 }}>{title}</div>
      <div style={{ fontSize: 12.5, color: colors.textSecondary, marginTop: 2 }}>
        Sample {sample || '—'} · {subtitle}
      </div>
    </div>
  )
}
