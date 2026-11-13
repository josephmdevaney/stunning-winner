import { colors } from '../theme'

export function EpigeneticsTopBar({ title, subtitle }: { title: string; subtitle: string }) {
  return (
    <div
      style={{
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'space-between',
        padding: '16px 28px',
        borderBottom: `1px solid ${colors.border}`,
        background: 'white',
      }}
    >
      <div>
        <div style={{ fontSize: 19, fontWeight: 600 }}>{title}</div>
        <div style={{ fontSize: 12.5, color: colors.textSecondary, marginTop: 2 }}>{subtitle}</div>
      </div>
      <div style={{ fontSize: 11.5, color: colors.textTertiary, fontWeight: 600, letterSpacing: '0.03em' }}>
        RESEARCH USE ONLY
      </div>
    </div>
  )
}
