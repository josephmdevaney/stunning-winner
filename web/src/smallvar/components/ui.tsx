import type { ReactNode } from 'react'
import { colors } from '../theme'

export function SectionTitle({ title, sub }: { title: string; sub?: ReactNode }) {
  return (
    <div style={{ marginBottom: 12 }}>
      <div style={{ fontSize: 15, fontWeight: 600 }}>{title}</div>
      {sub && <div style={{ fontSize: 12, color: colors.textSecondary, marginTop: 2, lineHeight: 1.5 }}>{sub}</div>}
    </div>
  )
}

export function Pill({ children, color, bg, title }: { children: ReactNode; color: string; bg: string; title?: string }) {
  return (
    <span
      title={title}
      style={{ display: 'inline-block', padding: '2px 9px', borderRadius: 999, fontSize: 11.5, fontWeight: 600, color, background: bg, whiteSpace: 'nowrap' }}
    >
      {children}
    </span>
  )
}

export function Stat({ label, value, sub }: { label: string; value: ReactNode; sub?: ReactNode }) {
  return (
    <div style={{ minWidth: 120 }}>
      <div style={{ fontSize: 11, fontWeight: 600, color: colors.textSecondary, textTransform: 'uppercase', letterSpacing: '0.04em' }}>{label}</div>
      <div style={{ fontSize: 22, fontWeight: 700, marginTop: 2 }}>{value}</div>
      {sub && <div style={{ fontSize: 11.5, color: colors.textTertiary, marginTop: 1 }}>{sub}</div>}
    </div>
  )
}

export function LoadState({ data, error }: { data: unknown; error: string | null }) {
  if (error) return <div style={{ padding: 28, color: colors.statusPathogenic }}>Failed to load: {error}</div>
  if (!data) return <div style={{ padding: 28, color: colors.textSecondary }}>Loading…</div>
  return null
}

export function Button({
  children, onClick, primary = false, disabled = false, title,
}: { children: ReactNode; onClick: () => void; primary?: boolean; disabled?: boolean; title?: string }) {
  return (
    <button
      onClick={onClick}
      disabled={disabled}
      title={title}
      style={{
        border: `1px solid ${primary ? colors.accentButtonBg : colors.borderInput}`,
        background: primary ? colors.accentButtonBg : 'white',
        color: primary ? 'white' : colors.textPrimary,
        borderRadius: 6,
        padding: '6px 12px',
        fontSize: 12.5,
        fontWeight: 600,
        cursor: disabled ? 'default' : 'pointer',
        opacity: disabled ? 0.5 : 1,
        flexShrink: 0,
      }}
    >
      {children}
    </button>
  )
}
