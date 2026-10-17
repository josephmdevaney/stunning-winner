import type { ReactNode } from 'react'
import { callBg, callColor, colors } from '../theme'
import type { CallStatus } from '../api'
import { cardStyle } from './styles'

export function StatusPill({ status, label }: { status: CallStatus; label: string }) {
  return (
    <span
      style={{
        display: 'inline-block',
        padding: '3px 10px',
        borderRadius: 999,
        fontSize: 11.5,
        fontWeight: 600,
        color: callColor(status),
        background: callBg(status),
        whiteSpace: 'nowrap',
      }}
    >
      {label}
    </span>
  )
}

export function SectionTitle({ title, sub }: { title: string; sub?: ReactNode }) {
  return (
    <div style={{ marginBottom: 12 }}>
      <div style={{ fontSize: 15, fontWeight: 600 }}>{title}</div>
      {sub && <div style={{ fontSize: 12, color: colors.textSecondary, marginTop: 2 }}>{sub}</div>}
    </div>
  )
}

// The loading/error placeholder every Kivvi screen shows until its data arrives.
export function LoadState({ data, error }: { data: unknown; error: string | null }) {
  if (error) return <div style={{ padding: 28, color: colors.statusPathogenic }}>Failed to load: {error}</div>
  if (!data) return <div style={{ padding: 28, color: colors.textSecondary }}>Loading…</div>
  return null
}

export function Notes({ notes }: { notes: string[] }) {
  if (!notes.length) return null
  return (
    <div style={{ ...cardStyle, padding: '12px 18px', background: colors.interpBg.intermediate, fontSize: 12.5, color: colors.interpColor.intermediate }}>
      {notes.map((n) => (
        <div key={n}>{n}</div>
      ))}
    </div>
  )
}
