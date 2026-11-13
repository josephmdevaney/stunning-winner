import type { ReactNode } from 'react'
import { colors } from '../theme'

// A small colored label. The tint is derived from the text color so one prop
// drives the whole chip.
export function Chip({ color, children, title }: { color: string; children: ReactNode; title?: string }) {
  return (
    <span
      title={title}
      style={{
        display: 'inline-block',
        padding: '2px 8px',
        borderRadius: 20,
        fontSize: 10.5,
        fontWeight: 700,
        letterSpacing: '0.02em',
        color,
        border: `1px solid ${color}`,
        background: 'white',
        whiteSpace: 'nowrap',
      }}
    >
      {children}
    </span>
  )
}

export function Card({ children, title, sub }: { children: ReactNode; title?: string; sub?: string }) {
  return (
    <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, overflow: 'hidden' }}>
      {title && (
        <div style={{ padding: '16px 20px', borderBottom: `1px solid ${colors.borderSoft}` }}>
          <div style={{ fontSize: 15, fontWeight: 600 }}>{title}</div>
          {sub && <div style={{ fontSize: 12, color: colors.textSecondary, marginTop: 2 }}>{sub}</div>}
        </div>
      )}
      {children}
    </div>
  )
}
