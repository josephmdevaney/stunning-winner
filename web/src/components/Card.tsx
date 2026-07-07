import type { CSSProperties, ReactNode } from 'react'
import { colors } from '../theme'

export function Card({
  children,
  style,
  padded = true,
}: {
  children: ReactNode
  style?: CSSProperties
  padded?: boolean
}) {
  return (
    <div
      style={{
        background: 'white',
        border: `1px solid ${colors.border}`,
        borderRadius: 12,
        ...(padded ? { padding: '20px 24px' } : {}),
        ...style,
      }}
    >
      {children}
    </div>
  )
}

export function CardHeading({ title, sub }: { title: string; sub?: string }) {
  return (
    <>
      <div style={{ fontSize: 14.5, fontWeight: 600, marginBottom: 2 }}>{title}</div>
      {sub && <div style={{ fontSize: 12, color: colors.textTertiary, marginBottom: 14 }}>{sub}</div>}
    </>
  )
}
