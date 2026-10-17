import type { CSSProperties } from 'react'
import { colors } from '../theme'

export const cardStyle: CSSProperties = {
  background: 'white',
  border: `1px solid ${colors.border}`,
  borderRadius: 12,
  padding: '20px 24px',
}

export const pageStyle: CSSProperties = { padding: 28, display: 'flex', flexDirection: 'column', gap: 20 }
