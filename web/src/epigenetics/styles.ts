import { colors } from '../theme'

// Table cell styles shared by the Epigenetics screens.
export const th = {
  padding: '10px 14px',
  fontSize: 10.5,
  fontWeight: 600,
  color: colors.textSecondary,
  letterSpacing: '0.04em',
  textTransform: 'uppercase' as const,
  textAlign: 'left' as const,
  borderBottom: `1px solid ${colors.borderSofter}`,
  whiteSpace: 'nowrap' as const,
}

export const td = {
  padding: '10px 14px',
  fontSize: 13,
  borderBottom: `1px solid ${colors.borderSoftest}`,
  verticalAlign: 'middle' as const,
}
