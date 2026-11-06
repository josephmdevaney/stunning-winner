import { colors } from '../theme'

// Column-header row shared by the Paraphase tables.
export const headerRowStyle = {
  padding: '10px 20px',
  fontSize: 10.5,
  fontWeight: 600,
  color: colors.textSecondary,
  letterSpacing: '0.04em',
  textTransform: 'uppercase' as const,
  borderBottom: `1px solid ${colors.borderSofter}`,
}
