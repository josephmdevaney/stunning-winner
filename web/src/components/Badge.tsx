import { statusColor } from '../theme'

export function StatusBadge({ status }: { status: string }) {
  return (
    <span
      style={{
        padding: '4px 10px',
        borderRadius: 20,
        fontSize: 11.5,
        fontWeight: 600,
        color: 'white',
        background: statusColor(status),
      }}
    >
      {status}
    </span>
  )
}

export function FlagBadge({ flag }: { flag: string }) {
  const bg = flag === 'Pass' ? statusColor('Normal') : flag === 'No call' ? 'oklch(0.6 0.02 260)' : statusColor('Intermediate')
  return (
    <span
      style={{
        padding: '4px 10px',
        borderRadius: 20,
        fontSize: 11.5,
        fontWeight: 600,
        color: 'white',
        background: bg,
      }}
    >
      {flag}
    </span>
  )
}
