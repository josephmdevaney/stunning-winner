import { colors } from '../theme'

export function Pileup({ reads }: { reads: { value: number; cluster: 1 | 2 }[] }) {
  if (reads.length === 0) {
    return <div style={{ fontSize: 13, color: colors.textTertiary, padding: '20px 0' }}>No genotype call at this locus.</div>
  }

  const sorted = [...reads].sort((a, b) => a.value - b.value)
  const maxV = Math.max(...sorted.map((r) => r.value), 1)
  const W = 860
  const rowH = 10
  const gap = 3
  const leftPad = 50
  const H = sorted.length * (rowH + gap) + 20
  const scaleX = (W - leftPad - 20) / maxV

  const ticks = [0, 0.25, 0.5, 0.75, 1].map((f) => Math.round(maxV * f))

  return (
    <svg width="100%" viewBox={`0 0 ${W} ${H + 22}`} style={{ display: 'block' }}>
      {ticks.map((t) => {
        const x = leftPad + t * scaleX
        return (
          <g key={t}>
            <line x1={x} y1={0} x2={x} y2={H} stroke={colors.gridLine} strokeWidth={1} />
            <text x={x} y={H + 14} fontSize={10} fill={colors.textSecondary} textAnchor="middle" fontFamily="IBM Plex Mono">
              {t}
            </text>
          </g>
        )
      })}
      <text
        x={4}
        y={H / 2}
        fontSize={10}
        fill={colors.textSecondary}
        transform={`rotate(-90, 4, ${H / 2})`}
        textAnchor="middle"
        fontFamily="IBM Plex Mono"
      >
        reads
      </text>
      {sorted.map((r, i) => {
        const y = 10 + i * (rowH + gap)
        const w = Math.max(2, r.value * scaleX)
        return (
          <rect
            key={i}
            x={leftPad}
            y={y}
            width={w}
            height={rowH}
            rx={2}
            fill={r.cluster === 1 ? colors.cluster1 : colors.cluster2}
            opacity={0.85}
          />
        )
      })}
    </svg>
  )
}
