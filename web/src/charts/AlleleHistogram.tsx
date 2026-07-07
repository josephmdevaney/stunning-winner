import { colors } from '../theme'

export function AlleleHistogram({
  reads,
  a1,
  a2,
}: {
  reads: { value: number; cluster: 1 | 2 }[]
  a1: number | null
  a2: number | null
}) {
  if (reads.length === 0) {
    return <div style={{ fontSize: 13, color: colors.textTertiary, padding: '20px 0' }}>No genotype call at this locus.</div>
  }

  const maxV = Math.max(...reads.map((r) => r.value), a1 ?? 0, a2 ?? 0, 1)
  const binSize = Math.max(1, Math.round(maxV / 16))
  const bins: Record<number, number> = {}
  reads.forEach((r) => {
    const b = Math.floor(r.value / binSize) * binSize
    bins[b] = (bins[b] ?? 0) + 1
  })
  const binKeys = Object.keys(bins).map(Number).sort((a, b) => a - b)
  const maxCount = Math.max(...binKeys.map((b) => bins[b]), 1)

  const hW = 360
  const hH = 160
  const hPad = 26
  const barW = (hW - hPad * 2) / (binKeys.length || 1)

  return (
    <svg width="100%" viewBox={`0 0 ${hW} ${hH}`} style={{ display: 'block' }}>
      {binKeys.map((b, i) => {
        const count = bins[b]
        const barH = (count / maxCount) * (hH - 30)
        const x = hPad + i * barW
        const y = hH - 20 - barH
        const isA1 = a1 != null && a2 != null ? Math.abs(b - a1) < Math.abs(b - a2) : true
        return <rect key={b} x={x + 1} y={y} width={Math.max(0, barW - 2)} height={barH} fill={isA1 ? colors.cluster1 : colors.cluster2} opacity={0.85} rx={1} />
      })}
      <line x1={hPad} y1={hH - 20} x2={hW - hPad} y2={hH - 20} stroke={colors.axisLine} strokeWidth={1} />
      {binKeys
        .filter((_, i) => i % Math.ceil(binKeys.length / 6 || 1) === 0)
        .map((b) => {
          const idx = binKeys.indexOf(b)
          const x = hPad + idx * barW + barW / 2
          return (
            <text key={b} x={x} y={hH - 6} fontSize={9.5} fill={colors.textSecondary} textAnchor="middle" fontFamily="IBM Plex Mono">
              {b}
            </text>
          )
        })}
    </svg>
  )
}
