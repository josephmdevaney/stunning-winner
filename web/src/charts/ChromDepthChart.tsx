import { colors } from '../theme'

export function ChromDepthChart({
  data,
}: {
  data: { chrom: string; mean_repeat_region_depth: number | null }[]
}) {
  const W = 1040
  const H = 190
  const padL = 44
  const padR = 14
  const padBottom = 28
  const values = data.map((d) => d.mean_repeat_region_depth ?? 0)
  const maxD = Math.max(...values, 10) * 1.15
  const barW = (W - padL - padR) / data.length

  const yTicks = [0, 0.25, 0.5, 0.75, 1].map((f) => Math.round(maxD * f))

  return (
    <svg width="100%" viewBox={`0 0 ${W} ${H}`} style={{ display: 'block' }}>
      {yTicks.map((v) => {
        const y = H - padBottom - (v / maxD) * (H - padBottom - 16)
        return (
          <g key={v}>
            <line x1={padL} y1={y} x2={W - padR} y2={y} stroke={colors.gridLine} />
            <text x={padL - 8} y={y + 3} fontSize={9} fill={colors.textSecondary} textAnchor="end" fontFamily="IBM Plex Mono">
              {v}x
            </text>
          </g>
        )
      })}
      {data.map((d, i) => {
        const v = d.mean_repeat_region_depth ?? 0
        const barH = (v / maxD) * (H - padBottom - 16)
        const x = padL + i * barW
        const y = H - padBottom - barH
        return (
          <rect
            key={d.chrom}
            x={x + 2}
            y={y}
            width={Math.max(0, barW - 4)}
            height={barH}
            fill={d.mean_repeat_region_depth == null ? colors.trackBg : colors.cluster1}
            rx={2}
          />
        )
      })}
      <line x1={padL} y1={H - padBottom} x2={W - padR} y2={H - padBottom} stroke={colors.axisLine} />
      {data.map((d, i) => {
        const x = padL + i * barW + barW / 2
        return (
          <text key={d.chrom} x={x} y={H - padBottom + 12} fontSize={9} fill={colors.textSecondary} textAnchor="middle" fontFamily="IBM Plex Mono">
            {d.chrom.replace('chr', '')}
          </text>
        )
      })}
    </svg>
  )
}
