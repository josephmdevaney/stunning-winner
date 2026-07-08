import { mito, mono } from '../theme'

// Bar histogram of HiFi read lengths (kb), bell-shaped around the full-length
// ~16.5 kb mtDNA amplicon.
export function ReadLengthHistogram({ bins }: { bins: { kb: number; count: number }[] }) {
  const hW = 380
  const hH = 170
  const hPad = 28
  const maxB = Math.max(...bins.map((b) => b.count), 1)
  const barW = (hW - hPad * 2) / bins.length
  const axisTicks = [0, 8, 17]

  return (
    <svg width="100%" viewBox={`0 0 ${hW} ${hH}`} style={{ display: 'block' }}>
      {bins.map((b, i) => {
        const barH = (b.count / maxB) * (hH - 30)
        const x = hPad + i * barW
        const y = hH - 20 - barH
        return <rect key={i} x={x + 1} y={y} width={barW - 2} height={barH} fill={mito.wildType} opacity={0.85} rx={1} />
      })}
      <line x1={hPad} y1={hH - 20} x2={hW - hPad} y2={hH - 20} stroke="oklch(0.85 0.01 90)" />
      {axisTicks.map((i) => {
        const x = hPad + i * barW
        const kb = bins[i]?.kb ?? 2 + i * 1.1
        return (
          <text key={i} x={x} y={hH - 6} fontSize={9.5} fill="oklch(0.5 0.02 260)" textAnchor="middle" fontFamily={mono}>
            {kb.toFixed(0)}kb
          </text>
        )
      })}
    </svg>
  )
}
