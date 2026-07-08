import { mito, mono } from '../theme'

// Filled line/area chart of read depth vs. genomic position (1–16,569). The
// Y-axis max is configurable per report; gridlines/labels are drawn at 0/50/100%
// of it and the plotted depth line is clipped to it.
export function CoverageTrack({
  points,
  genomeLength,
  yMax,
}: {
  points: { pos: number; depth: number }[]
  genomeLength: number
  yMax: number
}) {
  const cW = 380
  const cH = 170
  const cPad = 34
  const nPts = points.length
  const stepX = (cW - cPad - 10) / Math.max(1, nPts - 1)
  const yFor = (v: number) => cH - 20 - (Math.min(v, yMax) / yMax) * (cH - 34)

  const linePts = points.map((p, i) => `${cPad + i * stepX},${yFor(p.depth)}`).join(' ')
  const areaPts = `${cPad},${cH - 20} ${linePts} ${cPad + (nPts - 1) * stepX},${cH - 20}`
  const xLabels = [0, Math.round(genomeLength / 2), genomeLength]

  return (
    <svg width="100%" viewBox={`0 0 ${cW} ${cH}`} style={{ display: 'block' }}>
      {[0, 0.5, 1].map((f, i) => {
        const val = Math.round(yMax * f)
        const y = yFor(val)
        return (
          <g key={i}>
            <line x1={cPad} y1={y} x2={cW - 10} y2={y} stroke="oklch(0.93 0.005 90)" strokeWidth={1} />
            <text x={cPad - 4} y={y + 3} fontSize={8.5} fill="oklch(0.5 0.02 260)" textAnchor="end" fontFamily={mono}>
              {val}
            </text>
          </g>
        )
      })}
      <polygon points={areaPts} fill={mito.wildType} opacity={0.25} />
      <polyline points={linePts} fill="none" stroke={mito.wildType} strokeWidth={1.8} />
      <line x1={cPad} y1={cH - 20} x2={cW - 10} y2={cH - 20} stroke="oklch(0.85 0.01 90)" />
      {xLabels.map((p, i) => {
        const x = cPad + (i * (cW - cPad - 10)) / 2
        return (
          <text key={i} x={x} y={cH - 6} fontSize={9} fill="oklch(0.5 0.02 260)" textAnchor="middle" fontFamily={mono}>
            {p}
          </text>
        )
      })}
    </svg>
  )
}
