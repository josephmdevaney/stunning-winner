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
  const plotW = cW - cPad - 10
  // Place each point at its real genomic position, so a sparse depth profile
  // (e.g. depth at called sites) sits at the correct coordinates.
  const xFor = (pos: number) => cPad + (pos / genomeLength) * plotW
  const yFor = (v: number) => cH - 20 - (Math.min(v, yMax) / yMax) * (cH - 34)

  const xy = points.map((p) => ({ x: xFor(p.pos), y: yFor(p.depth) }))
  const linePts = xy.map((p) => `${p.x},${p.y}`).join(' ')
  const firstX = xy[0]?.x ?? cPad
  const lastX = xy[xy.length - 1]?.x ?? cPad
  const areaPts = `${firstX},${cH - 20} ${linePts} ${lastX},${cH - 20}`
  const sparse = points.length <= 30
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
      {/* For a sparse profile, mark each measured site so the line doesn't imply
          continuous coverage between called positions. */}
      {sparse && xy.map((p, i) => <circle key={i} cx={p.x} cy={p.y} r={2.4} fill={mito.wildType} />)}
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
