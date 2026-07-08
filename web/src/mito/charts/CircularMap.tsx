import { mito, mono } from '../theme'
import type { GeneRegion, Variant } from '../api'

// Circular ideogram of the 16,569 bp mtDNA genome: a thin outer reference
// circle, colored arcs per annotated region, and a marker for the variant —
// a radial line + dot for point variants, or a thick arc for a large deletion.
export function CircularMap({
  geneMap,
  genomeLength,
  variant,
}: {
  geneMap: GeneRegion[]
  genomeLength: number
  variant: Variant
}) {
  const cx = 150
  const cy = 150
  const rOuter = 118
  const rInner = 96
  const rMid = (rOuter + rInner) / 2

  const arcPath = (startPos: number, stopPos: number, r: number) => {
    const a0 = (startPos / genomeLength) * 2 * Math.PI - Math.PI / 2
    const a1 = (stopPos / genomeLength) * 2 * Math.PI - Math.PI / 2
    const x0 = cx + r * Math.cos(a0)
    const y0 = cy + r * Math.sin(a0)
    const x1 = cx + r * Math.cos(a1)
    const y1 = cy + r * Math.sin(a1)
    const large = stopPos - startPos > genomeLength / 2 ? 1 : 0
    return `M ${x0} ${y0} A ${r} ${r} 0 ${large} 1 ${x1} ${y1}`
  }

  let marker
  if (variant.isDeletion && variant.stop != null) {
    marker = (
      <path d={arcPath(variant.pos, variant.stop, rMid)} stroke={mito.deletion} strokeWidth={rOuter - rInner + 10} fill="none" />
    )
  } else {
    const a = (variant.pos / genomeLength) * 2 * Math.PI - Math.PI / 2
    const x1 = cx + (rInner - 14) * Math.cos(a)
    const y1 = cy + (rInner - 14) * Math.sin(a)
    const x2 = cx + (rOuter + 16) * Math.cos(a)
    const y2 = cy + (rOuter + 16) * Math.sin(a)
    marker = (
      <>
        <line x1={x1} y1={y1} x2={x2} y2={y2} stroke={mito.deletion} strokeWidth={2.5} />
        <circle cx={x2} cy={y2} r={4} fill={mito.deletion} />
      </>
    )
  }

  return (
    <svg width="100%" viewBox="0 0 300 300" style={{ display: 'block' }}>
      <circle cx={cx} cy={cy} r={rOuter + 6} fill="none" stroke="oklch(0.9 0.01 90)" strokeWidth={1} />
      {geneMap.map((g, i) => (
        <path
          key={i}
          d={arcPath(g.start, g.stop, rMid)}
          stroke={mito.region[g.type] ?? 'oklch(0.8 0.01 90)'}
          strokeWidth={rOuter - rInner}
          fill="none"
          opacity={0.85}
        />
      ))}
      {marker}
      <text x={cx} y={cy - 4} fontSize={11} fill="oklch(0.5 0.02 260)" textAnchor="middle" fontFamily={mono}>
        {genomeLength.toLocaleString()} bp
      </text>
      <text x={cx} y={cy + 12} fontSize={10} fill="oklch(0.55 0.02 260)" textAnchor="middle" fontFamily={mono}>
        rCRS
      </text>
    </svg>
  )
}
