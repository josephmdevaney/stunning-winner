import { mito } from '../theme'
import type { PileupRead } from '../api'

// IGV-style pileup: one row per (subsampled) HiFi read drawn as a light-gray
// span, with a colored tick at the variant column — red if the read carries the
// mutant allele, blue if wild-type. The red:blue ratio shows heteroplasmy.
export function Pileup({ reads }: { reads: PileupRead[] }) {
  const pW = 360
  const rowH = 9
  const gap = 2
  const pPad = 10
  const pH = reads.length * (rowH + gap) + pPad

  return (
    <svg width="100%" viewBox={`0 0 ${pW} ${pH}`} style={{ display: 'block' }}>
      {reads.map((r, i) => {
        const y = pPad / 2 + i * (rowH + gap)
        const x = pPad + r.x
        const midX = x + r.len / 2
        const color = r.kind === 'mut' ? mito.mutant : mito.wildType
        return (
          <g key={i}>
            <rect x={x} y={y} width={r.len} height={rowH} rx={2} fill="oklch(0.92 0.01 90)" />
            <rect x={midX - 3} y={y} width={6} height={rowH} fill={color} />
          </g>
        )
      })}
    </svg>
  )
}
