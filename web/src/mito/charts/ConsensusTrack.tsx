import { mito, mono } from '../theme'

// Two-row sequence track: top = reference (rCRS) bases, bottom = sample
// consensus. Every column matches except the variant column, where the sample
// tile is highlighted and shows the alternate base.
const REF_BASES = 'ATCGGATCCGTAAGCTTGACCTGA'.split('')

export function ConsensusTrack({ altBase }: { altBase: string }) {
  const variantIdx = Math.floor(REF_BASES.length / 2)
  const cW = 720
  const cH = 90
  const cellW = cW / REF_BASES.length

  return (
    <svg width="100%" viewBox={`0 0 ${cW} ${cH}`} style={{ display: 'block' }}>
      {REF_BASES.map((b, i) => {
        const x = i * cellW
        const isVar = i === variantIdx
        return (
          <g key={`ref-${i}`}>
            <rect
              x={x + 1}
              y={6}
              width={cellW - 2}
              height={28}
              rx={3}
              fill={isVar ? 'oklch(0.93 0.05 220)' : 'oklch(0.97 0.005 90)'}
            />
            <text x={x + cellW / 2} y={25} fontSize={13} fontWeight={600} textAnchor="middle" fill="oklch(0.3 0.02 260)" fontFamily={mono}>
              {b}
            </text>
          </g>
        )
      })}
      {REF_BASES.map((b, i) => {
        const x = i * cellW
        const isVar = i === variantIdx
        const shown = isVar ? altBase : b
        return (
          <g key={`sample-${i}`}>
            <rect x={x + 1} y={50} width={cellW - 2} height={28} rx={3} fill={isVar ? mito.mutant : 'oklch(0.96 0.01 90)'} />
            <text
              x={x + cellW / 2}
              y={69}
              fontSize={13}
              fontWeight={700}
              textAnchor="middle"
              fill={isVar ? 'white' : 'oklch(0.3 0.02 260)'}
              fontFamily={mono}
            >
              {shown}
            </text>
          </g>
        )
      })}
    </svg>
  )
}
