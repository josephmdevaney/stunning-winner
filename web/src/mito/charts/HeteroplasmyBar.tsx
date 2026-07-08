import { mito, mono } from '../theme'

// Single horizontal stacked bar: blue wild-type fraction, red mutant fraction
// (width ∝ heteroplasmy), with a dashed 5% reporting-threshold marker.
export function HeteroplasmyBar({ het }: { het: number }) {
  const bw = 700
  const bh = 46
  const mutW = (het / 100) * bw
  const threshX = (5 / 100) * bw

  return (
    <svg width="100%" viewBox={`0 0 ${bw} ${bh + 22}`} style={{ display: 'block' }}>
      <rect x={0} y={0} width={bw} height={bh} rx={8} fill={mito.wildType} />
      <rect x={0} y={0} width={mutW} height={bh} rx={8} fill={mito.mutant} />
      <line x1={threshX} y1={-4} x2={threshX} y2={bh + 4} stroke={mito.threshold} strokeWidth={1.5} strokeDasharray="4,3" />
      <text x={10} y={bh / 2 + 5} fontSize={13} fontWeight={700} fill="white" fontFamily={mono}>
        {het}%
      </text>
      <text x={bw - 10} y={bh / 2 + 5} fontSize={13} fontWeight={700} fill="white" textAnchor="end" fontFamily={mono}>
        {(100 - het).toFixed(1)}%
      </text>
      <text x={threshX} y={bh + 18} fontSize={10} fill="oklch(0.5 0.02 260)" textAnchor="middle" fontFamily={mono}>
        5% threshold
      </text>
    </svg>
  )
}
