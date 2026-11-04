import { useEffect, useState } from 'react'
import { structvarApi, type RegionGene, type RegionResponse, type Sv } from '../api'
import { colors, fmtSize, mono, TYPE_COLOR } from '../theme'

const W = 940
const LEFT = 8
const ROW = 30
const MAX_REGION = 20_000_000

// The SV drawn over the genes around it: each gene as its MANE (or fallback) transcript, with exons as
// boxes (taller where coding), and the proband's other calls as ticks.
export function RegionView({ sample, sv }: { sample: string; sv: Sv }) {
  const span = Math.max(sv.end - sv.pos, 0)
  const pad = Math.max(span * 0.3, 20_000)
  const start = Math.max(1, Math.round(sv.pos - pad))
  const end = Math.round(sv.end + pad)
  const tooBig = end - start > MAX_REGION
  const [data, setData] = useState<RegionResponse | null>(null)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    if (tooBig) return
    let current = true
    setData(null)
    setError(null)
    structvarApi.region(sample, sv.chrom, start, end).then((r) => current && setData(r)).catch((e) => current && setError(String(e)))
    return () => {
      current = false
    }
  }, [sample, sv.chrom, start, end, tooBig])

  if (tooBig) return <div style={{ fontSize: 12.5, color: colors.textSecondary }}>The event is too large to draw ({fmtSize(span)}).</div>
  if (error) return <div style={{ fontSize: 12.5, color: colors.statusPathogenic }}>Could not load the region: {error}</div>
  if (!data) return <div style={{ fontSize: 12.5, color: colors.textSecondary }}>Loading region…</div>

  const x = (p: number) => LEFT + ((p - start) / (end - start)) * (W - 2 * LEFT)
  // pack genes into rows so overlapping ones do not draw on top of each other
  const rows: number[] = []
  const placed = data.genes.map((g) => {
    let r = rows.findIndex((last) => last < g.start - (end - start) * 0.02)
    if (r < 0) {
      r = rows.length
      rows.push(0)
    }
    rows[r] = g.end
    return { g, r }
  })
  const nRows = Math.min(Math.max(rows.length, 1), 8)
  const top = 34
  const H = top + nRows * ROW + 40
  const color = TYPE_COLOR[sv.svtype] ?? colors.textSecondary
  const ticks = [0, 0.25, 0.5, 0.75, 1].map((f) => Math.round(start + f * (end - start)))
  const point = sv.svtype === 'INS' || sv.svtype === 'BND'

  return (
    <div style={{ overflowX: 'auto' }}>
      <svg width={W} height={H} style={{ display: 'block', fontFamily: mono }}>
        {/* the SV */}
        {point ? (
          <g>
            <line x1={x(sv.pos)} x2={x(sv.pos)} y1={6} y2={H - 28} stroke={color} strokeWidth={1.5} strokeDasharray="3 3" />
            <polygon points={`${x(sv.pos) - 6},6 ${x(sv.pos) + 6},6 ${x(sv.pos)},18`} fill={color} />
          </g>
        ) : (
          <g>
            <rect x={x(sv.pos)} y={8} width={Math.max(2, x(sv.end) - x(sv.pos))} height={12} rx={2} fill={color} opacity={sv.imprecise ? 0.55 : 0.9} />
            <line x1={x(sv.pos)} x2={x(sv.pos)} y1={20} y2={H - 28} stroke={color} strokeWidth={1} strokeDasharray="3 3" opacity={0.7} />
            <line x1={x(sv.end)} x2={x(sv.end)} y1={20} y2={H - 28} stroke={color} strokeWidth={1} strokeDasharray="3 3" opacity={0.7} />
          </g>
        )}
        <text x={Math.min(Math.max(x(sv.pos), 70), W - 70)} y={30} fontSize={10.5} fill={colors.textSecondary} textAnchor="middle">
          {sv.svtype} {fmtSize(sv.size)}
        </text>

        {placed.slice(0, 40).map(({ g, r }) => r < nRows && <GeneGlyph key={g.gene_id} g={g} y={top + r * ROW} x={x} />)}

        {/* the proband's other calls */}
        {data.svs.filter((s) => s.id !== sv.id).map((s) => (
          <rect key={s.id} x={x(s.pos)} y={H - 24} width={Math.max(1, x(s.end) - x(s.pos))} height={5} fill={TYPE_COLOR[s.svtype] ?? colors.textSecondary} opacity={0.5} />
        ))}
        <rect x={x(sv.pos)} y={H - 24} width={Math.max(2, x(sv.end) - x(sv.pos))} height={5} fill={color} />

        {/* axis */}
        <line x1={LEFT} x2={W - LEFT} y1={H - 12} y2={H - 12} stroke={colors.axisLine} />
        {ticks.map((t) => (
          <g key={t}>
            <line x1={x(t)} x2={x(t)} y1={H - 12} y2={H - 8} stroke={colors.axisLine} />
            <text x={x(t)} y={H - 0} fontSize={9.5} fill={colors.textSecondary} textAnchor={t === ticks[0] ? 'start' : t === ticks[4] ? 'end' : 'middle'}>
              {sv.chrom}:{t.toLocaleString()}
            </text>
          </g>
        ))}
      </svg>
      <div style={{ fontSize: 11.5, color: colors.textTertiary, marginTop: 4 }}>
        {data.n_genes} gene{data.n_genes === 1 ? '' : 's'} in view{data.n_genes > 40 ? ' (first 40 drawn)' : ''}. Exons are those of the MANE Select transcript where the gene has one (taller where coding).
        Faint ticks along the bottom are the proband&rsquo;s other PASS calls.
      </div>
    </div>
  )
}

function GeneGlyph({ g, y, x }: { g: RegionGene; y: number; x: (p: number) => number }) {
  const x0 = x(g.start)
  const x1 = x(g.end)
  const mid = y + 10
  const arrow = g.strand === '-' ? '‹' : '›'
  const chevrons: number[] = []
  for (let px = x0 + 14; px < x1 - 6; px += 40) chevrons.push(px)
  return (
    <g>
      <line x1={x0} x2={x1} y1={mid} y2={mid} stroke={colors.textSecondary} strokeWidth={1} />
      {chevrons.map((px) => (
        <text key={px} x={px} y={mid + 3.5} fontSize={11} fill={colors.textSecondary} textAnchor="middle">
          {arrow}
        </text>
      ))}
      {g.exons.map(([s, e]) => {
        const coding = g.cds != null && e >= g.cds[0] && s <= g.cds[1]
        const h = coding ? 12 : 6
        return <rect key={s} x={x(s)} y={mid - h / 2} width={Math.max(1.5, x(e) - x(s))} height={h} fill={coding ? colors.textBody : colors.textTertiary} opacity={0.85} />
      })}
      <text x={Math.max(x0, 4)} y={y + 26} fontSize={10.5} fontWeight={600} fill={colors.textPrimary}>
        {g.symbol}
      </text>
    </g>
  )
}
