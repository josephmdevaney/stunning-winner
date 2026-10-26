import { useState } from 'react'
import { Card, CardHeading } from '../../components/Card'
import { colors, mono } from '../../theme'
import type { ArrayCondition, ArrayOverview, ArrayResponse } from '../api'
import { EPI_CALL_COLOR, EPI_COLOR, fmtPos } from '../labels'
import { Callout, Chip } from '../ui'

export function ArraySignatures({ data, error }: { data: ArrayResponse | null; error: string | null }) {
  if (error && !data) return <Callout tone="caution" title="Array-derived classifiers">{error}</Callout>
  if (!data) {
    return (
      <Card>
        <CardHeading title="Array-derived classifiers (34 conditions)" sub="Scoring: each genome retrains 34 SVMs, which takes a few seconds." />
      </Card>
    )
  }
  if (!data.installed) {
    return (
      <Card>
        <CardHeading title="Array-derived classifiers (34 conditions)" sub="Not set up. The setup command also installs scikit-learn into a small Python of its own." />
        <code style={{ display: 'block', padding: '10px 14px', background: colors.pageBg, borderRadius: 8, fontFamily: mono, fontSize: 12.5 }}>{data.setup_command}</code>
      </Card>
    )
  }
  return <ArrayResults data={data} />
}

function ArrayResults({ data }: { data: ArrayOverview }) {
  const [picked, setPicked] = useState<string | null>(null)
  const selected = data.conditions.find((c) => c.id === (picked ?? data.matches[0] ?? data.conditions[0].id)) ?? data.conditions[0]
  const lo = Math.min(-0.5, ...data.conditions.map((c) => c.score))
  const hi = Math.max(1.2, ...data.conditions.map((c) => c.score))
  const x = (v: number) => ((v - lo) / (hi - lo)) * 100

  const headline =
    data.call === 'control'
      ? `No condition reaches ${data.threshold}`
      : data.call === 'one'
        ? `${data.matches[0]} reaches ${data.threshold}`
        : `${data.matches.length} conditions reach ${data.threshold}: ${data.matches.join(', ')}`

  return (
    <>
      <Card>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', gap: 16 }}>
          <CardHeading title="Array-derived classifiers (34 conditions)" sub="Support-vector machines trained on published median array methylation, one profile per condition, over each condition's signature CpGs." />
          <Chip color={data.call === 'control' ? EPI_CALL_COLOR.negative : EPI_CALL_COLOR.positive} solid>
            {headline}
          </Chip>
        </div>
        <div style={{ fontSize: 12, color: colors.textSecondary, marginBottom: 14 }}>
          {data.cpgs.found.toLocaleString()} of {data.cpgs.total.toLocaleString()} signature CpGs had a row in this pileup; a CpG without one takes the reference mean.
        </div>

        {data.warnings.map((w) => (
          <div key={w.id} style={{ fontSize: 12, color: colors.statusIntermediate, lineHeight: 1.5, marginBottom: 10 }}>
            {w.text}
          </div>
        ))}

        <div style={{ display: 'flex', flexDirection: 'column', gap: 3 }}>
          {data.conditions.map((c) => {
            const on = c.id === selected.id
            return (
              <div
                key={c.id}
                onClick={() => setPicked(c.id)}
                title={`${c.name} (${c.genes})`}
                style={{ display: 'grid', gridTemplateColumns: '92px 1fr 56px 74px', gap: 10, alignItems: 'center', padding: '3px 6px', borderRadius: 6, cursor: 'pointer', background: on ? colors.pageBg : 'transparent' }}
              >
                <span style={{ fontFamily: mono, fontSize: 12, fontWeight: c.matches ? 700 : 500 }}>{c.id}</span>
                <div style={{ position: 'relative', height: 12 }}>
                  <div style={{ position: 'absolute', left: `${x(0)}%`, top: 0, bottom: 0, width: 1, background: colors.axisLine }} />
                  <div style={{ position: 'absolute', left: `${x(data.threshold)}%`, top: -2, bottom: -2, borderLeft: `1.5px dashed ${colors.statusIntermediate}` }} />
                  <div
                    style={{
                      position: 'absolute',
                      top: 1,
                      height: 10,
                      borderRadius: 2,
                      background: c.matches ? EPI_CALL_COLOR.positive : colors.statusUnknown,
                      opacity: c.matches ? 1 : 0.55,
                      left: `${Math.min(x(0), x(c.score))}%`,
                      width: `${Math.abs(x(c.score) - x(0))}%`,
                    }}
                  />
                </div>
                <span style={{ textAlign: 'right', fontVariantNumeric: 'tabular-nums', fontSize: 12.5, fontWeight: c.matches ? 700 : 400 }}>{c.score.toFixed(2)}</span>
                <span style={{ fontSize: 11, color: c.n_found < c.n_cpgs ? colors.statusIntermediate : colors.textTertiary, fontVariantNumeric: 'tabular-nums' }}>
                  {c.n_found}/{c.n_cpgs} CpGs
                </span>
              </div>
            )
          })}
        </div>
        <div style={{ fontSize: 11.5, color: colors.textTertiary, marginTop: 10 }}>Dashed line: the {data.threshold} cut-off. Select a condition for its CpGs.</div>
      </Card>

      <Callout tone="caution" title="How far to trust these">
        Each classifier has a single positive example, a published median profile, so they were not tuned on individuals. In the preprint's long-read genomes they worked well for some conditions and poorly for others: the Kleefstra classifier
        called all 28 genomes of people with Kleefstra syndrome, the Kabuki classifier missed the two individuals with KDM6A variants, and the SETD1B classifier called only 1 of 3 PacBio genomes. A few controls scored above the cut-off on
        other conditions (Tatton-Brown-Rahman, SBBYSS). With 34 classifiers an occasional match is expected, and DNA from a cell line or another tissue can differ from blood for reasons unrelated to the condition. Treat a match as a lead.
        The cut-off here is 0.35, as in EpiSignaLR's code and methods (its results text says 0.3).
      </Callout>

      <ConditionDetail condition={selected} />
    </>
  )
}

const W = 900
const H = 250
const PAD = { left: 40, right: 14, top: 14, bottom: 30 }

// One condition's signature CpGs, ordered by how far the condition's reference profile sits from controls'. The
// reference profiles are the published array medians; the diamonds are this genome.
function ConditionDetail({ condition }: { condition: ArrayCondition }) {
  const { control, condition: ref, sample, coords } = condition.cpgs
  const order = coords
    .map((_, i) => i)
    .sort((a, b) => ((ref[a] ?? 0) - (control[a] ?? 0)) - ((ref[b] ?? 0) - (control[b] ?? 0)))
  const n = order.length
  const x = (rank: number) => PAD.left + ((rank + 0.5) / n) * (W - PAD.left - PAD.right)
  const y = (v: number) => PAD.top + (1 - v) * (H - PAD.top - PAD.bottom)
  const dot = (values: (number | null)[], fill: string, r: number, key: string) =>
    order.map((i, rank) => (values[i] == null ? null : <circle key={`${key}${rank}`} cx={x(rank)} cy={y(values[i] as number)} r={r} fill={fill} opacity={0.75} />))

  return (
    <Card>
      <CardHeading title={`${condition.id}: ${condition.name}`} sub={`${condition.genes} · score ${condition.score.toFixed(2)} · ${condition.n_found} of ${condition.n_cpgs} signature CpGs in this pileup`} />
      <svg viewBox={`0 0 ${W} ${H}`} style={{ width: '100%', height: 'auto', display: 'block' }} role="img" aria-label={`${condition.id} signature CpGs: reference controls, reference condition and this genome`}>
        {[0, 0.25, 0.5, 0.75, 1].map((t) => (
          <g key={t}>
            <line x1={PAD.left} x2={W - PAD.right} y1={y(t)} y2={y(t)} stroke={colors.gridLine} />
            <text x={PAD.left - 8} y={y(t) + 4} textAnchor="end" fontSize={11} fill={colors.textTertiary}>
              {Math.round(t * 100)}%
            </text>
          </g>
        ))}
        {dot(control, EPI_COLOR.control, 2.6, 'c')}
        {dot(ref, EPI_COLOR.case, 2.6, 'r')}
        {order.map((i, rank) =>
          sample[i] == null ? (
            <line key={`m${rank}`} x1={x(rank)} x2={x(rank)} y1={H - PAD.bottom} y2={H - PAD.bottom + 6} stroke={colors.statusIntermediate} strokeWidth={1.5}>
              <title>{`No pileup row at ${coords[i][0]}:${fmtPos(coords[i][1])}`}</title>
            </line>
          ) : (
            <polygon key={`s${rank}`} points={`${x(rank)},${y(sample[i] as number) - 4.5} ${x(rank) + 3.6},${y(sample[i] as number)} ${x(rank)},${y(sample[i] as number) + 4.5} ${x(rank) - 3.6},${y(sample[i] as number)}`} fill={EPI_COLOR.sample} stroke="white" strokeWidth={0.8}>
              <title>{`${coords[i][0]}:${fmtPos(coords[i][1])} · this genome ${((sample[i] as number) * 100).toFixed(1)}%`}</title>
            </polygon>
          ),
        )}
        <text x={PAD.left} y={H - 8} fontSize={10.5} fill={colors.textTertiary}>
          CpGs ordered by where the condition's reference profile sits relative to controls (left: furthest below)
        </text>
      </svg>
      <div style={{ display: 'flex', gap: 18, flexWrap: 'wrap', fontSize: 12, color: colors.textSecondary, marginTop: 8 }}>
        <Key color={EPI_COLOR.control} label="Reference controls (array median)" />
        <Key color={EPI_COLOR.case} label={`Reference ${condition.id} (array median)`} />
        <Key color={EPI_COLOR.sample} label="This genome" diamond />
        <span style={{ color: colors.statusIntermediate }}>| no pileup row</span>
      </div>
    </Card>
  )
}

function Key({ color, label, diamond }: { color: string; label: string; diamond?: boolean }) {
  return (
    <span style={{ display: 'inline-flex', alignItems: 'center', gap: 5 }}>
      <span style={{ width: 8, height: 8, background: color, borderRadius: diamond ? 0 : '50%', transform: diamond ? 'rotate(45deg)' : undefined }} />
      {label}
    </span>
  )
}
