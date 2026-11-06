import { colors } from '../theme'
import type { EpiCall, EpiFeature, EpiModel } from './api'
import { EPI_CALL_COLOR, EPI_CALL_LABEL, EPI_COLOR } from './labels'
import { Chip } from './ui'

export function CallChip({ call }: { call: EpiCall }) {
  return (
    <Chip color={EPI_CALL_COLOR[call]} solid={call === 'positive' || call === 'negative'}>
      {EPI_CALL_LABEL[call]}
    </Chip>
  )
}

const BAR_W = 300
const BAR_PAD = 8

// The case probability on a 0-1 axis: the model's indeterminate zone shaded, its threshold marked, the authors'
// leave-one-out scores for their own training samples as small dots above, and this genome as the diamond.
export function ScoreBar({ model, height = 64 }: { model: EpiModel; height?: number }) {
  const x = (p: number) => BAR_PAD + p * (BAR_W - 2 * BAR_PAD)
  const [lo, hi] = model.indeterminate_zone
  const trackY = 36
  const p = model.probability
  return (
    <svg
      viewBox={`0 0 ${BAR_W} 64`}
      style={{ width: '100%', height, display: 'block' }}
      role="img"
      aria-label={`Case probability ${p == null ? 'not available' : p.toFixed(2)} on a 0 to 1 scale, indeterminate between ${lo} and ${hi}`}
    >
      <rect x={x(0)} y={trackY} width={x(lo) - x(0)} height={10} fill="oklch(0.94 0.03 230)" />
      <rect x={x(lo)} y={trackY} width={x(hi) - x(lo)} height={10} fill="oklch(0.94 0.05 80)" />
      <rect x={x(hi)} y={trackY} width={x(1) - x(hi)} height={10} fill="oklch(0.94 0.04 25)" />
      <line x1={x(model.threshold)} x2={x(model.threshold)} y1={trackY - 3} y2={trackY + 13} stroke={colors.axisLine} strokeWidth={1.5} />
      {model.loo.map((s, i) => (
        <circle key={i} cx={x(s.score)} cy={s.group ? 24 : 15} r={2.6} fill={s.group ? EPI_COLOR.case : EPI_COLOR.control} opacity={0.6}>
          <title>{`Training ${s.group ? 'case' : 'control'}, leave-one-out score ${s.score.toFixed(2)}`}</title>
        </circle>
      ))}
      {p != null && <polygon points={`${x(p)},${trackY - 1} ${x(p) + 6},${trackY + 5} ${x(p)},${trackY + 11} ${x(p) - 6},${trackY + 5}`} fill={EPI_COLOR.sample} stroke="white" strokeWidth={1.2} />}
      {[0, 0.5, 1].map((t) => (
        <text key={t} x={x(t)} y={61} textAnchor={t === 0 ? 'start' : t === 1 ? 'end' : 'middle'} fontSize={9.5} fill={colors.textTertiary}>
          {t}
        </text>
      ))}
    </svg>
  )
}

export function ScoreLegend() {
  const dot = (color: string, label: string) => (
    <span style={{ display: 'inline-flex', alignItems: 'center', gap: 5 }}>
      <span style={{ width: 8, height: 8, borderRadius: '50%', background: color, opacity: 0.7 }} />
      {label}
    </span>
  )
  return (
    <div style={{ display: 'flex', gap: 14, flexWrap: 'wrap', fontSize: 11.5, color: colors.textSecondary }}>
      <span style={{ display: 'inline-flex', alignItems: 'center', gap: 5 }}>
        <span style={{ width: 8, height: 8, background: EPI_COLOR.sample, transform: 'rotate(45deg)' }} />
        This genome
      </span>
      {dot(EPI_COLOR.control, 'Training controls (leave-one-out)')}
      {dot(EPI_COLOR.case, 'Training cases (leave-one-out)')}
      <span>Shaded band: indeterminate</span>
    </div>
  )
}

const STRIP_W = 200

// One DMR's methylation in the training controls and cases, with this genome's value. The scale is fitted to the
// values shown, so a DMR whose controls and cases differ by a few points still separates visibly.
export function DotStrip({ feature }: { feature: EpiFeature }) {
  const all = [...feature.controls, ...feature.cases, ...(feature.value == null ? [] : [feature.value])]
  const lo = Math.min(...all)
  const hi = Math.max(...all)
  const pad = Math.max((hi - lo) * 0.08, 0.005)
  const x = (v: number) => 4 + ((v - (lo - pad)) / (hi - lo + 2 * pad)) * (STRIP_W - 8)
  return (
    <svg viewBox={`0 0 ${STRIP_W} 26`} style={{ width: STRIP_W, height: 26, display: 'block' }} role="img" aria-label={`${feature.dmr} methylation, training controls and cases`}>
      <line x1={4} x2={STRIP_W - 4} y1={13} y2={13} stroke={colors.gridLine} />
      <line x1={x(feature.training_mean)} x2={x(feature.training_mean)} y1={4} y2={22} stroke={colors.axisLine} strokeDasharray="2 2" />
      {feature.controls.map((v, i) => (
        <circle key={`c${i}`} cx={x(v)} cy={9} r={2.4} fill={EPI_COLOR.control} opacity={0.6} />
      ))}
      {feature.cases.map((v, i) => (
        <circle key={`k${i}`} cx={x(v)} cy={17} r={2.4} fill={EPI_COLOR.case} opacity={0.6} />
      ))}
      {feature.value != null && (
        <polygon points={`${x(feature.value)},5 ${x(feature.value) + 5},13 ${x(feature.value)},21 ${x(feature.value) - 5},13`} fill={EPI_COLOR.sample} stroke="white" strokeWidth={1} />
      )}
    </svg>
  )
}

export function ModelFacts({ model }: { model: EpiModel }) {
  const p = model.performance
  const pct = (v: number) => `${Math.round(v * 100)}%`
  return (
    <div style={{ fontSize: 12, color: colors.textSecondary, lineHeight: 1.55 }}>
      Trained on {model.training.cases} cases and {model.training.controls} controls. Leave-one-out: AUC {p.auc.toFixed(2)}, sensitivity {pct(p.sensitivity)}, specificity {pct(p.specificity)}.
    </div>
  )
}
