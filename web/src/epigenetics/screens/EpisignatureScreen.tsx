import { Card, CardHeading } from '../../components/Card'
import { colors, mono } from '../../theme'
import type { EpiModel, EpisignatureOverview } from '../api'
import { CallChip, DotStrip, ModelFacts, ScoreBar, ScoreLegend } from '../episignatureParts'
import { EPI_COLOR, fmtPos, mean } from '../labels'
import { td, th } from '../styles'
import { Callout } from '../ui'

export function EpisignatureScreen({ data, modelId, onSelectModel }: { data: EpisignatureOverview; modelId: string; onSelectModel: (id: string) => void }) {
  const model = data.models.find((m) => m.id === modelId) ?? data.models[0]
  const ranked = [...model.features].sort((a, b) => Math.abs(b.contribution) - Math.abs(a.contribution))
  const reach = Math.max(0.01, ...ranked.map((f) => Math.abs(f.contribution)), Math.abs(model.intercept))

  return (
    <div style={{ padding: 28, display: 'flex', flexDirection: 'column', gap: 20 }}>
      <Card>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', gap: 20 }}>
          <div>
            <select
              value={model.id}
              onChange={(e) => onSelectModel(e.target.value)}
              style={{ fontSize: 16, fontWeight: 600, padding: '4px 8px', borderRadius: 6, border: `1px solid ${colors.borderInput}`, background: 'white' }}
            >
              {data.models.map((m) => (
                <option key={m.id} value={m.id}>
                  {m.name}
                </option>
              ))}
            </select>
            <div style={{ fontSize: 12, color: colors.textSecondary, marginTop: 8 }}>
              {model.gene}-anchored signature · {model.dmrs.selected} of {model.dmrs.panel} DMRs selected by the model · {model.dmrs.with_data} of {model.dmrs.panel} had data in this genome
            </div>
          </div>
          <div style={{ display: 'flex', alignItems: 'center', gap: 14 }}>
            <div style={{ fontSize: 30, fontWeight: 700, fontVariantNumeric: 'tabular-nums' }}>{model.probability == null ? '—' : model.probability.toFixed(2)}</div>
            <CallChip call={model.call} />
          </div>
        </div>
        <div style={{ maxWidth: 520, marginTop: 12 }}>
          <ScoreBar model={model} height={72} />
        </div>
        <div style={{ marginTop: 6 }}>
          <ScoreLegend />
        </div>
        <div style={{ marginTop: 10 }}>
          <ModelFacts model={model} />
        </div>
      </Card>

      {model.warnings.map((w) => (
        <Callout key={w.id} tone={w.severity === 'caution' ? 'caution' : 'note'} title={w.id.startsWith('pileup') ? 'Pileup mode' : 'Missing DMRs'}>
          {w.text}
        </Callout>
      ))}

      <Card>
        <CardHeading
          title="What drives the score"
          sub="The SVM's decision value is its intercept plus one term per selected DMR. Terms to the right push toward the case signature, terms to the left toward control."
        />
        <Contributions model={model} ranked={ranked} reach={reach} />
      </Card>

      <Card padded={false}>
        <div style={{ padding: '16px 20px', borderBottom: `1px solid ${colors.borderSoft}` }}>
          <div style={{ fontSize: 15, fontWeight: 600 }}>Selected DMRs</div>
          <div style={{ fontSize: 12, color: colors.textSecondary, marginTop: 2 }}>
            This genome (diamond) against the model's training {' '}
            <span style={{ color: EPI_COLOR.control, fontWeight: 600 }}>controls</span> and <span style={{ color: EPI_COLOR.case, fontWeight: 600 }}>cases</span>; each row has its own scale and the dashed line is the training mean.
          </div>
        </div>
        <div style={{ overflowX: 'auto' }}>
          <table style={{ width: '100%', borderCollapse: 'collapse' }}>
            <thead>
              <tr>
                <th style={th}>DMR</th>
                <th style={th}>Location</th>
                <th style={{ ...th, textAlign: 'right' }}>CpGs</th>
                <th style={{ ...th, textAlign: 'right' }}>This genome</th>
                <th style={{ ...th, textAlign: 'right', color: EPI_COLOR.control }}>Controls</th>
                <th style={{ ...th, textAlign: 'right', color: EPI_COLOR.case }}>Cases</th>
                <th style={th}>Distribution</th>
                <th style={{ ...th, textAlign: 'right' }}>Term</th>
              </tr>
            </thead>
            <tbody>
              {ranked.map((f) => {
                const c = mean(f.controls)
                const k = mean(f.cases)
                return (
                  <tr key={f.dmr} style={{ opacity: f.imputed ? 0.65 : 1 }}>
                    <td style={{ ...td, fontFamily: mono, fontSize: 12.5, fontWeight: 600 }}>{f.dmr}</td>
                    <td style={{ ...td, fontFamily: mono, fontSize: 11.5, color: colors.textSecondary }}>{f.region ? `${f.region.chrom}:${fmtPos(f.region.start)}` : '—'}</td>
                    <td style={{ ...td, textAlign: 'right', fontVariantNumeric: 'tabular-nums' }}>
                      {f.n_cpgs}/{f.n_cpgs_listed}
                    </td>
                    <td style={{ ...td, textAlign: 'right', fontVariantNumeric: 'tabular-nums', fontWeight: 600 }}>{f.value == null ? 'no data' : pct(f.value)}</td>
                    <td style={{ ...td, textAlign: 'right', fontVariantNumeric: 'tabular-nums' }}>{c == null ? '—' : pct(c)}</td>
                    <td style={{ ...td, textAlign: 'right', fontVariantNumeric: 'tabular-nums' }}>{k == null ? '—' : pct(k)}</td>
                    <td style={td}>
                      <DotStrip feature={f} />
                    </td>
                    <td style={{ ...td, textAlign: 'right', fontVariantNumeric: 'tabular-nums', color: f.contribution > 0 ? EPI_COLOR.case : EPI_COLOR.control }}>{signed(f.contribution)}</td>
                  </tr>
                )
              })}
            </tbody>
          </table>
        </div>
      </Card>

      <Card>
        <CardHeading title="How well the model did on its own data" sub={model.performance.design} />
        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(150px, 1fr))', gap: 16, fontSize: 13 }}>
          <Fact k="AUC" v={model.performance.auc.toFixed(2)} />
          <Fact k="Sensitivity" v={`${Math.round(model.performance.sensitivity * 100)}%`} />
          <Fact k="Specificity" v={`${Math.round(model.performance.specificity * 100)}%`} />
          <Fact k="Balanced accuracy" v={`${Math.round(model.performance.balanced_accuracy * 100)}%`} />
          <Fact k="Permutation p" v={model.performance.permutation_p.toFixed(3)} />
          <Fact k="Training set" v={`${model.training.cases} cases, ${model.training.controls} controls`} />
        </div>
        <div style={{ fontSize: 12, color: colors.textSecondary, marginTop: 14, lineHeight: 1.55 }}>
          With so few training genomes these are rough estimates. {model.performance.specificity < 1 ? `A specificity of ${Math.round(model.performance.specificity * 100)}% means that in the authors' own leave-one-out test some controls were called cases.` : ''}
        </div>
      </Card>
    </div>
  )
}

const pct = (v: number) => `${(v * 100).toFixed(1)}%`
const signed = (v: number) => `${v > 0 ? '+' : v < 0 ? '−' : ''}${Math.abs(v).toFixed(2)}`

function Fact({ k, v }: { k: string; v: string }) {
  return (
    <div>
      <div style={{ fontSize: 10.5, fontWeight: 600, color: colors.textSecondary, letterSpacing: '0.04em', textTransform: 'uppercase' }}>{k}</div>
      <div style={{ fontSize: 15, fontWeight: 600, marginTop: 3 }}>{v}</div>
    </div>
  )
}

// Diverging bars about zero: the intercept first, then each selected DMR's term, largest first.
function Contributions({ model, ranked, reach }: { model: EpiModel; ranked: EpiModel['features']; reach: number }) {
  const rows = [{ key: 'intercept', label: 'Intercept', value: model.intercept, note: 'baseline' }, ...ranked.map((f) => ({ key: f.dmr, label: f.dmr, value: f.contribution, note: f.imputed ? 'no data: training mean' : '' }))]
  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: 5 }}>
      {rows.map((r) => (
        <div key={r.key} style={{ display: 'grid', gridTemplateColumns: '110px 1fr 56px 150px', alignItems: 'center', gap: 10, fontSize: 12.5 }}>
          <span style={{ fontFamily: mono, fontSize: 12 }}>{r.label}</span>
          <div style={{ position: 'relative', height: 12 }}>
            <div style={{ position: 'absolute', left: '50%', top: 0, bottom: 0, width: 1, background: colors.axisLine }} />
            <div
              style={{
                position: 'absolute',
                top: 1,
                height: 10,
                borderRadius: 2,
                background: r.key === 'intercept' ? colors.statusUnknown : r.value > 0 ? EPI_COLOR.case : EPI_COLOR.control,
                left: r.value >= 0 ? '50%' : `${50 - (Math.abs(r.value) / reach) * 50}%`,
                width: `${(Math.abs(r.value) / reach) * 50}%`,
              }}
            />
          </div>
          <span style={{ textAlign: 'right', fontVariantNumeric: 'tabular-nums' }}>{signed(r.value)}</span>
          <span style={{ fontSize: 11.5, color: colors.textTertiary }}>{r.note}</span>
        </div>
      ))}
      <div style={{ display: 'grid', gridTemplateColumns: '110px 1fr 56px 150px', gap: 10, fontSize: 12.5, fontWeight: 600, borderTop: `1px solid ${colors.borderSoft}`, paddingTop: 8, marginTop: 4 }}>
        <span>Decision value</span>
        <span style={{ color: colors.textSecondary, fontWeight: 400 }}>sum of the above</span>
        <span style={{ textAlign: 'right', fontVariantNumeric: 'tabular-nums' }}>{signed(model.decision)}</span>
        <span />
      </div>
    </div>
  )
}
