import { Card, CardHeading } from '../../components/Card'
import { colors, mono } from '../../theme'
import type { EpiModel, EpiWarning, EpisignatureResponse } from '../api'
import { CallChip, ModelFacts, ScoreBar, ScoreLegend } from '../episignatureParts'
import { Callout, Message } from '../ui'

export function EpisignaturesScreen({
  data,
  error,
  onOpenModel,
}: {
  data: EpisignatureResponse | null
  error: string | null
  onOpenModel: (id: string) => void
}) {
  if (error && !data) return <Message tone="error">Episignatures could not be scored: {error}</Message>
  if (!data) return <Message>Scoring episignatures…</Message>

  if (!data.installed) {
    return (
      <div style={{ padding: 28 }}>
        <Card>
          <CardHeading
            title="Episignature models are not set up"
            sub="They come from the EpiSignaLR project, which publishes no licence, so they are fetched on this machine rather than shipped with the app."
          />
          <div style={{ fontSize: 13, marginBottom: 10 }}>Run this once from the repository root, then reload:</div>
          <code style={{ display: 'block', padding: '10px 14px', background: colors.pageBg, borderRadius: 8, fontFamily: mono, fontSize: 12.5 }}>{data.setup_command}</code>
          <div style={{ fontSize: 12.5, color: colors.textSecondary, marginTop: 16, lineHeight: 1.6 }}>
            It then scores a genome for {data.available.map((a) => `${a.name} (${a.gene})`).join(', ')}.
          </div>
        </Card>
      </div>
    )
  }

  const cases = data.models.map((m) => m.training.cases)
  const caseRange = Math.min(...cases) === Math.max(...cases) ? String(cases[0]) : `${Math.min(...cases)} to ${Math.max(...cases)}`

  // The pileup-mode warning is the same on every model, so show it once at the top.
  const shared = new Map<string, EpiWarning>()
  data.models.forEach((m) => m.warnings.filter((w) => w.id.startsWith('pileup')).forEach((w) => shared.set(w.id, w)))

  return (
    <div style={{ padding: 28, display: 'flex', flexDirection: 'column', gap: 20 }}>
      <Callout tone="caution" title="Read these scores as leads, not results">
        The models were trained on nanopore genomes of blood DNA from a few cases ({caseRange} per condition). PacBio genomes sit nearer the case group than nanopore controls do, and in the authors' own data PacBio libraries made most of the
        indeterminate and false-positive calls. DNA from other sources, such as a cell line, can differ for reasons unrelated to the condition. A match or an indeterminate score is a reason to look further; a negative score does not rule a
        condition out. EpiSignaLR is a preprint, for research use only.
      </Callout>
      {[...shared.values()].map((w) => (
        <Callout key={w.id} tone="caution" title="This pileup is not the kind the models were built on">
          {w.text}
        </Callout>
      ))}

      <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(300px, 1fr))', gap: 16 }}>
        {data.models.map((m) => (
          <ModelCard key={m.id} model={m} onOpen={() => onOpenModel(m.id)} />
        ))}
      </div>
      <ScoreLegend />

      <Card>
        <CardHeading title="About these scores" />
        <div style={{ fontSize: 12.5, color: colors.textBody, lineHeight: 1.65 }}>
          <p style={{ margin: '0 0 8px' }}>
            For each condition the genome's methylation is averaged over a panel of DMRs (the coverage-weighted fraction of modified calls across the DMR's CpGs), then a linear SVM scores the handful of DMRs the model selected.
            A probability of 0.5 or more is a match and 0.4 to 0.6 is flagged indeterminate. A DMR with no covered CpG is replaced by the training mean, and no call is made if that applies to more than a quarter of the selected DMRs.
          </p>
          <p style={{ margin: '0 0 8px' }}>
            The scoring is reimplemented here and matches the authors' own scripts: on the same pileup the DMR values are identical and the probabilities agree to four decimal places.
          </p>
          <p style={{ margin: 0, color: colors.textSecondary }}>
            Source: {data.source.repo} at commit <span style={{ fontFamily: mono }}>{data.source.commit.slice(0, 7)}</span>, fetched {data.source.fetched_at.slice(0, 10)}. Licence: {data.source.licence}.
          </p>
        </div>
      </Card>
    </div>
  )
}

function ModelCard({ model, onOpen }: { model: EpiModel; onOpen: () => void }) {
  const own = model.warnings.filter((w) => !w.id.startsWith('pileup'))
  return (
    <Card style={{ display: 'flex', flexDirection: 'column', gap: 10 }}>
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', gap: 10 }}>
        <div>
          <div style={{ fontSize: 14.5, fontWeight: 600 }}>{model.name}</div>
          <div style={{ fontSize: 11.5, color: colors.textTertiary, fontFamily: mono }}>{model.gene}</div>
        </div>
        <CallChip call={model.call} />
      </div>
      <div style={{ fontSize: 30, fontWeight: 700, fontVariantNumeric: 'tabular-nums' }}>
        {model.probability == null ? '—' : model.probability.toFixed(2)}
        <span style={{ fontSize: 12, fontWeight: 500, color: colors.textTertiary, marginLeft: 8 }}>case probability</span>
      </div>
      <ScoreBar model={model} />
      <ModelFacts model={model} />
      {own.map((w) => (
        <div key={w.id} style={{ fontSize: 12, color: colors.statusIntermediate, lineHeight: 1.5 }}>
          {w.text}
        </div>
      ))}
      <button
        onClick={onOpen}
        style={{ alignSelf: 'flex-start', marginTop: 2, padding: '6px 14px', borderRadius: 6, border: `1px solid ${colors.borderInput}`, background: 'white', fontSize: 12.5, fontWeight: 600, cursor: 'pointer' }}
      >
        See the DMRs behind this score
      </button>
    </Card>
  )
}
