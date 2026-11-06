import type { ExpressionEffects, ScEtrEffect, ScEtrLead, ScEtrTarget } from '../api'
import { colors, mono } from '../theme'
import { Card, CardHeading } from './Card'

// Direction is carried by the glyph as well as the color, and the colors are deliberately not the
// status palette: red/green already mean pathogenic/normal elsewhere in this app, and an expression
// change is neither.
const UP = colors.cluster1
const DOWN = colors.methylated

const TITLE = 'Expression effects in immune cells'
const SUB = 'TenK10K · common repeat-length variation vs. gene expression · research annotation, not a pathogenicity call'

function cellName(cellType: string): string {
  return cellType.replace(/_/g, ' ')
}

function fmtQ(q: number): string {
  return q === 0 ? '<1e-300' : q.toExponential(1)
}

function effectTitle(e: ScEtrEffect): string {
  return `${cellName(e.cell_type)}: ${e.coeff > 0 ? '+' : ''}${e.coeff} per extra repeat unit (SE ${e.se}), gene-level q ${fmtQ(e.qval)}`
}

function directionSummary(effects: ScEtrEffect[]): { text: string; color: string } {
  const up = effects.filter((e) => e.coeff > 0).length
  const down = effects.length - up
  if (down === 0) return { text: '▲ higher expression with longer repeats', color: UP }
  if (up === 0) return { text: '▼ lower expression with longer repeats', color: DOWN }
  return { text: `▲ ${up} · ▼ ${down} across cell types`, color: colors.textSecondary }
}

export function ExpressionEffectsCard({ data, motif }: { data: ExpressionEffects; motif: string }) {
  if (!data.available) {
    return (
      <Card>
        <CardHeading title={TITLE} sub={SUB} />
        <div style={{ fontSize: 13, color: colors.statusUncertain }}>
          The expression annotation file could not be read: {data.error}
        </div>
      </Card>
    )
  }

  const atRepeat = data.leads.filter((l) => l.relation === 'overlaps')
  const nearby = data.leads.filter((l) => l.relation === 'adjacent')
  const paper = data.source.paper

  return (
    <Card>
      <CardHeading title={TITLE} sub={SUB} />

      {data.leads.length === 0 && (
        <div style={{ fontSize: 13, color: colors.textBody, lineHeight: 1.55 }}>
          {data.polymorphic
            ? 'No gene’s strongest expression association in this study is a repeat at this locus. The repeat was variable in the cohort, but not being a lead association is not evidence of no effect.'
            : 'No repeat at this locus is in the study’s polymorphic catalog (repeats found variable in a ~200-person subset), so it was probably not tested.'}
        </div>
      )}

      {atRepeat.length > 0 && (
        <Group label="At this repeat">
          {atRepeat.map((l) => (
            <LeadBlock key={`${l.start}-${l.motif}`} lead={l} diseaseMotif={motif} />
          ))}
        </Group>
      )}

      {nearby.length > 0 && (
        <Group label="Nearby repeat · a different repeat, not the disease repeat" muted>
          {nearby.map((l) => (
            <LeadBlock key={`${l.start}-${l.motif}`} lead={l} diseaseMotif={motif} muted />
          ))}
        </Group>
      )}

      <div style={{ marginTop: 16, paddingTop: 12, borderTop: `1px solid ${colors.borderSofter}`, fontSize: 11.5, lineHeight: 1.55, color: colors.textTertiary }}>
        Effect is the change in normalized expression per extra repeat unit (both alleles summed), in common-range repeat
        lengths in {paper.cohort}; {data.source.significance}. Only each gene&rsquo;s lead repeat per cell type is released,
        so a repeat that is not listed has not been shown to have no effect. Source: {paper.authors},{' '}
        <a href={`https://doi.org/${paper.doi}`} target="_blank" rel="noreferrer" style={{ color: colors.accentButtonBg }}>
          bioRxiv {paper.doi}
        </a>{' '}
        ({paper.version}); data{' '}
        <a href={`https://doi.org/${paper.data_doi}`} target="_blank" rel="noreferrer" style={{ color: colors.accentButtonBg }}>
          {paper.data_doi}
        </a>
        ; {paper.license}.
      </div>
    </Card>
  )
}

function Group({ label, muted = false, children }: { label: string; muted?: boolean; children: React.ReactNode }) {
  return (
    <div style={{ marginTop: 6 }}>
      <div
        style={{
          fontSize: 11.5,
          fontWeight: 600,
          textTransform: 'uppercase',
          letterSpacing: '0.04em',
          color: muted ? colors.textTertiary : colors.textSecondary,
          margin: '10px 0 8px',
        }}
      >
        {label}
      </div>
      <div style={{ display: 'flex', flexDirection: 'column', gap: 12 }}>{children}</div>
    </div>
  )
}

function LeadBlock({ lead, diseaseMotif, muted = false }: { lead: ScEtrLead; diseaseMotif: string; muted?: boolean }) {
  const place = lead.relation === 'adjacent' ? (lead.distance_bp === 0 ? 'adjacent' : `${lead.distance_bp.toLocaleString()} bp away`) : null
  return (
    <div
      style={{
        border: `1px solid ${colors.borderSofter}`,
        borderRadius: 8,
        padding: '10px 14px',
        background: muted ? colors.borderSoftest : 'white',
      }}
    >
      <div style={{ display: 'flex', alignItems: 'baseline', gap: 10, flexWrap: 'wrap', fontSize: 12.5 }}>
        <span style={{ fontFamily: mono, color: colors.textSecondary }}>
          {lead.chrom}:{lead.start.toLocaleString()}-{lead.stop.toLocaleString()}
        </span>
        <span style={{ fontFamily: mono, fontWeight: 700, color: colors.textPrimary }}>({lead.motif})n</span>
        {place && <span style={{ color: colors.textTertiary }}>{place}</span>}
        {lead.relation === 'overlaps' && !lead.same_motif && (
          <span style={{ color: colors.statusUncertain }}>
            motif differs from the disease motif {diseaseMotif}: a different repeat unit in the same tract
          </span>
        )}
      </div>
      <div style={{ display: 'flex', flexDirection: 'column', gap: 10, marginTop: 8 }}>
        {lead.targets.map((t) => (
          <TargetRow key={t.gene_id} target={t} />
        ))}
      </div>
    </div>
  )
}

function TargetRow({ target }: { target: ScEtrTarget }) {
  const summary = directionSummary(target.effects)
  const bestQ = Math.min(...target.effects.map((e) => e.qval))
  const n = target.effects.length
  return (
    <div>
      <div style={{ display: 'flex', alignItems: 'baseline', gap: 10, flexWrap: 'wrap', fontSize: 13 }}>
        <span style={{ fontWeight: 700 }} title={target.gene_id}>
          {target.gene ?? target.gene_id}
        </span>
        <span style={{ color: colors.textSecondary }}>
          {n} cell type{n === 1 ? '' : 's'}
        </span>
        <span style={{ color: summary.color, fontWeight: 600 }}>{summary.text}</span>
        <span style={{ color: colors.textTertiary, fontFamily: mono, fontSize: 12 }}>best q {fmtQ(bestQ)}</span>
      </div>
      <div style={{ display: 'flex', flexWrap: 'wrap', gap: 6, marginTop: 6 }}>
        {target.effects.map((e) => {
          const color = e.coeff > 0 ? UP : DOWN
          return (
            <span
              key={e.cell_type}
              title={effectTitle(e)}
              style={{
                fontSize: 11.5,
                padding: '2px 8px',
                borderRadius: 10,
                background: colors.chipBg,
                color: colors.textBody,
                border: `1px solid ${colors.borderSoft}`,
                whiteSpace: 'nowrap',
              }}
            >
              {cellName(e.cell_type)}{' '}
              <span style={{ color, fontWeight: 600, fontFamily: mono }}>
                {e.coeff > 0 ? '▲' : '▼'} {Math.abs(e.coeff).toFixed(3)}
              </span>
            </span>
          )
        })}
      </div>
    </div>
  )
}
