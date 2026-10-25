import { useMemo, useState } from 'react'
import { Card, CardHeading } from '../../components/Card'
import { StatCard } from '../../components/StatCard'
import { colors, mono } from '../../theme'
import type { Category, DisorderSummary, ImprintingOverview, RegionResult } from '../api'
import { DISORDER_STATUS_COLOR, DISORDER_STATUS_LABEL, STATUS_COLOR, TRACK_COLOR, fmtPct } from '../labels'
import { td, th } from '../styles'
import { CategoryChip, Chip, MethBar, RegionStatusChip } from '../ui'

type CategoryFilter = 'all' | Category
type StatusFilter = 'all' | 'disturbed' | 'attention'

export function ImprintingScreen({ overview, onOpenRegion }: { overview: ImprintingOverview; onOpenRegion: (id: string) => void }) {
  const clinical = overview.regions.filter((r) => r.category === 'CA')
  const assessed = clinical.filter((r) => r.status !== 'no_data' && r.status !== 'low_coverage')
  const intact = clinical.filter((r) => r.status === 'imprinted')
  const disturbed = clinical.filter((r) => r.status === 'loss' || r.status === 'gain')
  const attention = clinical.filter((r) => r.status === 'borderline' || r.status === 'low_coverage' || r.status === 'no_data')

  return (
    <div style={{ padding: 28, display: 'flex', flexDirection: 'column', gap: 20 }}>
      <div style={{ display: 'grid', gridTemplateColumns: 'repeat(4,1fr)', gap: 16 }}>
        <StatCard label="Clinical DMRs assessed" value={`${assessed.length} / ${clinical.length}`} sub="of the diagnostic imprinting-disorder panel" />
        <StatCard label="Imprint intact" value={String(intact.length)} sub="one allele methylated, one not" dotColor={STATUS_COLOR.imprinted} />
        <StatCard
          label="Disturbed"
          value={String(disturbed.length)}
          sub="loss or gain of methylation"
          dotColor={disturbed.length ? STATUS_COLOR.loss : STATUS_COLOR.imprinted}
          valueStyle={disturbed.length ? { color: STATUS_COLOR.loss } : undefined}
        />
        <StatCard label="Needs a look" value={String(attention.length)} sub="borderline or too little coverage" dotColor={attention.length ? STATUS_COLOR.borderline : undefined} />
      </div>

      <Card>
        <CardHeading
          title="Imprinting disorders"
          sub="Patterns in the clinically associated DMRs that are characteristic of each disorder. A pattern is a lead, not a diagnosis."
        />
        <div style={{ display: 'flex', flexDirection: 'column' }}>
          {overview.disorders.map((d) => (
            <DisorderRow key={d.id} d={d} />
          ))}
          <MlidRow mlid={overview.mlid} />
        </div>
        <div style={{ fontSize: 11.5, color: colors.textTertiary, marginTop: 14, lineHeight: 1.5 }}>{overview.caveats}</div>
      </Card>

      <RegionTable regions={overview.regions} onOpenRegion={onOpenRegion} />
    </div>
  )
}

function DisorderRow({ d }: { d: DisorderSummary }) {
  const color = DISORDER_STATUS_COLOR[d.status]
  return (
    <div style={{ display: 'grid', gridTemplateColumns: '230px 130px 1fr', gap: 16, padding: '11px 0', borderTop: `1px solid ${colors.borderSoftest}`, alignItems: 'center' }}>
      <div>
        <div style={{ fontSize: 13.5, fontWeight: 600 }}>{d.name}</div>
        <div style={{ fontSize: 11.5, color: colors.textTertiary, fontFamily: mono }}>{d.locus}</div>
      </div>
      <div>
        <Chip color={color} solid={d.status === 'pattern_present'} title={d.meaning}>
          {DISORDER_STATUS_LABEL[d.status]}
        </Chip>
      </div>
      <div style={{ display: 'flex', flexWrap: 'wrap', gap: '4px 16px', fontSize: 12, color: colors.textBody }}>
        {d.triggers.map((t) => (
          <span key={t.dmr} title={t.primary ? 'Primary criterion' : 'Supportive criterion'}>
            <span style={{ fontFamily: mono }}>{t.name}</span>{' '}
            <span style={{ color: t.hit ? STATUS_COLOR.loss : colors.textTertiary }}>
              {fmtPct(t.methylation)} · {t.hit ? `${t.expected} ✓` : t.observed === 'not_assessed' ? 'not assessed' : `expected ${t.expected}`}
            </span>
            {!t.primary && <span style={{ color: colors.textTertiary }}> (supportive)</span>}
          </span>
        ))}
      </div>
    </div>
  )
}

function MlidRow({ mlid }: { mlid: ImprintingOverview['mlid'] }) {
  return (
    <div style={{ display: 'grid', gridTemplateColumns: '230px 130px 1fr', gap: 16, padding: '11px 0', borderTop: `1px solid ${colors.borderSoftest}`, alignItems: 'center' }}>
      <div>
        <div style={{ fontSize: 13.5, fontWeight: 600 }}>Multi-locus imprinting disturbance</div>
        <div style={{ fontSize: 11.5, color: colors.textTertiary }}>aberrant DMRs on 2+ chromosomes</div>
      </div>
      <div>
        <Chip color={mlid.possible ? STATUS_COLOR.loss : STATUS_COLOR.imprinted} solid={mlid.possible}>
          {mlid.possible ? 'Pattern present' : 'No pattern'}
        </Chip>
      </div>
      <div style={{ fontSize: 12, color: colors.textBody }}>
        {mlid.aberrant.length
          ? `${mlid.aberrant.length} aberrant on ${mlid.chromosomes.join(', ')}: ${mlid.aberrant.map((a) => a.name).join(', ')}`
          : `none of ${mlid.assessed} assessed ${mlid.scope} is aberrant`}
      </div>
    </div>
  )
}

function RegionTable({ regions, onOpenRegion }: { regions: RegionResult[]; onOpenRegion: (id: string) => void }) {
  const [category, setCategory] = useState<CategoryFilter>('all')
  const [statusFilter, setStatusFilter] = useState<StatusFilter>('all')

  const shown = useMemo(
    () =>
      regions.filter((r) => {
        if (category !== 'all' && r.category !== category) return false
        if (statusFilter === 'disturbed') return r.interpretable && (r.status === 'loss' || r.status === 'gain')
        if (statusFilter === 'attention') return r.interpretable && (r.status === 'borderline' || r.status === 'low_coverage' || r.status === 'no_data')
        return true
      }),
    [regions, category, statusFilter],
  )

  const select = (value: string, onChange: (v: string) => void, options: [string, string][]) => (
    <select
      value={value}
      onChange={(e) => onChange(e.target.value)}
      style={{ fontSize: 12.5, padding: '5px 8px', borderRadius: 6, border: `1px solid ${colors.borderInput}`, background: 'white', color: colors.textPrimary }}
    >
      {options.map(([v, l]) => (
        <option key={v} value={v}>
          {l}
        </option>
      ))}
    </select>
  )

  return (
    <Card padded={false}>
      <div style={{ padding: '16px 20px', display: 'flex', alignItems: 'center', justifyContent: 'space-between', borderBottom: `1px solid ${colors.borderSoft}` }}>
        <div>
          <div style={{ fontSize: 15, fontWeight: 600 }}>Imprinted regions</div>
          <div style={{ fontSize: 12, color: colors.textSecondary, marginTop: 2 }}>
            {shown.length} of {regions.length} DMRs · select a row for the CpG-level view
          </div>
        </div>
        <div style={{ display: 'flex', gap: 8 }}>
          {select(category, (v) => setCategory(v as CategoryFilter), [
            ['all', 'All categories'],
            ['CA', 'Clinical'],
            ['NC', 'Research'],
            ['X', 'Polymorphic'],
          ])}
          {select(statusFilter, (v) => setStatusFilter(v as StatusFilter), [
            ['all', 'Any status'],
            ['disturbed', 'Disturbed'],
            ['attention', 'Needs a look'],
          ])}
        </div>
      </div>
      <div style={{ overflowX: 'auto' }}>
        <table style={{ width: '100%', borderCollapse: 'collapse' }}>
          <thead>
            <tr>
              <th style={th}>DMR</th>
              <th style={th}>Category</th>
              <th style={th}>Combined</th>
              <th style={{ ...th, color: TRACK_COLOR.hap1 }}>Hap 1</th>
              <th style={{ ...th, color: TRACK_COLOR.hap2 }}>Hap 2</th>
              <th style={th}>Status</th>
              <th style={{ ...th, textAlign: 'right' }}>CpGs</th>
              <th style={{ ...th, textAlign: 'right' }}>Depth</th>
            </tr>
          </thead>
          <tbody>
            {shown.map((r) => (
              <tr
                key={r.id}
                onClick={() => onOpenRegion(r.id)}
                style={{ cursor: 'pointer', opacity: r.interpretable ? 1 : 0.7 }}
                onMouseEnter={(e) => (e.currentTarget.style.background = colors.pageBg)}
                onMouseLeave={(e) => (e.currentTarget.style.background = 'transparent')}
              >
                <td style={td}>
                  <div style={{ fontWeight: 600, fontFamily: mono, fontSize: 12.5 }}>{r.name}</div>
                  <div style={{ fontSize: 11.5, color: colors.textTertiary }}>
                    {r.chrom}
                    {r.disorders.length ? ` · ${r.disorders.join(', ')}` : ''}
                    {r.coordinate_notes ? ' · coordinates disputed' : ''}
                  </div>
                </td>
                <td style={td}>
                  <CategoryChip category={r.category} />
                </td>
                <td style={td}>
                  <MethBar value={r.combined?.methylation} />
                </td>
                <td style={td}>
                  <MethBar value={r.hap1?.methylation} color={TRACK_COLOR.hap1} width={56} />
                </td>
                <td style={td}>
                  <MethBar value={r.hap2?.methylation} color={TRACK_COLOR.hap2} width={56} />
                </td>
                <td style={td}>
                  <RegionStatusChip region={r} />
                </td>
                <td style={{ ...td, textAlign: 'right', fontVariantNumeric: 'tabular-nums' }}>{r.combined?.n_cpgs ?? '—'}</td>
                <td style={{ ...td, textAlign: 'right', fontVariantNumeric: 'tabular-nums' }}>{r.combined ? `${r.combined.mean_coverage.toFixed(0)}×` : '—'}</td>
              </tr>
            ))}
            {!shown.length && (
              <tr>
                <td colSpan={8} style={{ ...td, textAlign: 'center', color: colors.textSecondary, padding: 32 }}>
                  No regions match these filters.
                </td>
              </tr>
            )}
          </tbody>
        </table>
      </div>
    </Card>
  )
}
