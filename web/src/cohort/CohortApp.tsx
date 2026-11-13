import { useEffect, useMemo, useState } from 'react'
import type { Product } from '../components/ProductSwitcher'
import { ProductSwitcher } from '../components/ProductSwitcher'
import { Card, CardHeading } from '../components/Card'
import { StatCard, StatGrid } from '../components/StatCard'
import { StatusBadge, FlagBadge } from '../components/Badge'
import { colors, mono } from '../theme'
import { cohortApi, type CohortMode, type CohortResponse } from './api'

const MODES: { key: CohortMode; label: string }[] = [
  { key: 'carrier', label: 'Carrier' },
  { key: 'puretarget_v2', label: 'PureTarget V2' },
]

export default function CohortApp({
  product,
  onSwitchProduct,
}: {
  product: Product
  onSwitchProduct: (p: Product) => void
}) {
  const [mode, setMode] = useState<CohortMode>('carrier')
  const [data, setData] = useState<CohortResponse | null>(null)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    let current = true
    setData(null)
    setError(null)
    cohortApi
      .get(mode)
      .then((d) => current && setData(d))
      .catch((e) => current && setError(String(e)))
    return () => {
      current = false
    }
  }, [mode])

  return (
    <div style={{ display: 'flex', width: '100%', minHeight: '100vh', background: colors.pageBg, color: colors.textPrimary }}>
      <div style={{ width: 232, flexShrink: 0, background: colors.sidebarBg, display: 'flex', flexDirection: 'column', padding: '20px 14px', gap: 4 }}>
        <ProductSwitcher product={product} onChange={onSwitchProduct} />
        <div style={{ color: 'white', fontWeight: 600, fontSize: 15, padding: '8px 8px 16px' }}>Cohort</div>
        <div style={{ color: colors.sidebarTextMuted, fontSize: 11, fontWeight: 600, letterSpacing: '0.08em', textTransform: 'uppercase', padding: '4px 10px 6px' }}>
          Assay
        </div>
        {MODES.map((m) => (
          <div
            key={m.key}
            onClick={() => setMode(m.key)}
            style={{
              padding: '10px 12px',
              borderRadius: 8,
              cursor: 'pointer',
              fontSize: 13.5,
              fontWeight: 500,
              color: colors.sidebarTextItem,
              background: mode === m.key ? colors.sidebarHover : 'transparent',
            }}
          >
            {m.label}
          </div>
        ))}
      </div>

      <div style={{ flex: 1, display: 'flex', flexDirection: 'column', minWidth: 0 }}>
        <div style={{ padding: '16px 28px', borderBottom: `1px solid ${colors.border}`, background: 'white' }}>
          <div style={{ fontSize: 19, fontWeight: 600 }}>Cohort Overview</div>
          <div style={{ fontSize: 12.5, color: colors.textSecondary, marginTop: 2 }}>
            {data ? `${data.label} · ${data.stats.samples} samples` : MODES.find((m) => m.key === mode)?.label}
          </div>
        </div>
        <div style={{ flex: 1, overflowY: 'auto', padding: 28 }}>
          {error && <Card>Failed to load cohort: {error}</Card>}
          {!error && !data && <Card>Loading cohort…</Card>}
          {data && <CohortOverview data={data} />}
        </div>
      </div>
    </div>
  )
}

function CohortOverview({ data }: { data: CohortResponse }) {
  const { stats } = data
  const [onlyReportable, setOnlyReportable] = useState(true)
  const rows = useMemo(
    () => data.samples.filter((s) => !onlyReportable || s.reportable.length > 0),
    [data, onlyReportable],
  )
  const maxBin = Math.max(1, ...data.depth_histogram.map((b) => b.count))
  const pct = (n: number) => (stats.samples ? `${Math.round((100 * n) / stats.samples)}%` : '—')

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: 20 }}>
      <StatGrid>
        <StatCard label="Samples" value={String(stats.samples)} sub={`${data.runs.length} run${data.runs.length === 1 ? '' : 's'}`} />
        <StatCard
          label="QC pass"
          value={`${stats.qc_pass} / ${stats.samples}`}
          sub={`All loci ≥ ${data.min_locus_depth}x · ${stats.qc_fail} flagged`}
          dotColor={stats.qc_fail ? colors.statusIntermediate : colors.statusNormal}
        />
        <StatCard
          label="Reportable samples"
          value={String(stats.samples_with_reportable)}
          sub={`${pct(stats.samples_with_reportable)} of cohort · ${stats.reportable_findings} findings`}
          dotColor={stats.samples_with_reportable ? colors.statusPathogenic : colors.statusNormal}
        />
        <StatCard
          label="Mean locus depth"
          value={stats.mean_depth != null ? `${stats.mean_depth}x` : '—'}
          sub={stats.median_depth != null ? `median ${stats.median_depth}x per sample` : undefined}
        />
      </StatGrid>

      <div style={{ display: 'grid', gridTemplateColumns: '2fr 1fr', gap: 20 }}>
        <Card>
          <CardHeading title="Coverage distribution" sub="Mean depth across genotyped loci, samples per bin" />
          <div style={{ display: 'flex', alignItems: 'flex-end', gap: 6, height: 140 }}>
            {data.depth_histogram.map((b) => (
              <div key={b.label} style={{ flex: 1, display: 'flex', flexDirection: 'column', alignItems: 'center', justifyContent: 'flex-end', height: '100%' }}>
                <div style={{ fontSize: 11, color: colors.textSecondary }}>{b.count || ''}</div>
                <div
                  style={{
                    width: '100%',
                    height: `${(b.count / maxBin) * 100}%`,
                    background: b.lo < data.min_locus_depth ? colors.statusIntermediate : colors.accentButtonBg,
                    borderRadius: '3px 3px 0 0',
                    minHeight: b.count ? 2 : 0,
                  }}
                />
                <div style={{ fontSize: 10.5, color: colors.textTertiary, marginTop: 4, whiteSpace: 'nowrap' }}>{b.label}</div>
              </div>
            ))}
          </div>
        </Card>
        <Card>
          <CardHeading title="Composition" />
          <Breakdown title="Sex" items={data.sex.map((s) => ({ label: s.label === 'M' ? 'Male' : s.label === 'F' ? 'Female' : s.label, count: s.count }))} total={stats.samples} />
          <Breakdown title="Run" items={data.runs.map((r) => ({ label: r.run, count: r.count }))} total={stats.samples} />
        </Card>
      </div>

      <Card padded={false}>
        <div style={{ padding: '18px 24px 8px', display: 'flex', justifyContent: 'space-between', alignItems: 'baseline' }}>
          <div>
            <CardHeading title="Reportable findings" sub={`${data.genes.length} gene${data.genes.length === 1 ? '' : 's'} with findings`} />
          </div>
          <label style={{ fontSize: 12.5, color: colors.textSecondary, cursor: 'pointer' }}>
            <input type="checkbox" checked={onlyReportable} onChange={(e) => setOnlyReportable(e.target.checked)} /> Only samples with findings
          </label>
        </div>
        {data.genes.length > 0 && (
          <div style={{ padding: '0 24px 12px', display: 'flex', flexWrap: 'wrap', gap: 8 }}>
            {data.genes.map((g) => (
              <span key={g.gene} title={g.disease} style={{ background: colors.chipBg, borderRadius: 20, padding: '4px 10px', fontSize: 12, fontWeight: 600 }}>
                {g.gene} <span style={{ color: colors.textSecondary, fontWeight: 500 }}>· {g.samples}</span>
              </span>
            ))}
          </div>
        )}
        <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: 13 }}>
          <thead>
            <tr style={{ textAlign: 'left', color: colors.textSecondary, fontSize: 11.5 }}>
              {['Sample', 'QC', 'Depth', 'Sex', 'Findings'].map((h) => (
                <th key={h} style={{ padding: '8px 24px', borderTop: `1px solid ${colors.borderSoft}`, fontWeight: 600 }}>{h}</th>
              ))}
            </tr>
          </thead>
          <tbody>
            {rows.length === 0 && (
              <tr><td colSpan={5} style={{ padding: '24px', textAlign: 'center', color: colors.textSecondary }}>No samples to show.</td></tr>
            )}
            {rows.map((s) => (
              <tr key={s.sample} style={{ borderTop: `1px solid ${colors.borderSofter}` }}>
                <td style={{ padding: '10px 24px', fontFamily: mono, fontSize: 12 }}>{s.sample}</td>
                <td style={{ padding: '10px 24px' }}><FlagBadge flag={s.qc === 'No data' ? 'No call' : s.qc === 'Pass' ? 'Pass' : 'Low depth'} /></td>
                <td style={{ padding: '10px 24px' }}>{s.mean_depth != null ? `${s.mean_depth}x` : '—'}</td>
                <td style={{ padding: '10px 24px' }}>{s.sex ?? '—'}</td>
                <td style={{ padding: '10px 24px' }}>
                  {s.reportable.length === 0 ? (
                    <span style={{ color: colors.textTertiary }}>None</span>
                  ) : (
                    <div style={{ display: 'flex', flexDirection: 'column', gap: 4 }}>
                      {s.reportable.map((f, i) => (
                        <div key={i} style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
                          <b>{f.gene}</b>
                          <StatusBadge status={f.status} />
                          {f.genotype && <span style={{ fontFamily: mono, fontSize: 11.5, color: colors.textSecondary }}>{f.genotype}</span>}
                        </div>
                      ))}
                    </div>
                  )}
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </Card>
    </div>
  )
}

function Breakdown({ title, items, total }: { title: string; items: { label: string; count: number }[]; total: number }) {
  return (
    <div style={{ marginBottom: 14 }}>
      <div style={{ fontSize: 11, fontWeight: 600, color: colors.textSecondary, textTransform: 'uppercase', letterSpacing: '0.06em', marginBottom: 4 }}>{title}</div>
      {items.slice(0, 6).map((i) => (
        <div key={i.label} style={{ display: 'flex', justifyContent: 'space-between', fontSize: 12.5, padding: '2px 0' }}>
          <span style={{ overflow: 'hidden', textOverflow: 'ellipsis' }}>{i.label}</span>
          <span style={{ color: colors.textSecondary }}>{i.count} ({total ? Math.round((100 * i.count) / total) : 0}%)</span>
        </div>
      ))}
      {items.length > 6 && <div style={{ fontSize: 11.5, color: colors.textTertiary }}>+{items.length - 6} more</div>}
    </div>
  )
}
