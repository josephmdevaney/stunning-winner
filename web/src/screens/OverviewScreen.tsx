import { useEffect, useState } from 'react'
import { api, NON_REPEAT_ID_PREFIX, type OverviewResponse } from '../api'
import { colors, mono } from '../theme'
import { StatCard, StatGrid } from '../components/StatCard'
import { StatusBadge } from '../components/Badge'

const PILLS = ['All loci', 'Pathogenic panel', 'Novel candidates']

export function OverviewScreen({
  mode, sample, onSelectLocus, onSelectGeneFinding,
}: {
  mode: string; sample: string; onSelectLocus: (id: string) => void; onSelectGeneFinding: (gene: string) => void
}) {
  const [data, setData] = useState<OverviewResponse | null>(null)

  useEffect(() => {
    setData(null)
    api.overview(mode, sample).then(setData)
  }, [mode, sample])

  if (!data) return <div style={{ padding: 28 }}>Loading…</div>

  const { stats, loci } = data

  return (
    <div style={{ padding: 28, display: 'flex', flexDirection: 'column', gap: 22 }}>
      <StatGrid>
        <StatCard label="Loci analyzed" value={String(stats.loci_analyzed)} sub="Disease-associated panel" />
        <StatCard
          label="Pathogenic findings"
          value={String(stats.pathogenic_findings)}
          sub="Above clinical threshold"
          valueStyle={{ color: colors.statusPathogenic }}
        />
        <StatCard
          label="Mean read depth"
          value={stats.mean_read_depth != null ? `${stats.mean_read_depth}x` : '—'}
          sub="HiFi consensus reads, genotyped loci"
        />
        <StatCard
          label="Sample QC"
          value={stats.qc_status}
          sub={`${stats.loci_genotyped}/${stats.loci_analyzed} loci genotyped`}
          valueStyle={{ color: stats.qc_status === 'Pass' ? colors.statusNormal : colors.statusIntermediate }}
        />
      </StatGrid>

      <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, overflow: 'hidden' }}>
        <div
          style={{
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'space-between',
            padding: '16px 20px',
            borderBottom: `1px solid ${colors.borderSoft}`,
          }}
        >
          <div style={{ fontSize: 15, fontWeight: 600 }}>Repeat expansion panel</div>
          <div style={{ display: 'flex', gap: 8 }}>
            {PILLS.map((p, i) => (
              <div
                key={p}
                style={{
                  padding: '6px 12px',
                  borderRadius: 20,
                  background: i === 0 ? colors.accentButtonBg : colors.chipBg,
                  color: i === 0 ? 'white' : colors.textLabel,
                  fontSize: 12,
                  fontWeight: 600,
                }}
              >
                {p}
              </div>
            ))}
          </div>
        </div>

        <div
          style={{
            display: 'grid',
            gridTemplateColumns: '1fr 1.6fr 0.8fr 0.7fr 0.7fr 1fr 0.6fr 1.3fr',
            gap: 12,
            padding: '10px 20px',
            fontSize: 11,
            fontWeight: 600,
            color: colors.textSecondary,
            letterSpacing: '0.04em',
            textTransform: 'uppercase',
            borderBottom: `1px solid ${colors.borderSofter}`,
          }}
        >
          <div>Gene</div>
          <div>Region (GRCh38)</div>
          <div>Motif</div>
          <div>Allele 1</div>
          <div>Allele 2</div>
          <div>Zygosity</div>
          <div>Depth</div>
          <div>Status</div>
        </div>

        {loci.map((l) => (
          <div
            key={l.id}
            onClick={
              l.clickable
                ? () => (l.id.startsWith(NON_REPEAT_ID_PREFIX) ? onSelectGeneFinding(l.gene) : onSelectLocus(l.id))
                : undefined
            }
            style={{
              display: 'grid',
              gridTemplateColumns: '1fr 1.6fr 0.8fr 0.7fr 0.7fr 1fr 0.6fr 1.3fr',
              gap: 12,
              padding: '14px 20px',
              alignItems: 'center',
              fontSize: 13.5,
              cursor: l.clickable ? 'pointer' : 'default',
              opacity: l.clickable ? 1 : 0.7,
              borderBottom: `1px solid ${colors.borderSoftest}`,
            }}
            onMouseEnter={(e) => {
              if (l.clickable) e.currentTarget.style.background = 'oklch(0.975 0.005 220)'
            }}
            onMouseLeave={(e) => (e.currentTarget.style.background = 'transparent')}
          >
            <div style={{ fontWeight: 600 }}>{l.gene}</div>
            <div style={{ fontFamily: mono, fontSize: 12, color: colors.textSecondary }}>{l.region}</div>
            <div style={{ fontFamily: mono, fontWeight: 600 }}>{l.motif}</div>
            <div style={{ fontFamily: mono }}>{l.a1 ?? '—'}</div>
            <div style={{ fontFamily: mono, fontWeight: 700 }}>{l.a2 ?? '—'}</div>
            <div style={{ color: colors.textSecondary }}>{l.zygosity ?? '—'}</div>
            <div style={{ fontFamily: mono, color: colors.textSecondary }}>{l.depth != null ? `${l.depth}x` : '—'}</div>
            <div>
              <StatusBadge status={l.status} />
            </div>
          </div>
        ))}
      </div>
    </div>
  )
}
