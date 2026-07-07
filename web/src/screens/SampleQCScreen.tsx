import { useEffect, useState } from 'react'
import { api, type SampleQCResponse } from '../api'
import { colors, mono } from '../theme'
import { StatCard, StatGrid } from '../components/StatCard'
import { FlagBadge } from '../components/Badge'
import { UnavailablePanel } from '../components/Unavailable'

export function SampleQCScreen({ mode, sample }: { mode: string; sample: string }) {
  const [data, setData] = useState<SampleQCResponse | null>(null)

  useEffect(() => {
    setData(null)
    api.sampleQC(mode, sample).then(setData)
  }, [mode, sample])

  if (!data) return <div style={{ padding: 28 }}>Loading…</div>

  const flagged = data.rows.filter((r) => r.flag !== 'Pass').length

  return (
    <div style={{ padding: 28, display: 'flex', flexDirection: 'column', gap: 20 }}>
      <StatGrid>
        <StatCard
          label="Overall status"
          value={data.overall_status}
          sub={`${data.loci_passing}/${data.loci_total} loci sufficiently covered`}
          dotColor={data.overall_status === 'Pass' ? colors.statusNormal : colors.statusIntermediate}
        />
        <StatCard
          label="Mean locus depth"
          value={data.mean_locus_depth != null ? `${data.mean_locus_depth}x` : '—'}
          sub="STRchive panel, genotyped loci"
          dotColor={colors.cluster1}
        />
        <StatCard label="Loci genotyped" value={`${data.rows.filter((r) => r.flag !== 'No call').length}/${data.loci_total}`} sub="TRGT calls present" dotColor={colors.cluster1} />
        <StatCard
          label="Flagged loci"
          value={String(flagged)}
          sub={flagged ? 'Below 10x minimum or no call' : 'None'}
          dotColor={flagged ? colors.statusIntermediate : colors.statusNormal}
        />
      </StatGrid>

      <UnavailablePanel metrics={data.unavailable_metrics} />

      <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, overflow: 'hidden' }}>
        <div style={{ padding: '16px 20px', borderBottom: `1px solid ${colors.borderSoft}`, fontSize: 15, fontWeight: 600 }}>
          Per-locus coverage adequacy
        </div>
        <div
          style={{
            display: 'grid',
            gridTemplateColumns: '1fr 0.8fr 1fr 1.2fr',
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
          <div>Depth</div>
          <div>Min required</div>
          <div>Flag</div>
        </div>
        {data.rows.map((r) => (
          <div
            key={r.id}
            style={{
              display: 'grid',
              gridTemplateColumns: '1fr 0.8fr 1fr 1.2fr',
              padding: '12px 20px',
              alignItems: 'center',
              fontSize: 13.5,
              borderBottom: `1px solid ${colors.borderSoftest}`,
            }}
          >
            <div style={{ fontWeight: 600 }}>{r.gene}</div>
            <div style={{ fontFamily: mono }}>{r.depth != null ? `${r.depth}x` : '—'}</div>
            <div style={{ fontFamily: mono, color: colors.textSecondary }}>{r.min_required}x</div>
            <div>
              <FlagBadge flag={r.flag} />
            </div>
          </div>
        ))}
      </div>
    </div>
  )
}

