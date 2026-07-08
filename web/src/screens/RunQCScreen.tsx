import { useEffect, useState } from 'react'
import { api, type RunQCResponse } from '../api'
import { colors } from '../theme'
import { StatCard, StatGrid } from '../components/StatCard'
import { UnavailablePanel } from '../components/Unavailable'
import { Card, CardHeading } from '../components/Card'
import { ChromDepthChart } from '../charts/ChromDepthChart'

export function RunQCScreen({ mode, sample }: { mode: string; sample: string }) {
  const [data, setData] = useState<RunQCResponse | null>(null)

  useEffect(() => {
    setData(null)
    api.runQC(mode, sample).then(setData)
  }, [mode, sample])

  if (!data) return <div style={{ padding: 28 }}>Loading…</div>

  const pct = data.repeat_regions_total ? Math.round((data.repeat_regions_genotyped / data.repeat_regions_total) * 100) : 0

  return (
    <div style={{ padding: 28, display: 'flex', flexDirection: 'column', gap: 20 }}>
      <StatGrid>
        <StatCard
          label={data.scoped_to_panel ? 'Repeat regions in panel' : 'Repeat regions in catalog'}
          value={data.repeat_regions_total.toLocaleString()}
          sub={data.scoped_to_panel ? "This assay's own repeat loci" : 'Genome-wide TRGT catalog'}
          dotColor={colors.cluster1}
        />
        <StatCard
          label="Repeat regions genotyped"
          value={data.repeat_regions_genotyped.toLocaleString()}
          sub={`${pct}% of ${data.scoped_to_panel ? 'panel' : 'catalog'} called`}
          dotColor={pct > 90 ? colors.statusNormal : colors.statusIntermediate}
        />
        <StatCard
          label="Mean repeat-region depth"
          value={data.repeat_region_mean_depth != null ? `${data.repeat_region_mean_depth}x` : '—'}
          sub="From TRGT VCF SD field"
          dotColor={colors.cluster1}
        />
        <StatCard
          label="Chromosomes covered"
          value={String(data.mean_depth_by_chromosome.length)}
          sub={data.scoped_to_panel ? 'Present in this panel' : 'Autosomes + X/Y/M present in catalog'}
          dotColor={colors.cluster1}
        />
      </StatGrid>

      <div style={{ fontSize: 12.5, color: colors.textTertiary, fontStyle: 'italic' }}>{data.note}</div>

      <UnavailablePanel metrics={data.unavailable_metrics} />

      <Card>
        <CardHeading
          title={`Mean repeat-region depth by chromosome`}
          sub={
            data.scoped_to_panel
              ? "From this sample's TRGT VCF, this assay's own panel loci only -- not whole-chromosome coverage"
              : "From this sample's TRGT VCF, cataloged repeat regions only -- not whole-chromosome coverage"
          }
        />
        <ChromDepthChart data={data.mean_depth_by_chromosome} />
      </Card>
    </div>
  )
}
