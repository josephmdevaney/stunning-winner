import { useEffect, useMemo, useState } from 'react'
import { colors, mono } from '../../theme'
import { StatCard } from '../../components/StatCard'
import { paraphaseApi, type JsonRegionRow, type JsonTable } from '../api'
import { FLAG_COLOR, FLAG_LABEL } from '../flags'
import { headerRowStyle } from '../styles'
import { Card, Chip } from '../ui'

type Filter = 'medical' | 'flagged' | 'all'
type SortKey = 'region' | 'copies' | 'depth'

const GRID = '1fr 1.5fr 0.8fr 0.8fr 1fr 2fr 1.6fr'

export function ParaphaseRegionsScreen({
  sample,
  onSelectRegion,
}: {
  sample: string
  onSelectRegion: (region: string) => void
}) {
  const [data, setData] = useState<JsonTable | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [filter, setFilter] = useState<Filter>('medical')
  const [query, setQuery] = useState('')
  const [sort, setSort] = useState<{ key: SortKey; dir: 1 | -1 }>({ key: 'region', dir: 1 })

  useEffect(() => {
    setData(null)
    setError(null)
    paraphaseApi.jsonTable(sample).then(setData).catch((e) => setError(String(e)))
  }, [sample])

  const counts = useMemo(() => {
    const rows = data?.regions ?? []
    return {
      medical: rows.filter((r) => r.medical).length,
      flagged: rows.filter((r) => r.flags.length > 0).length,
      all: rows.length,
    }
  }, [data])

  const shown = useMemo(() => {
    if (!data) return []
    const q = query.trim().toLowerCase()
    const value = (r: JsonRegionRow): number | string =>
      sort.key === 'region' ? r.region.toLowerCase() : sort.key === 'copies' ? r.total_cn ?? -1 : r.depth_median ?? -1
    return data.regions
      .filter((r) => (filter === 'medical' ? r.medical : filter === 'flagged' ? r.flags.length > 0 : true))
      .filter((r) => !q || r.region.toLowerCase().includes(q) || r.genes.some((g) => g.toLowerCase().includes(q)))
      .slice()
      .sort((a, b) => {
        const av = value(a), bv = value(b)
        return (av < bv ? -1 : av > bv ? 1 : 0) * sort.dir
      })
  }, [data, filter, query, sort])

  if (error) return <div style={{ padding: 28, color: colors.statusPathogenic }}>Failed to load: {error}</div>
  if (!data) return <div style={{ padding: 28, color: colors.textSecondary }}>Loading…</div>

  const toggleSort = (key: SortKey) =>
    setSort((s) => (s.key === key ? { key, dir: (s.dir * -1) as 1 | -1 } : { key, dir: 1 }))
  const arrow = (key: SortKey) => (sort.key === key ? (sort.dir === 1 ? ' ↑' : ' ↓') : '')

  return (
    <div style={{ padding: 28, display: 'flex', flexDirection: 'column', gap: 20 }}>
      <div style={{ display: 'grid', gridTemplateColumns: 'repeat(4,1fr)', gap: 16 }}>
        <StatCard
          label="Regions in this JSON"
          value={String(data.n_regions)}
          sub={`${counts.medical} medical · schema ${data.source_schema ?? '—'} as written`}
        />
        <StatCard
          label="Flagged regions"
          value={String(counts.flagged)}
          sub="Fusions, SVs, unresolved copy numbers, no data"
          valueStyle={counts.flagged ? { color: colors.statusUncertain } : undefined}
        />
        <StatCard
          label="Made by Paraphase"
          value={data.paraphase_version ?? 'Not recorded'}
          sub={data.paraphase_version ? 'From the per-region VCF header' : 'No per-region VCFs alongside this JSON'}
        />
        <StatCard label="Sample sex" value={data.sample_sex ?? '—'} sub="Inferred by Paraphase from X/Y coverage" />
      </div>

      <Card
        title="Regions"
        sub="Read straight from Paraphase's JSON — copy numbers, depth and haplotype counts as Paraphase reported them. No annotation run needed."
      >
        <div style={{ display: 'flex', alignItems: 'center', gap: 10, padding: '12px 20px', borderBottom: `1px solid ${colors.borderSofter}` }}>
          {(['medical', 'flagged', 'all'] as Filter[]).map((f) => (
            <div
              key={f}
              onClick={() => setFilter(f)}
              style={{
                padding: '5px 12px',
                borderRadius: 20,
                fontSize: 12,
                fontWeight: 600,
                cursor: 'pointer',
                background: filter === f ? colors.accentButtonBg : colors.chipBg,
                color: filter === f ? 'white' : colors.textLabel,
              }}
            >
              {f === 'medical' ? 'Medical' : f === 'flagged' ? 'Flagged' : 'All'} ({counts[f]})
            </div>
          ))}
          <input
            value={query}
            onChange={(e) => setQuery(e.target.value)}
            placeholder="Search region or gene"
            style={{
              marginLeft: 'auto',
              width: 240,
              padding: '6px 10px',
              fontSize: 12.5,
              border: `1px solid ${colors.borderInput}`,
              borderRadius: 8,
              outline: 'none',
            }}
          />
        </div>

        <div style={{ display: 'grid', gridTemplateColumns: GRID, gap: 10, ...headerRowStyle }}>
          <div style={{ cursor: 'pointer' }} onClick={() => toggleSort('region')}>Region{arrow('region')}</div>
          <div>Genes</div>
          <div style={{ cursor: 'pointer' }} onClick={() => toggleSort('copies')}>Copies{arrow('copies')}</div>
          <div>Haplotypes</div>
          <div style={{ cursor: 'pointer' }} onClick={() => toggleSort('depth')}>Depth{arrow('depth')}</div>
          <div>Headline</div>
          <div>Flags</div>
        </div>

        {shown.length === 0 && (
          <div style={{ padding: '40px 20px', textAlign: 'center', color: colors.textSecondary, fontSize: 13.5 }}>
            No regions match.
          </div>
        )}

        {shown.map((r) => (
          <div
            key={r.region}
            onClick={() => onSelectRegion(r.region)}
            style={{
              display: 'grid',
              gridTemplateColumns: GRID,
              gap: 10,
              padding: '12px 20px',
              alignItems: 'center',
              fontSize: 12.5,
              cursor: 'pointer',
              borderBottom: `1px solid ${colors.borderSoftest}`,
            }}
            onMouseEnter={(e) => (e.currentTarget.style.filter = 'brightness(0.98)')}
            onMouseLeave={(e) => (e.currentTarget.style.filter = 'none')}
          >
            <div style={{ fontWeight: 600, fontFamily: mono }}>{r.region}</div>
            <div style={{ color: colors.textSecondary, overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }}>
              {r.genes.join(', ') || '—'}
            </div>
            <div style={{ fontFamily: mono }}>
              {r.total_cn ?? '—'}
              {r.gene_cn != null && <span style={{ color: colors.textTertiary }}> ({r.gene_cn} gene)</span>}
            </div>
            <div style={{ fontFamily: mono }}>
              {r.n_haplotypes}
              {r.two_copy_haplotypes > 0 && (
                <span style={{ color: colors.textTertiary }} title="Haplotypes Paraphase inferred to exist in two identical copies, from depth">
                  {' '}· {r.two_copy_haplotypes} at 2 copies
                </span>
              )}
            </div>
            <div style={{ fontFamily: mono }}>
              {r.depth_median != null ? `${r.depth_median}x` : '—'}
              {r.depth_vs_genome != null && (
                <span style={{ color: colors.textTertiary }}> · {r.depth_vs_genome}× genome</span>
              )}
            </div>
            <div style={{ color: colors.textBody }}>{r.headline ?? '—'}</div>
            <div style={{ display: 'flex', gap: 6, flexWrap: 'wrap' }}>
              {r.flags.map((f) => (
                <Chip key={f.kind} color={FLAG_COLOR[f.kind]} title={f.text}>
                  {FLAG_LABEL[f.kind]}
                </Chip>
              ))}
            </div>
          </div>
        ))}
      </Card>
    </div>
  )
}
