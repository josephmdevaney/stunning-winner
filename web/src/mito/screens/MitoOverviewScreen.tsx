import { useEffect, useMemo, useState } from 'react'
import { colors, mono, mito, filterColor } from '../theme'
import { mitoApi, type OverviewResponse, type Variant } from '../api'

const GRID = '0.9fr 0.8fr 0.9fr 1.3fr 0.85fr 0.55fr 0.7fr 0.8fr 1.3fr'

type Origin = 'all' | 'pathogenic' | 'novel'
type HetBand = 'all' | 'homoplasmic' | 'het' | 'below'

export function MitoOverviewScreen({ onSelectVariant }: { onSelectVariant: (id: string) => void }) {
  const [data, setData] = useState<OverviewResponse | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [origin, setOrigin] = useState<Origin>('all')
  const [consequence, setConsequence] = useState('all')
  const [hetBand, setHetBand] = useState<HetBand>('all')

  useEffect(() => {
    mitoApi.overview().then(setData).catch((e) => setError(String(e)))
  }, [])

  const filtered = useMemo(() => {
    if (!data) return []
    return data.variants
      .filter((v) => {
        if (origin === 'pathogenic' && v.origin !== 'pathogenic') return false
        if (origin === 'novel' && v.origin !== 'novel') return false
        if (consequence !== 'all' && v.consequence !== consequence) return false
        if (hetBand === 'homoplasmic' && v.heteroplasmy < 95) return false
        if (hetBand === 'het' && (v.heteroplasmy < 5 || v.heteroplasmy >= 95)) return false
        if (hetBand === 'below' && v.heteroplasmy >= 5) return false
        return true
      })
      .slice()
      .sort((a, b) => a.pos - b.pos)
  }, [data, origin, consequence, hetBand])

  if (error) return <div style={{ padding: 28, color: colors.statusPathogenic }}>Failed to load: {error}</div>
  if (!data) return <div style={{ padding: 28, color: colors.textSecondary }}>Loading…</div>

  const pathoCount = data.stats.pathogenic_findings

  return (
    <div style={{ padding: 28, display: 'flex', flexDirection: 'column', gap: 22 }}>
      <div style={{ display: 'grid', gridTemplateColumns: 'repeat(4,1fr)', gap: 16 }}>
        <MitoStat label="Variants called" value={String(data.stats.variants_called)} sub="Genome-wide, vs rCRS" />
        <MitoStat
          label="Pathogenic findings"
          value={String(pathoCount)}
          sub="Above heteroplasmy reporting threshold"
          valueColor={pathoCount ? colors.statusPathogenic : undefined}
        />
        <MitoStat label="Haplogroup" value={data.haplogroup.name} sub={`${data.haplogroup.confidence}% confidence`} valueMono />
        <MitoStat label="Sample QC" value={data.stats.qc_status} sub="Genome-wide coverage sufficient" valueColor={colors.statusNormal} />
      </div>

      <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, overflow: 'hidden' }}>
        <div
          style={{
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'space-between',
            padding: '16px 20px',
            borderBottom: `1px solid ${colors.borderSoft}`,
            flexWrap: 'wrap',
            gap: 12,
          }}
        >
          <div>
            <div style={{ fontSize: 15, fontWeight: 600 }}>mtDNA variant calls</div>
            <div style={{ fontSize: 12, color: colors.textSecondary, marginTop: 2 }}>
              Showing {filtered.length} of {data.variants.length} calls · sorted by position
            </div>
          </div>
          <div style={{ display: 'flex', gap: 16, alignItems: 'center', flexWrap: 'wrap' }}>
            <div style={{ display: 'flex', gap: 8 }}>
              {(
                [
                  { key: 'all', label: 'All' },
                  { key: 'pathogenic', label: 'Known pathogenic' },
                  { key: 'novel', label: 'Novel / benign' },
                ] as { key: Origin; label: string }[]
              ).map((c) => (
                <div
                  key={c.key}
                  onClick={() => setOrigin(c.key)}
                  style={{
                    padding: '6px 12px',
                    borderRadius: 20,
                    fontSize: 12,
                    fontWeight: 600,
                    cursor: 'pointer',
                    background: origin === c.key ? colors.accentButtonBg : mito.chipInactiveBg,
                    color: origin === c.key ? 'white' : colors.textLabel,
                  }}
                >
                  {c.label}
                </div>
              ))}
            </div>
            <select value={consequence} onChange={(e) => setConsequence(e.target.value)} style={selectStyle}>
              <option value="all">all</option>
              {data.consequence_options.map((o) => (
                <option key={o} value={o}>
                  {o}
                </option>
              ))}
            </select>
            <select value={hetBand} onChange={(e) => setHetBand(e.target.value as HetBand)} style={selectStyle}>
              <option value="all">All heteroplasmy</option>
              <option value="homoplasmic">Homoplasmic (≥95%)</option>
              <option value="het">Heteroplasmic (5–95%)</option>
              <option value="below">Below threshold (&lt;5%)</option>
            </select>
          </div>
        </div>

        <div
          style={{
            display: 'grid',
            gridTemplateColumns: GRID,
            gap: 10,
            padding: '10px 20px',
            fontSize: 10.5,
            fontWeight: 600,
            color: colors.textSecondary,
            letterSpacing: '0.04em',
            textTransform: 'uppercase',
            borderBottom: `1px solid ${colors.borderSofter}`,
          }}
        >
          <div>Position</div>
          <div>Change</div>
          <div>Locus</div>
          <div>Consequence</div>
          <div>Heteroplasmy</div>
          <div>QUAL</div>
          <div>Filter</div>
          <div>Pop freq</div>
          <div>Annotation</div>
        </div>

        {filtered.map((v) => (
          <VariantRow key={v.id} v={v} onClick={() => onSelectVariant(v.id)} />
        ))}
      </div>

      <div
        style={{
          background: 'white',
          border: `1px solid ${colors.border}`,
          borderRadius: 12,
          padding: '20px 24px',
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'space-between',
          gap: 24,
        }}
      >
        <div>
          <div style={{ fontSize: 13, fontWeight: 600, color: colors.textSecondary, textTransform: 'uppercase', letterSpacing: '0.04em' }}>
            Haplogroup assignment
          </div>
          <div style={{ display: 'flex', alignItems: 'baseline', gap: 10, marginTop: 8 }}>
            <div style={{ fontSize: 24, fontWeight: 700, fontFamily: mono }}>{data.haplogroup.name}</div>
            <div style={{ fontSize: 13, color: colors.textSecondary }}>
              {data.haplogroup.confidence}% confidence · {data.haplogroup.region}
            </div>
          </div>
        </div>
        <div style={{ display: 'flex', gap: 8, flexWrap: 'wrap', maxWidth: 420, justifyContent: 'flex-end' }}>
          {data.haplogroup.markers.map((m) => (
            <div
              key={m}
              style={{ padding: '5px 10px', borderRadius: 6, background: mito.markerBg, fontFamily: mono, fontSize: 12, color: mito.markerText, fontWeight: 600 }}
            >
              {m}
            </div>
          ))}
        </div>
      </div>
    </div>
  )
}

const selectStyle: React.CSSProperties = {
  padding: '6px 10px',
  borderRadius: 8,
  border: `1px solid ${colors.borderInputAlt}`,
  fontSize: 12,
  color: colors.textLabel,
  background: 'white',
}

function VariantRow({ v, onClick }: { v: Variant; onClick: () => void }) {
  const hetLabel = v.isDeletion ? `${v.heteroplasmy}% del` : `${v.heteroplasmy}%`
  return (
    <div
      onClick={onClick}
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
      onMouseEnter={(e) => (e.currentTarget.style.background = 'oklch(0.975 0.005 220)')}
      onMouseLeave={(e) => (e.currentTarget.style.background = 'transparent')}
    >
      <div style={{ fontFamily: mono, fontWeight: 600 }}>{v.position}</div>
      <div style={{ fontFamily: mono, color: colors.textSecondary }}>{v.change}</div>
      <div>
        <div style={{ fontWeight: 600 }}>{v.locus}</div>
        {v.disease && <div style={{ fontSize: 10.5, color: colors.textSecondary, marginTop: 1 }}>{v.disease}</div>}
      </div>
      <div style={{ color: colors.textLabel }}>{v.consequence}</div>
      <div style={{ fontFamily: mono, fontWeight: 700 }}>{hetLabel}</div>
      <div style={{ fontFamily: mono, color: colors.textSecondary }}>{v.qual}</div>
      <div style={{ fontFamily: mono, color: filterColor(v.filter), fontWeight: 600 }}>{v.filter}</div>
      <div style={{ fontFamily: mono, color: colors.textSecondary }}>{v.gnomadFreq}</div>
      <div>
        <span style={{ padding: '4px 10px', borderRadius: 20, fontSize: 11, fontWeight: 600, color: 'white', background: v.statusColor }}>
          {v.status}
        </span>
      </div>
    </div>
  )
}

function MitoStat({
  label,
  value,
  sub,
  valueColor,
  valueMono,
}: {
  label: string
  value: string
  sub: string
  valueColor?: string
  valueMono?: boolean
}) {
  return (
    <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, padding: '18px 20px' }}>
      <div style={{ fontSize: 12, fontWeight: 600, color: colors.textSecondary, letterSpacing: '0.02em' }}>{label}</div>
      <div style={{ fontSize: 28, fontWeight: 700, marginTop: 8, color: valueColor, fontFamily: valueMono ? mono : undefined }}>{value}</div>
      <div style={{ fontSize: 12, color: colors.textTertiary, marginTop: 4 }}>{sub}</div>
    </div>
  )
}
