import { useEffect, useState } from 'react'
import { colors, mono } from '../../theme'
import { StatCard } from '../../components/StatCard'
import { paraphaseApi, type ToolVersion, type ToolVersions } from '../api'
import { GAP_LABEL, gapColor, STATUS_LABEL, statusColor } from '../flags'
import { headerRowStyle } from '../styles'
import { Card, Chip } from '../ui'

const GRID = '1.5fr 2.6fr 1.1fr 1fr'

const SOURCE_LABEL = {
  live: 'live from GitHub',
  cache: 'GitHub, cached',
  snapshot: 'offline snapshot',
} as const

function fmtWhen(iso: string | null): string {
  if (!iso) return 'unknown date'
  const d = new Date(iso)
  return Number.isNaN(d.getTime()) ? iso : d.toLocaleString(undefined, { dateStyle: 'medium', timeStyle: 'short' })
}

export function ParaphaseVersionsScreen() {
  const [data, setData] = useState<ToolVersions | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [checking, setChecking] = useState(false)

  const load = (refresh: boolean) => {
    setChecking(true)
    paraphaseApi
      .toolVersions(refresh)
      .then((d) => {
        setData(d)
        setError(null)
      })
      .catch((e) => setError(String(e)))
      .finally(() => setChecking(false))
  }

  useEffect(() => load(false), [])

  if (error && !data) return <div style={{ padding: 28, color: colors.statusPathogenic }}>Failed to load: {error}</div>
  if (!data) return <div style={{ padding: 28, color: colors.textSecondary }}>Checking the newest releases…</div>

  const used = data.tools.filter((t) => t.used)
  const others = data.tools.filter((t) => !t.used)
  const errorCount = Object.keys(data.errors).length

  return (
    <div style={{ padding: 28, display: 'flex', flexDirection: 'column', gap: 20 }}>
      <div style={{ display: 'grid', gridTemplateColumns: 'repeat(4,1fr)', gap: 16 }}>
        <StatCard label="PacBio tools in your data flow" value={String(data.summary.used)} sub="Read, run or vendored by this suite" />
        <StatCard label="On the newest release" value={String(data.summary.current)} sub="Every sample's version matches" dotColor={colors.statusNormal} />
        <StatCard
          label="Behind the newest release"
          value={String(data.summary.behind)}
          sub="At least one sample's data came from an older version"
          valueStyle={data.summary.behind ? { color: colors.statusIntermediate } : undefined}
          dotColor={data.summary.behind ? colors.statusIntermediate : undefined}
        />
        <StatCard label="Version not recorded" value={String(data.summary.not_recorded)} sub="The files this suite reads don't say" />
      </div>

      <div style={{ display: 'flex', alignItems: 'center', gap: 14, fontSize: 12.5, color: colors.textSecondary }}>
        <div>
          Newest releases: <b>{SOURCE_LABEL[data.source]}</b> · {fmtWhen(data.fetched_at)}
          {data.source === 'snapshot' && ' — GitHub could not be reached, so this may be out of date.'}
          {errorCount > 0 && data.source !== 'snapshot' && ` · ${errorCount} tool${errorCount === 1 ? '' : 's'} could not be checked.`}
        </div>
        <div
          onClick={checking ? undefined : () => load(true)}
          style={{
            marginLeft: 'auto',
            padding: '6px 14px',
            borderRadius: 20,
            fontSize: 12,
            fontWeight: 600,
            cursor: checking ? 'default' : 'pointer',
            background: colors.accentButtonBg,
            color: 'white',
            opacity: checking ? 0.6 : 1,
          }}
        >
          {checking ? 'Checking…' : 'Check GitHub now'}
        </div>
      </div>

      <Card
        title="Tools in this suite's data flow"
        sub="The version each tool wrote into your files, against its newest GitHub release. Versions come from the files themselves (VCF headers, recorded caller strings), never assumed."
      >
        <div style={{ display: 'grid', gridTemplateColumns: GRID, gap: 14, ...headerRowStyle }}>
          <div>Tool</div>
          <div>In your data</div>
          <div>Newest release</div>
          <div>Status</div>
        </div>
        {used.map((t) => (
          <ToolRow key={t.key} tool={t} />
        ))}
      </Card>

      <Card
        title="Other PacBio tools"
        sub="Not read by this suite, listed so the whole toolchain's newest releases are in one place."
      >
        {others.map((t) => (
          <div
            key={t.key}
            style={{ display: 'grid', gridTemplateColumns: '1.5fr 3.7fr 1fr', gap: 14, padding: '12px 20px', fontSize: 12.5, alignItems: 'center', borderBottom: `1px solid ${colors.borderSoftest}` }}
          >
            <a href={t.url} target="_blank" rel="noreferrer" style={{ fontWeight: 600, color: colors.textPrimary }}>
              {t.name}
            </a>
            <div style={{ color: colors.textSecondary }}>{t.role}</div>
            <ReleaseCell tool={t} />
          </div>
        ))}
      </Card>
    </div>
  )
}

function ReleaseCell({ tool }: { tool: ToolVersion }) {
  if (!tool.latest) return <div style={{ color: colors.textTertiary }}>unknown</div>
  return (
    <div>
      <a href={tool.latest.url ?? tool.url} target="_blank" rel="noreferrer" style={{ fontFamily: mono, fontWeight: 600, color: colors.textPrimary }}>
        {tool.latest.tag}
      </a>
      {tool.latest.published && <div style={{ fontSize: 11, color: colors.textTertiary }}>{tool.latest.published}</div>}
    </div>
  )
}

function ToolRow({ tool }: { tool: ToolVersion }) {
  return (
    <div
      style={{
        display: 'grid',
        gridTemplateColumns: GRID,
        gap: 14,
        padding: '14px 20px',
        fontSize: 12.5,
        alignItems: 'start',
        borderBottom: `1px solid ${colors.borderSoftest}`,
      }}
    >
      <div>
        <a href={tool.url} target="_blank" rel="noreferrer" style={{ fontWeight: 600, color: colors.textPrimary }}>
          {tool.name}
        </a>
        <div style={{ fontSize: 11.5, color: colors.textSecondary, marginTop: 3, lineHeight: 1.45 }}>{tool.suite_use}</div>
      </div>

      <div style={{ display: 'flex', flexDirection: 'column', gap: 7 }}>
        {tool.observed.length === 0 && (
          <span style={{ color: colors.textTertiary }}>No file this suite reads records this tool's version.</span>
        )}
        {tool.observed.map((o) => (
          <div key={`${o.label}-${o.raw}`} style={{ display: 'flex', alignItems: 'baseline', gap: 8, flexWrap: 'wrap' }}>
            <span style={{ fontFamily: mono, fontWeight: 700 }}>{o.version ? `v${o.version}` : 'not recorded'}</span>
            <span style={{ color: colors.textSecondary }}>
              {o.label}
              {o.count > 1 && ` · ${o.count} samples`}
            </span>
            {o.gap && o.gap !== 'none' && (
              <Chip color={gapColor(o.gap)}>{GAP_LABEL[o.gap]}</Chip>
            )}
          </div>
        ))}
      </div>

      <ReleaseCell tool={tool} />

      <div>
        <Chip color={statusColor(tool.status, tool.worst_gap)}>{STATUS_LABEL[tool.status]}</Chip>
      </div>
    </div>
  )
}
