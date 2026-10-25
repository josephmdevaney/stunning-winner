import { useEffect, useState } from 'react'
import { Card, CardHeading } from '../../components/Card'
import { colors, mono } from '../../theme'
import { epigeneticsApi, type ImprintingOverview, type RegionResult, type RegionTrack, type TrackStats, type Window } from '../api'
import { CATEGORY_HELP, TRACK_COLOR, fmtPct, fmtPos } from '../labels'
import { td, th } from '../styles'
import { CategoryChip, RegionStatusChip } from '../ui'

const W = 900
const H = 280
const PAD = { left: 46, right: 16, top: 14, bottom: 34 }

export function RegionScreen({
  sample,
  window,
  overview,
  regionId,
  onSelectRegion,
}: {
  sample: string
  window: Window
  overview: ImprintingOverview
  regionId: string
  onSelectRegion: (id: string) => void
}) {
  const region = overview.regions.find((r) => r.id === regionId) ?? overview.regions[0]
  const [track, setTrack] = useState<RegionTrack | null>(null)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    let current = true
    setTrack(null)
    setError(null)
    epigeneticsApi
      .region(sample, region.id, window)
      .then((t) => current && setTrack(t))
      .catch((e) => current && setError(String(e.message ?? e)))
    return () => {
      current = false
    }
  }, [sample, region.id, window])

  return (
    <div style={{ padding: 28, display: 'flex', flexDirection: 'column', gap: 20 }}>
      <Card>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', gap: 20 }}>
          <div>
            <select
              value={region.id}
              onChange={(e) => onSelectRegion(e.target.value)}
              style={{ fontSize: 16, fontWeight: 600, fontFamily: mono, padding: '4px 8px', borderRadius: 6, border: `1px solid ${colors.borderInput}`, background: 'white' }}
            >
              {overview.regions.map((r) => (
                <option key={r.id} value={r.id}>
                  {r.name} ({r.chrom})
                </option>
              ))}
            </select>
            <div style={{ fontSize: 12, color: colors.textSecondary, marginTop: 8, fontFamily: mono }}>
              {region.dmr.chrom}:{fmtPos(region.dmr.start)}-{fmtPos(region.dmr.end)} · GRCh38 · {region.n_cpgs_published} CpGs in the published DMR
            </div>
          </div>
          <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
            <CategoryChip category={region.category} />
            <RegionStatusChip region={region} />
          </div>
        </div>
        {region.pattern && region.interpretable && (
          <div style={{ marginTop: 14, fontSize: 13.5 }}>
            <strong>{region.pattern}</strong>
            {region.allelic_difference != null && ` · hap 1 and hap 2 differ by ${region.allelic_difference.toFixed(1)} points`}
          </div>
        )}
        {!region.interpretable && <div style={{ marginTop: 14, fontSize: 12.5, color: colors.textSecondary }}>{CATEGORY_HELP.X}.</div>}
        {region.coordinate_notes && (
          <div style={{ marginTop: 12, fontSize: 12.5, color: colors.statusIntermediate, lineHeight: 1.5 }}>Coordinates: {region.coordinate_notes}</div>
        )}
      </Card>

      <Card>
        <CardHeading
          title="CpG methylation"
          sub={`Each point is one CpG site (pb-CpG-tools). The shaded band is the DMR${track && track.analysis_window.which === 'core' ? '; the darker band is the core window the call is made on' : ''}.`}
        />
        {error && <div style={{ color: colors.statusPathogenic, fontSize: 13 }}>{error}</div>}
        {!error && !track && <div style={{ color: colors.textSecondary, fontSize: 13, padding: '40px 0', textAlign: 'center' }}>Loading CpGs…</div>}
        {track && <TrackPlot track={track} />}
      </Card>

      <Card padded={false}>
        <div style={{ padding: '16px 20px', borderBottom: `1px solid ${colors.borderSoft}`, fontSize: 15, fontWeight: 600 }}>What the call is made on</div>
        <table style={{ width: '100%', borderCollapse: 'collapse' }}>
          <thead>
            <tr>
              <th style={th}>Track</th>
              <th style={{ ...th, textAlign: 'right' }}>Methylation</th>
              <th style={{ ...th, textAlign: 'right' }}>CpGs</th>
              <th style={{ ...th, textAlign: 'right' }}>Mean depth</th>
            </tr>
          </thead>
          <tbody>
            <StatsRow label="Combined" color={TRACK_COLOR.combined} stats={region.combined} />
            <StatsRow label="Haplotype 1" color={TRACK_COLOR.hap1} stats={region.hap1} />
            <StatsRow label="Haplotype 2" color={TRACK_COLOR.hap2} stats={region.hap2} />
          </tbody>
        </table>
        <Facts region={region} />
      </Card>
    </div>
  )
}

function StatsRow({ label, color, stats }: { label: string; color: string; stats: TrackStats | null }) {
  return (
    <tr>
      <td style={td}>
        <span style={{ display: 'inline-block', width: 9, height: 9, borderRadius: '50%', background: color, marginRight: 8 }} />
        {label}
      </td>
      <td style={{ ...td, textAlign: 'right', fontVariantNumeric: 'tabular-nums' }}>{fmtPct(stats?.methylation)}</td>
      <td style={{ ...td, textAlign: 'right', fontVariantNumeric: 'tabular-nums' }}>{stats?.n_cpgs ?? '—'}</td>
      <td style={{ ...td, textAlign: 'right', fontVariantNumeric: 'tabular-nums' }}>{stats ? `${stats.mean_coverage.toFixed(1)}×` : '—'}</td>
    </tr>
  )
}

function Facts({ region }: { region: RegionResult }) {
  const rows: [string, string][] = [
    ['Imprint origin', region.imprint_origin || '—'],
    ['Methylated allele', region.methylated_allele ?? 'not established'],
    ['Phasing', region.phasing === 'phased' ? 'haplotype-resolved' : region.phasing === 'unphased' ? 'unphased (combined only)' : 'no data'],
    ['Disorders', region.disorders.length ? region.disorders.join(', ') : '—'],
  ]
  return (
    <div style={{ padding: '14px 20px 18px', display: 'grid', gridTemplateColumns: 'repeat(4,1fr)', gap: 16, borderTop: `1px solid ${colors.borderSoft}` }}>
      {rows.map(([k, v]) => (
        <div key={k}>
          <div style={{ fontSize: 10.5, fontWeight: 600, color: colors.textSecondary, letterSpacing: '0.04em', textTransform: 'uppercase' }}>{k}</div>
          <div style={{ fontSize: 13, marginTop: 3 }}>{v}</div>
        </div>
      ))}
      {region.notes.length > 0 && (
        <div style={{ gridColumn: '1 / -1', fontSize: 12.5, color: colors.textBody, lineHeight: 1.5 }}>{region.notes.join(' ')}</div>
      )}
    </div>
  )
}

function TrackPlot({ track }: { track: RegionTrack }) {
  const { start, end } = track.view
  const span = Math.max(1, end - start)
  const x = (pos: number) => PAD.left + ((pos - start) / span) * (W - PAD.left - PAD.right)
  const y = (pct: number) => PAD.top + (1 - pct / 100) * (H - PAD.top - PAD.bottom)

  const ticks = [0, 25, 50, 75, 100]
  const xTicks = [0, 0.25, 0.5, 0.75, 1].map((f) => Math.round(start + f * span))
  const haps = (['hap1', 'hap2'] as const).filter((k) => track.tracks[k].length > 0)
  const hasAny = track.tracks.combined.length > 0

  const band = (s: number, e: number, fill: string) => <rect x={x(s)} y={PAD.top} width={Math.max(1, x(e) - x(s))} height={H - PAD.top - PAD.bottom} fill={fill} />

  return (
    <div>
      <svg viewBox={`0 0 ${W} ${H}`} style={{ width: '100%', height: 'auto', display: 'block' }} role="img" aria-label={`CpG methylation across ${track.name}`}>
        {band(track.dmr.start, track.dmr.end, 'oklch(0.94 0.03 220)')}
        {track.analysis_window.which === 'core' && band(track.analysis_window.start, track.analysis_window.end, 'oklch(0.89 0.05 220)')}
        {ticks.map((t) => (
          <g key={t}>
            <line x1={PAD.left} x2={W - PAD.right} y1={y(t)} y2={y(t)} stroke={colors.gridLine} />
            <text x={PAD.left - 8} y={y(t) + 4} textAnchor="end" fontSize={11} fill={colors.textTertiary}>
              {t}%
            </text>
          </g>
        ))}
        {xTicks.map((p, i) => (
          <text key={p} x={x(p)} y={H - 12} textAnchor={i === 0 ? 'start' : i === xTicks.length - 1 ? 'end' : 'middle'} fontSize={10.5} fill={colors.textTertiary} fontFamily={mono}>
            {fmtPos(p)}
          </text>
        ))}
        <line x1={PAD.left} x2={W - PAD.right} y1={H - PAD.bottom} y2={H - PAD.bottom} stroke={colors.axisLine} />

        {hasAny &&
          [...track.tracks.combined].map((p) => <circle key={`c${p.pos}`} cx={x(p.pos)} cy={y(p.methylation)} r={2.2} fill={TRACK_COLOR.combined} opacity={0.35} />)}
        {haps.map((k) => (
          <g key={k}>
            <polyline
              fill="none"
              stroke={TRACK_COLOR[k]}
              strokeWidth={1.2}
              opacity={0.55}
              points={track.tracks[k].map((p) => `${x(p.pos)},${y(p.methylation)}`).join(' ')}
            />
            {track.tracks[k].map((p) => (
              <circle key={`${k}${p.pos}`} cx={x(p.pos)} cy={y(p.methylation)} r={3} fill={TRACK_COLOR[k]}>
                <title>{`${k} · ${fmtPos(p.pos)} · ${p.methylation.toFixed(1)}% · ${p.coverage}×`}</title>
              </circle>
            ))}
          </g>
        ))}
      </svg>
      <div style={{ display: 'flex', gap: 18, fontSize: 12, color: colors.textSecondary, marginTop: 8 }}>
        <Legend color={TRACK_COLOR.combined} label="Combined" />
        {track.tracks.hap1.length > 0 && <Legend color={TRACK_COLOR.hap1} label="Haplotype 1" />}
        {track.tracks.hap2.length > 0 && <Legend color={TRACK_COLOR.hap2} label="Haplotype 2" />}
        {!hasAny && <span>No CpGs in this window.</span>}
      </div>
    </div>
  )
}

function Legend({ color, label }: { color: string; label: string }) {
  return (
    <span style={{ display: 'inline-flex', alignItems: 'center', gap: 6 }}>
      <span style={{ width: 9, height: 9, borderRadius: '50%', background: color }} />
      {label}
    </span>
  )
}
