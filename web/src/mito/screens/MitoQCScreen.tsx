import { useEffect, useState } from 'react'
import { colors, mono } from '../theme'
import { mitoApi, type QCResponse } from '../api'
import { ReadLengthHistogram } from '../charts/ReadLengthHistogram'
import { CoverageTrack } from '../charts/CoverageTrack'

const REGION_GRID = '1fr 0.8fr 1fr 1fr 1.2fr'

const cardStyle: React.CSSProperties = {
  background: 'white',
  border: `1px solid ${colors.border}`,
  borderRadius: 12,
  padding: '20px 24px',
}

export function MitoQCScreen({
  coverageYMax,
  onChangeCoverageYMax,
}: {
  coverageYMax: number
  onChangeCoverageYMax: (v: number) => void
}) {
  const [data, setData] = useState<QCResponse | null>(null)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    mitoApi.qc().then(setData).catch((e) => setError(String(e)))
  }, [])

  if (error) return <div style={{ padding: 28, color: colors.statusPathogenic }}>Failed to load: {error}</div>
  if (!data) return <div style={{ padding: 28, color: colors.textSecondary }}>Loading…</div>

  return (
    <div style={{ padding: 28, display: 'flex', flexDirection: 'column', gap: 20 }}>
      <div style={{ display: 'grid', gridTemplateColumns: 'repeat(4,1fr)', gap: 16 }}>
        {data.stats.map((s) => (
          <div key={s.label} style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, padding: '18px 20px' }}>
            <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
              <div style={{ fontSize: 12, fontWeight: 600, color: colors.textSecondary }}>{s.label}</div>
              <div style={{ width: 8, height: 8, borderRadius: '50%', background: s.dotColor }} />
            </div>
            <div style={{ fontSize: 26, fontWeight: 700, marginTop: 8 }}>{s.value}</div>
            <div style={{ fontSize: 12, color: colors.textTertiary, marginTop: 4 }}>{s.sub}</div>
          </div>
        ))}
      </div>

      <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 20 }}>
        <div style={cardStyle}>
          <div style={{ fontSize: 14.5, fontWeight: 600, marginBottom: 2 }}>HiFi read length distribution</div>
          <div style={{ fontSize: 12, color: colors.textTertiary, marginBottom: 14 }}>Full-length mtDNA amplicon / enrichment reads</div>
          <ReadLengthHistogram bins={data.read_length_histogram} />
        </div>
        <div style={cardStyle}>
          <div style={{ display: 'flex', alignItems: 'flex-start', justifyContent: 'space-between', gap: 12 }}>
            <div>
              <div style={{ fontSize: 14.5, fontWeight: 600, marginBottom: 2 }}>Coverage across the mtDNA genome</div>
              <div style={{ fontSize: 12, color: colors.textTertiary, marginBottom: 14 }}>Linear depth track, position 1–{data.genome_length.toLocaleString()} (rCRS)</div>
            </div>
            <label style={{ display: 'flex', alignItems: 'center', gap: 6, fontSize: 11, color: colors.textSecondary, whiteSpace: 'nowrap' }}>
              Y-axis max
              <input
                type="number"
                min={1000}
                max={10000}
                step={250}
                value={coverageYMax}
                onChange={(e) => onChangeCoverageYMax(clamp(Number(e.target.value), 1000, 10000))}
                style={{
                  width: 72,
                  padding: '4px 6px',
                  borderRadius: 6,
                  border: `1px solid ${colors.borderInputAlt}`,
                  fontSize: 12,
                  fontFamily: mono,
                  color: colors.textLabel,
                }}
              />
            </label>
          </div>
          <CoverageTrack points={data.coverage_track} genomeLength={data.genome_length} yMax={coverageYMax} />
        </div>
      </div>

      <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, overflow: 'hidden' }}>
        <div style={{ padding: '16px 20px', borderBottom: `1px solid ${colors.borderSoft}`, fontSize: 15, fontWeight: 600 }}>
          Per-region coverage adequacy
        </div>
        <div
          style={{
            display: 'grid',
            gridTemplateColumns: REGION_GRID,
            padding: '10px 20px',
            fontSize: 11,
            fontWeight: 600,
            color: colors.textSecondary,
            letterSpacing: '0.04em',
            textTransform: 'uppercase',
            borderBottom: `1px solid ${colors.borderSofter}`,
          }}
        >
          <div>Region</div>
          <div>Depth</div>
          <div>Min required</div>
          <div>Strand balance</div>
          <div>Flag</div>
        </div>
        {data.regions.map((r) => (
          <div
            key={r.region}
            style={{
              display: 'grid',
              gridTemplateColumns: REGION_GRID,
              padding: '12px 20px',
              alignItems: 'center',
              fontSize: 13.5,
              borderBottom: `1px solid ${colors.borderSoftest}`,
            }}
          >
            <div style={{ fontWeight: 600 }}>{r.region}</div>
            <div style={{ fontFamily: mono }}>{r.depth}x</div>
            <div style={{ fontFamily: mono, color: colors.textSecondary }}>{data.min_required}x</div>
            <div style={{ fontFamily: mono }}>{r.strandBalance}</div>
            <div>
              <span
                style={{
                  padding: '4px 10px',
                  borderRadius: 20,
                  fontSize: 11.5,
                  fontWeight: 600,
                  color: 'white',
                  background: r.flag === 'Pass' ? colors.statusNormal : colors.statusIntermediate,
                }}
              >
                {r.flag}
              </span>
            </div>
          </div>
        ))}
      </div>

      <div style={{ ...cardStyle, display: 'grid', gridTemplateColumns: 'repeat(4,1fr)', gap: 24 }}>
        {data.instrument.map((it) => (
          <div key={it.label}>
            <div style={{ fontSize: 12, fontWeight: 600, color: colors.textSecondary, textTransform: 'uppercase', letterSpacing: '0.04em' }}>{it.label}</div>
            <div style={{ fontWeight: 600, fontSize: 14.5, marginTop: 6 }}>{it.value}</div>
          </div>
        ))}
      </div>

      {data.unavailable_metrics.length > 0 && (
        <div style={{ background: 'oklch(0.98 0.005 90)', border: `1px dashed ${colors.borderInput}`, borderRadius: 12, padding: '16px 20px' }}>
          <div style={{ fontSize: 12, fontWeight: 600, color: colors.textSecondary, textTransform: 'uppercase', letterSpacing: '0.04em', marginBottom: 10 }}>
            Not available from this pipeline
          </div>
          <div style={{ display: 'flex', flexDirection: 'column', gap: 8 }}>
            {data.unavailable_metrics.map((m) => (
              <div key={m.label} style={{ fontSize: 13 }}>
                <span style={{ fontWeight: 600 }}>{m.label}</span>
                <span style={{ color: colors.textTertiary }}> — {m.reason}</span>
              </div>
            ))}
          </div>
        </div>
      )}
    </div>
  )
}

function clamp(v: number, lo: number, hi: number): number {
  if (Number.isNaN(v)) return lo
  return Math.min(hi, Math.max(lo, v))
}
