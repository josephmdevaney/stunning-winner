import { useEffect, useRef, useState } from 'react'
import { colors, mono } from '../../theme'
import { paraphaseApi, type OverviewResponse } from '../api'

const STATE_LABEL: Record<string, string> = {
  not_run: 'Not yet run',
  blocked: 'Blocked',
  running: 'Running…',
  done: 'Done',
  failed: 'Failed',
}

export function ParaphaseOverviewScreen({
  sample,
  onSelectRegion,
}: {
  sample: string
  onSelectRegion: (region: string) => void
}) {
  const [data, setData] = useState<OverviewResponse | null>(null)
  const [error, setError] = useState<string | null>(null)
  const pollRef = useRef<ReturnType<typeof setInterval> | null>(null)

  const load = () => {
    paraphaseApi.overview(sample).then(setData).catch((e) => setError(String(e)))
  }

  useEffect(() => {
    setData(null)
    setError(null)
    load()
    return () => {
      if (pollRef.current) clearInterval(pollRef.current)
    }
  }, [sample])

  useEffect(() => {
    if (data?.state === 'running' && !pollRef.current) {
      pollRef.current = setInterval(load, 3000)
    }
    if (data?.state !== 'running' && pollRef.current) {
      clearInterval(pollRef.current)
      pollRef.current = null
    }
  }, [data?.state])

  const runAnnotation = () => {
    paraphaseApi.run(sample).then(load).catch((e) => setError(String(e)))
  }

  if (error) return <div style={{ padding: 28, color: colors.statusPathogenic }}>Failed to load: {error}</div>
  if (!data) return <div style={{ padding: 28, color: colors.textSecondary }}>Loading…</div>

  return (
    <div style={{ padding: 28, display: 'flex', flexDirection: 'column', gap: 20 }}>
      <StatusBanner data={data} onRun={runAnnotation} />

      {data.state === 'done' && (
        <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, overflow: 'hidden' }}>
          <div style={{ padding: '16px 20px', borderBottom: `1px solid ${colors.borderSoft}` }}>
            <div style={{ fontSize: 15, fontWeight: 600 }}>Medical regions</div>
            <div style={{ fontSize: 12, color: colors.textSecondary, marginTop: 2 }}>
              Per interpret_paraphase.py: zero haplotypes free of a pathogenic variant is the long-read
              analog of biallelic mutations for a recessive disease.
            </div>
          </div>

          <div
            style={{
              display: 'grid',
              gridTemplateColumns: '1fr 1fr 0.8fr 0.8fr 0.9fr 0.9fr 2fr',
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
            <div>Region</div>
            <div>Gene</div>
            <div>Depth</div>
            <div>Haplotypes</div>
            <div>With pathogenic</div>
            <div>Functional</div>
            <div>Info</div>
          </div>

          {data.regions.map((r, i) => {
            const clickable = r.total_haplotypes != null
            return (
              <div
                key={`${r.region}-${r.gene ?? i}`}
                onClick={clickable ? () => onSelectRegion(r.region) : undefined}
                style={{
                  display: 'grid',
                  gridTemplateColumns: '1fr 1fr 0.8fr 0.8fr 0.9fr 0.9fr 2fr',
                  gap: 10,
                  padding: '12px 20px',
                  alignItems: 'center',
                  fontSize: 12.5,
                  cursor: clickable ? 'pointer' : 'default',
                  borderBottom: `1px solid ${colors.borderSoftest}`,
                  background: r.zero_functional_haplotypes ? colors.interpBg.pathogenic : undefined,
                }}
                onMouseEnter={(e) => {
                  if (clickable) e.currentTarget.style.filter = 'brightness(0.98)'
                }}
                onMouseLeave={(e) => (e.currentTarget.style.filter = 'none')}
              >
                <div style={{ fontWeight: 600, fontFamily: mono }}>{r.region}</div>
                <div>{r.gene ?? '—'}</div>
                <div style={{ fontFamily: mono }}>{r.region_median_depth != null ? `${r.region_median_depth}x` : '—'}</div>
                <div style={{ fontFamily: mono }}>{r.total_haplotypes ?? '—'}</div>
                <div style={{ fontFamily: mono }}>{r.n_haplotypes_with_pathogenic ?? '—'}</div>
                <div
                  style={{
                    fontFamily: mono,
                    fontWeight: r.zero_functional_haplotypes ? 700 : 400,
                    color: r.zero_functional_haplotypes ? colors.statusPathogenic : undefined,
                  }}
                >
                  {r.n_haplotypes_without_pathogenic ?? '—'}
                </div>
                <div style={{ color: colors.textSecondary, overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }}>
                  {r.info ?? '—'}
                </div>
              </div>
            )
          })}
        </div>
      )}
    </div>
  )
}

function StatusBanner({ data, onRun }: { data: OverviewResponse; onRun: () => void }) {
  const bg =
    data.state === 'done' ? colors.borderSoftest
    : data.state === 'blocked' || data.state === 'failed' ? colors.interpBg.intermediate
    : colors.interpBg.uncertain
  const border = data.state === 'done' ? colors.border : data.state === 'blocked' || data.state === 'failed' ? colors.statusIntermediate : colors.borderInput

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: 8, padding: '12px 16px', borderRadius: 10, background: bg, border: `1px solid ${border}` }}>
      <div style={{ display: 'flex', alignItems: 'center', gap: 12 }}>
        <div style={{ fontSize: 12, fontWeight: 700 }}>Annotation status: {STATE_LABEL[data.state] ?? data.state}</div>
        {(data.state === 'not_run' || data.state === 'failed') && (
          <div
            onClick={onRun}
            style={{
              padding: '5px 12px',
              borderRadius: 20,
              fontSize: 11.5,
              fontWeight: 600,
              cursor: 'pointer',
              background: colors.accentButtonBg,
              color: 'white',
            }}
          >
            Run annotation
          </div>
        )}
      </div>
      {data.state === 'blocked' && <div style={{ fontSize: 11.5, color: colors.textSecondary, lineHeight: 1.5 }}>{data.reason}</div>}
      {data.state === 'running' && (
        <div style={{ fontSize: 11.5, color: colors.textSecondary, lineHeight: 1.5 }}>
          Running interpret_paraphase.py against Ensembl VEP + ClinVar (rate-limited to ~3 requests/sec) — this
          page refreshes automatically every few seconds.
        </div>
      )}
      {data.state === 'failed' && data.log_tail && (
        <pre style={{ fontSize: 10.5, color: colors.textSecondary, background: 'white', padding: 10, borderRadius: 6, overflowX: 'auto', maxHeight: 200, margin: 0 }}>
          {data.log_tail}
        </pre>
      )}
    </div>
  )
}
