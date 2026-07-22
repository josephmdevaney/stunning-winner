import { useMemo, useState } from 'react'
import { colors } from '../theme'
import { NON_REPEAT_ID_PREFIX, type LocusSummary } from '../api'

const MAX_RESULTS = 8

function matches(l: LocusSummary, query: string): boolean {
  const q = query.toLowerCase()
  return l.gene.toLowerCase().includes(q) || l.disease.toLowerCase().includes(q) || l.region.toLowerCase().includes(q)
}

export function TopBar({
  title, sample, modeLabel, sex, loci, onSelectLocus, onSelectGeneFinding,
}: {
  title: string
  sample: string
  modeLabel: string
  sex?: string | null
  loci: LocusSummary[]
  onSelectLocus: (id: string) => void
  onSelectGeneFinding: (gene: string) => void
}) {
  const [query, setQuery] = useState('')
  // Dropdown visibility can't just track native input focus: selecting a result
  // via onMouseDown+preventDefault (below) deliberately keeps the input's real
  // DOM focus so the click registers before a blur would close the dropdown --
  // which means a later click into the same already-focused input fires no new
  // 'focus' event at all. `dismissed` instead resets on every keystroke (and on
  // refocus), so the dropdown reopens whenever there's a query to match against,
  // regardless of whether the browser considers focus to have changed.
  const [dismissed, setDismissed] = useState(false)

  const results = useMemo(() => {
    if (!query.trim()) return []
    return loci.filter((l) => l.clickable && matches(l, query)).slice(0, MAX_RESULTS)
  }, [loci, query])

  const select = (l: LocusSummary) => {
    if (l.id.startsWith(NON_REPEAT_ID_PREFIX)) onSelectGeneFinding(l.gene)
    else onSelectLocus(l.id)
    setQuery('')
    setDismissed(true)
  }

  const open = !dismissed && query.trim().length > 0

  return (
    <div
      style={{
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'space-between',
        padding: '16px 28px',
        borderBottom: `1px solid ${colors.border}`,
        background: 'white',
      }}
    >
      <div>
        <div style={{ fontSize: 19, fontWeight: 600 }}>{title}</div>
        <div style={{ fontSize: 12.5, color: colors.textSecondary, marginTop: 2 }}>
          {modeLabel ? `${modeLabel} · ` : ''}Sample {sample}
          {sex ? ` · Sex ${sex}` : ''} · Run 2026-06-29 · GRCh38
        </div>
      </div>
      <div style={{ display: 'flex', alignItems: 'center', gap: 10 }}>
        <div style={{ position: 'relative' }}>
          <div
            style={{
              display: 'flex',
              alignItems: 'center',
              gap: 8,
              padding: '8px 14px',
              border: `1px solid ${colors.borderInputAlt}`,
              borderRadius: 8,
              fontSize: 13,
              color: colors.textSecondary,
              background: 'white',
            }}
          >
            <svg width={14} height={14} viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={2} style={{ flexShrink: 0 }}>
              <circle cx={11} cy={11} r={7} />
              <line x1={21} y1={21} x2={16.65} y2={16.65} />
            </svg>
            <input
              value={query}
              onChange={(e) => {
                setQuery(e.target.value)
                setDismissed(false)
              }}
              onFocus={() => setDismissed(false)}
              onBlur={() => setDismissed(true)}
              onKeyDown={(e) => {
                if (e.key === 'Enter' && results.length > 0) select(results[0])
                else if (e.key === 'Escape') {
                  setQuery('')
                  e.currentTarget.blur()
                }
              }}
              placeholder="Search gene or locus"
              style={{
                border: 'none',
                outline: 'none',
                font: 'inherit',
                color: colors.textPrimary,
                width: 200,
                background: 'transparent',
              }}
            />
          </div>

          {open && (
            <div
              style={{
                position: 'absolute',
                top: 'calc(100% + 6px)',
                right: 0,
                width: 320,
                background: 'white',
                border: `1px solid ${colors.border}`,
                borderRadius: 10,
                boxShadow: '0 8px 24px oklch(0.2 0.02 260 / 0.12)',
                overflow: 'hidden',
                zIndex: 20,
              }}
            >
              {results.length === 0 ? (
                <div style={{ padding: '14px 16px', fontSize: 13, color: colors.textTertiary }}>No matches</div>
              ) : (
                results.map((l) => (
                  // onMouseDown (not onClick) fires before the input's onBlur, so the
                  // click registers before `focused` flips false and the dropdown unmounts.
                  <div
                    key={l.id}
                    onMouseDown={(e) => {
                      e.preventDefault()
                      select(l)
                    }}
                    style={{
                      padding: '10px 16px',
                      cursor: 'pointer',
                      borderBottom: `1px solid ${colors.borderSoftest}`,
                    }}
                    onMouseEnter={(e) => (e.currentTarget.style.background = 'oklch(0.975 0.005 220)')}
                    onMouseLeave={(e) => (e.currentTarget.style.background = 'transparent')}
                  >
                    <div style={{ display: 'flex', alignItems: 'baseline', justifyContent: 'space-between', gap: 8 }}>
                      <span style={{ fontSize: 13.5, fontWeight: 600 }}>{l.gene}</span>
                      <span style={{ fontSize: 11, fontFamily: "'IBM Plex Mono'", color: colors.textTertiary }}>{l.region}</span>
                    </div>
                    <div style={{ fontSize: 12, color: colors.textSecondary, marginTop: 1 }}>{l.disease}</div>
                  </div>
                ))
              )}
            </div>
          )}
        </div>

        <div
          style={{
            padding: '9px 16px',
            borderRadius: 8,
            background: colors.accentButtonBg,
            color: 'white',
            fontSize: 13,
            fontWeight: 600,
            cursor: 'pointer',
          }}
        >
          Export Report
        </div>
      </div>
    </div>
  )
}
