import { useEffect, useState } from 'react'
import { api, type CatalogResponse } from '../api'
import { colors, mono } from '../theme'
import { StatCard, StatGrid } from '../components/StatCard'

const PAGE_SIZE = 50

export function RepeatCatalogScreen({ mode, sample }: { mode: string; sample: string }) {
  const [data, setData] = useState<CatalogResponse | null>(null)
  const [chrom, setChrom] = useState('all')
  const [motif, setMotif] = useState('')
  const [minLen, setMinLen] = useState('')
  const [page, setPage] = useState(0)

  useEffect(() => {
    setData(null)
    api
      .catalog({ mode, sample, chrom, motif, minLen: minLen ? Number(minLen) : undefined, page, pageSize: PAGE_SIZE })
      .then(setData)
  }, [mode, sample, chrom, motif, minLen, page])

  const clearFilters = () => {
    setChrom('all')
    setMotif('')
    setMinLen('')
    setPage(0)
  }

  if (!data) {
    return (
      <div style={{ padding: 28 }}>
        <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, padding: '60px 24px', textAlign: 'center', color: colors.textSecondary, fontSize: 14 }}>
          Loading genome-wide repeat catalog…
        </div>
      </div>
    )
  }

  const maxPage = Math.max(0, Math.ceil(data.total_matching / PAGE_SIZE) - 1)
  const rangeStart = data.total_matching === 0 ? 0 : page * PAGE_SIZE + 1
  const rangeEnd = Math.min(data.total_matching, page * PAGE_SIZE + PAGE_SIZE)

  return (
    <div style={{ padding: 28, display: 'flex', flexDirection: 'column', gap: 20 }}>
      <StatGrid>
        <StatCard label="Total repeat regions" value={data.total_all.toLocaleString()} sub="Genome-wide, TRGT calls for this sample" />
        <StatCard label="Unique motifs" value={data.unique_motifs.toLocaleString()} sub="Distinct motif sequences" />
        <StatCard label="Chromosomes covered" value={String(data.chromosomes.length - 1)} sub="Autosomes + X/Y" />
        <StatCard label="Currently matching" value={data.total_matching.toLocaleString()} sub="With active filters" />
      </StatGrid>

      <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, padding: '16px 20px', display: 'flex', alignItems: 'end', gap: 16, flexWrap: 'wrap' }}>
        <Field label="Chromosome">
          <select
            value={chrom}
            onChange={(e) => {
              setChrom(e.target.value)
              setPage(0)
            }}
            style={{ padding: '8px 10px', border: `1px solid ${colors.borderInput}`, borderRadius: 7, fontSize: 13, fontFamily: mono, minWidth: 140 }}
          >
            {data.chromosomes.map((c) => (
              <option key={c} value={c}>
                {c}
              </option>
            ))}
          </select>
        </Field>
        <Field label="Motif contains">
          <input
            value={motif}
            onChange={(e) => {
              setMotif(e.target.value)
              setPage(0)
            }}
            placeholder="e.g. CAG"
            style={{ padding: '8px 10px', border: `1px solid ${colors.borderInput}`, borderRadius: 7, fontSize: 13, fontFamily: mono, width: 140 }}
          />
        </Field>
        <Field label="Min length (bp)">
          <input
            type="number"
            value={minLen}
            onChange={(e) => {
              setMinLen(e.target.value)
              setPage(0)
            }}
            style={{ padding: '8px 10px', border: `1px solid ${colors.borderInput}`, borderRadius: 7, fontSize: 13, fontFamily: mono, width: 110 }}
          />
        </Field>
        <div
          onClick={clearFilters}
          style={{ padding: '9px 16px', borderRadius: 7, border: `1px solid ${colors.borderInput}`, fontSize: 13, fontWeight: 600, color: colors.textQuaternary, cursor: 'pointer' }}
        >
          Clear filters
        </div>
        <div style={{ marginLeft: 'auto', fontSize: 12.5, color: colors.textSecondary }}>
          {data.total_matching.toLocaleString()} of {data.total_all.toLocaleString()} regions match
        </div>
      </div>

      <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, overflow: 'hidden' }}>
        <div
          style={{
            display: 'grid',
            gridTemplateColumns: '0.9fr 1fr 1fr 0.7fr 0.9fr 0.6fr 0.6fr 0.6fr 1.3fr',
            padding: '10px 20px',
            fontSize: 11,
            fontWeight: 600,
            color: colors.textSecondary,
            letterSpacing: '0.04em',
            textTransform: 'uppercase',
            borderBottom: `1px solid ${colors.borderSofter}`,
          }}
        >
          <div>Chrom</div>
          <div>Start</div>
          <div>Stop</div>
          <div>Ref len</div>
          <div>Motif</div>
          <div>A1</div>
          <div>A2</div>
          <div>Depth</div>
          <div>Structure</div>
        </div>
        {data.rows.map((r) => (
          <div
            key={`${r.chrom}:${r.start}-${r.stop}`}
            style={{
              display: 'grid',
              gridTemplateColumns: '0.9fr 1fr 1fr 0.7fr 0.9fr 0.6fr 0.6fr 0.6fr 1.3fr',
              padding: '10px 20px',
              alignItems: 'center',
              fontSize: 13,
              fontFamily: mono,
              borderBottom: `1px solid ${colors.borderSoftest}`,
            }}
          >
            <div style={{ fontWeight: 600 }}>{r.chrom}</div>
            <div>{r.start.toLocaleString()}</div>
            <div>{r.stop.toLocaleString()}</div>
            <div>{r.length}bp</div>
            <div>{r.motifs}</div>
            <div style={{ fontWeight: 600 }}>{r.a1 ?? '—'}</div>
            <div style={{ fontWeight: 600 }}>{r.a2 ?? '—'}</div>
            <div>{r.depth != null ? `${r.depth}x` : '—'}</div>
            <div style={{ color: colors.textSecondary, overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }}>{r.struc}</div>
          </div>
        ))}
        {data.rows.length === 0 && (
          <div style={{ padding: '40px 20px', textAlign: 'center', color: colors.textSecondary, fontSize: 13.5 }}>No repeats match these filters.</div>
        )}
      </div>

      <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
        <div style={{ fontSize: 12.5, color: colors.textSecondary }}>
          {data.total_matching === 0 ? 'No results' : `Showing ${rangeStart.toLocaleString()}–${rangeEnd.toLocaleString()} of ${data.total_matching.toLocaleString()}`}
        </div>
        <div style={{ display: 'flex', gap: 8 }}>
          <PageButton label="Prev" disabled={page === 0} onClick={() => setPage((p) => Math.max(0, p - 1))} />
          <PageButton label="Next" disabled={page >= maxPage} onClick={() => setPage((p) => Math.min(maxPage, p + 1))} />
        </div>
      </div>
    </div>
  )
}

function Field({ label, children }: { label: string; children: React.ReactNode }) {
  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: 6 }}>
      <label style={{ fontSize: 11, fontWeight: 600, color: colors.textSecondary, textTransform: 'uppercase', letterSpacing: '0.04em' }}>{label}</label>
      {children}
    </div>
  )
}

function PageButton({ label, disabled, onClick }: { label: string; disabled: boolean; onClick: () => void }) {
  return (
    <div
      onClick={disabled ? undefined : onClick}
      style={{
        padding: '8px 16px',
        borderRadius: 7,
        border: `1px solid ${colors.borderInput}`,
        fontSize: 13,
        fontWeight: 600,
        cursor: disabled ? 'default' : 'pointer',
        opacity: disabled ? 0.4 : 1,
        pointerEvents: disabled ? 'none' : 'auto',
      }}
    >
      {label}
    </div>
  )
}
