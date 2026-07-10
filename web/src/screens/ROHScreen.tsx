import { useEffect, useState } from 'react'
import { rohApi, type RohOverviewResponse, type RohRegionsResponse } from '../api'
import { colors, mono } from '../theme'
import { StatCard, StatGrid } from '../components/StatCard'

const PAGE_SIZE = 50

export function ROHScreen({ sample }: { sample: string }) {
  const [overview, setOverview] = useState<RohOverviewResponse | null>(null)
  const [notFound, setNotFound] = useState(false)
  const [regions, setRegions] = useState<RohRegionsResponse | null>(null)
  const [chrom, setChrom] = useState('all')
  const [minLenMb, setMinLenMb] = useState('')
  const [page, setPage] = useState(0)

  useEffect(() => {
    setOverview(null)
    setNotFound(false)
    rohApi
      .overview(sample)
      .then(setOverview)
      .catch(() => setNotFound(true))
  }, [sample])

  useEffect(() => {
    if (notFound) return
    setRegions(null)
    rohApi
      .regions({ sample, chrom, minLenMb: minLenMb ? Number(minLenMb) : undefined, page, pageSize: PAGE_SIZE })
      .then(setRegions)
  }, [sample, chrom, minLenMb, page, notFound])

  if (notFound) {
    return (
      <div style={{ padding: 28 }}>
        <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, padding: '60px 24px', textAlign: 'center', color: colors.textSecondary, fontSize: 14 }}>
          No ROH (Runs of Homozygosity) run for sample "{sample}". Real bcftools roh output is only available for HG002.
        </div>
      </div>
    )
  }

  if (!overview) {
    return (
      <div style={{ padding: 28 }}>
        <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, padding: '60px 24px', textAlign: 'center', color: colors.textSecondary, fontSize: 14 }}>
          Loading ROH analysis…
        </div>
      </div>
    )
  }

  const s = overview.stats
  const maxChromTotal = Math.max(...overview.by_chromosome.map((c) => c.total_length_bp))

  return (
    <div style={{ padding: 28, display: 'flex', flexDirection: 'column', gap: 20 }}>
      <StatGrid>
        <StatCard
          label="Autosomal F (inbreeding coefficient)"
          value={`${s.f_value_pct}%`}
          sub={`ROH ≥${s.f_value_min_length_mb} Mb, ${overview.genome_build} autosomes`}
        />
        <StatCard label="Autosomal ROH total" value={`${(s.total_autosomal_roh_bp / 1_000_000).toFixed(1)} Mb`} sub={`${s.n_regions_autosomal.toLocaleString()} regions`} />
        <StatCard
          label="UPD screening flags"
          value={String(s.n_upd_flags)}
          sub="Heuristic, not diagnostic"
          valueStyle={s.n_upd_flags > 0 ? { color: colors.statusIntermediate } : undefined}
          dotColor={s.n_upd_flags > 0 ? colors.statusIntermediate : undefined}
        />
        <StatCard label="AR candidate genes" value={String(s.n_ar_candidate_genes)} sub={`Within ROH ≥${s.f_value_min_length_mb} Mb`} />
      </StatGrid>

      {overview.upd_flags.length > 0 && (
        <div style={{ background: colors.interpBg.intermediate, border: `1px solid ${colors.statusIntermediate}`, borderRadius: 12, padding: '16px 20px' }}>
          <div style={{ display: 'flex', alignItems: 'center', gap: 8, fontSize: 14, fontWeight: 600, color: colors.interpColor.intermediate }}>
            <div style={{ width: 8, height: 8, borderRadius: '50%', background: colors.statusIntermediate }} />
            Possible UPD — {overview.upd_flags.length} region{overview.upd_flags.length > 1 ? 's' : ''} flagged
          </div>
          {overview.upd_flags.map((f) => (
            <div key={`${f.chrom}:${f.start}`} style={{ fontSize: 12.5, color: colors.textBody, marginTop: 6, fontFamily: mono }}>
              {f.chrom}:{f.start.toLocaleString()}–{f.end.toLocaleString()} ({f.length_mb} Mb, {f.position}) — {f.upd_flag}
            </div>
          ))}
          <div style={{ fontSize: 11.5, color: colors.textSecondary, marginTop: 8, lineHeight: 1.5 }}>{overview.upd_methodology_note}</div>
        </div>
      )}

      <div style={{ display: 'grid', gridTemplateColumns: '1.1fr 0.9fr', gap: 16 }}>
        <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, padding: '18px 20px' }}>
          <div style={{ fontSize: 13, fontWeight: 600, marginBottom: 4 }}>Consanguinity context</div>
          <div style={{ fontSize: 11.5, color: colors.textSecondary, marginBottom: 12 }}>
            Theoretical inbreeding coefficient by pedigree relationship (F = (1/2)^n), vs. this sample's real ROH-derived estimate.
          </div>
          <RelationshipLadder values={overview.relationship_f_values} sampleF={s.f_value_pct} />
        </div>

        <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, padding: '18px 20px' }}>
          <div style={{ fontSize: 13, fontWeight: 600, marginBottom: 4 }}>AR candidate genes</div>
          <div style={{ fontSize: 11.5, color: colors.textSecondary, marginBottom: 12 }}>
            Autosomal-recessive OMIM genes overlapping an ROH ≥{s.f_value_min_length_mb} Mb (autozygosity-mapping candidates).
          </div>
          {overview.ar_candidate_genes.length === 0 ? (
            <div style={{ fontSize: 13, color: colors.textSecondary, padding: '12px 0' }}>None found for this sample.</div>
          ) : (
            <div style={{ display: 'flex', flexDirection: 'column', gap: 8 }}>
              {overview.ar_candidate_genes.map((g) => (
                <div key={g.symbol} style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', fontSize: 12.5, borderBottom: `1px solid ${colors.borderSoftest}`, paddingBottom: 8 }}>
                  <div>
                    <span style={{ fontWeight: 700, fontFamily: mono }}>{g.symbol}</span>
                    <span style={{ color: colors.textSecondary, marginLeft: 8 }}>{g.disease}</span>
                  </div>
                  <div style={{ color: colors.textTertiary, fontFamily: mono, whiteSpace: 'nowrap' }}>
                    {g.chrom} · {g.roh_length_mb} Mb
                  </div>
                </div>
              ))}
            </div>
          )}
        </div>
      </div>

      <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, padding: '18px 20px' }}>
        <div style={{ fontSize: 13, fontWeight: 600, marginBottom: 12 }}>ROH by chromosome</div>
        <div style={{ display: 'flex', flexDirection: 'column', gap: 6 }}>
          {overview.by_chromosome.map((c) => (
            <div key={c.chrom} style={{ display: 'grid', gridTemplateColumns: '60px 1fr 90px', alignItems: 'center', gap: 10 }}>
              <div style={{ fontSize: 11.5, fontFamily: mono, fontWeight: 600, display: 'flex', alignItems: 'center', gap: 5 }}>
                {c.chrom.replace('chr', '')}
                {c.imprinted && <div title="Imprinted chromosome" style={{ width: 5, height: 5, borderRadius: 2, background: colors.cluster2 }} />}
              </div>
              <div style={{ background: colors.trackBg, borderRadius: 4, height: 8, overflow: 'hidden' }}>
                <div
                  style={{
                    width: `${maxChromTotal ? (c.total_length_bp / maxChromTotal) * 100 : 0}%`,
                    height: '100%',
                    background: c.chrom === 'chrX' || c.chrom === 'chrY' || c.chrom === 'chrM' ? colors.sidebarDotInactive : colors.accent,
                  }}
                />
              </div>
              <div style={{ fontSize: 11, fontFamily: mono, color: colors.textSecondary, textAlign: 'right' }}>
                {(c.total_length_bp / 1_000_000).toFixed(1)} Mb ({c.n_regions})
              </div>
            </div>
          ))}
        </div>
        <div style={{ fontSize: 11, color: colors.textSecondary, marginTop: 12, lineHeight: 1.5 }}>{overview.non_autosomal_note}</div>
      </div>

      <RegionsTable
        regions={regions}
        chrom={chrom}
        minLenMb={minLenMb}
        onChangeChrom={(c) => {
          setChrom(c)
          setPage(0)
        }}
        onChangeMinLen={(v) => {
          setMinLenMb(v)
          setPage(0)
        }}
        page={page}
        onChangePage={setPage}
      />
    </div>
  )
}

function RelationshipLadder({ values, sampleF }: { values: { relationship: string; f_value_pct: number }[]; sampleF: number }) {
  const maxPct = values[0]?.f_value_pct ?? 25
  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: 7 }}>
      {values.map((v) => (
        <div key={v.relationship} style={{ display: 'grid', gridTemplateColumns: '1fr 60px', alignItems: 'center', gap: 10 }}>
          <div style={{ fontSize: 12, color: colors.textBody }}>{v.relationship}</div>
          <div style={{ position: 'relative', background: colors.trackBg, borderRadius: 4, height: 8, gridColumn: 'span 1' }}>
            <div style={{ width: `${(v.f_value_pct / maxPct) * 100}%`, height: '100%', borderRadius: 4, background: colors.borderInput }} />
          </div>
        </div>
      ))}
      <div style={{ display: 'grid', gridTemplateColumns: '1fr 60px', alignItems: 'center', gap: 10, marginTop: 4, paddingTop: 8, borderTop: `1px solid ${colors.borderSofter}` }}>
        <div style={{ fontSize: 12, fontWeight: 700, color: colors.accentButtonBg }}>This sample (real, from ROH)</div>
        <div style={{ position: 'relative', background: colors.trackBg, borderRadius: 4, height: 8 }}>
          <div style={{ width: `${Math.min(100, (sampleF / maxPct) * 100)}%`, height: '100%', borderRadius: 4, background: colors.accentButtonBg }} />
        </div>
      </div>
    </div>
  )
}

function RegionsTable({
  regions,
  chrom,
  minLenMb,
  onChangeChrom,
  onChangeMinLen,
  page,
  onChangePage,
}: {
  regions: RohRegionsResponse | null
  chrom: string
  minLenMb: string
  onChangeChrom: (c: string) => void
  onChangeMinLen: (v: string) => void
  page: number
  onChangePage: (p: number) => void
}) {
  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: 12 }}>
      <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, padding: '16px 20px', display: 'flex', alignItems: 'end', gap: 16, flexWrap: 'wrap' }}>
        <Field label="Chromosome">
          <select
            value={chrom}
            onChange={(e) => onChangeChrom(e.target.value)}
            style={{ padding: '8px 10px', border: `1px solid ${colors.borderInput}`, borderRadius: 7, fontSize: 13, fontFamily: mono, minWidth: 120 }}
          >
            {(regions?.chromosomes ?? ['all']).map((c) => (
              <option key={c} value={c}>
                {c}
              </option>
            ))}
          </select>
        </Field>
        <Field label="Min length (Mb)">
          <div style={{ display: 'flex', gap: 6 }}>
            {[0, 0.5, 1, 3, 5].map((v) => (
              <div
                key={v}
                onClick={() => onChangeMinLen(v ? String(v) : '')}
                style={{
                  padding: '7px 12px',
                  borderRadius: 20,
                  border: `1px solid ${colors.borderInput}`,
                  fontSize: 12,
                  fontWeight: 600,
                  cursor: 'pointer',
                  background: (minLenMb === '' && v === 0) || Number(minLenMb) === v ? colors.accentButtonBg : 'transparent',
                  color: (minLenMb === '' && v === 0) || Number(minLenMb) === v ? 'white' : colors.textQuaternary,
                }}
              >
                {v === 0 ? 'All' : `≥${v}`}
              </div>
            ))}
          </div>
        </Field>
        <div style={{ marginLeft: 'auto', fontSize: 12.5, color: colors.textSecondary }}>
          {regions ? `${regions.total_matching.toLocaleString()} of ${regions.total_all.toLocaleString()} regions match` : '—'}
        </div>
      </div>

      <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, overflow: 'hidden' }}>
        <div
          style={{
            display: 'grid',
            gridTemplateColumns: '0.6fr 1.4fr 0.7fr 1fr 0.7fr 0.6fr 1.6fr',
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
          <div>Region</div>
          <div>Length</div>
          <div>Position</div>
          <div>Markers</div>
          <div>Qual</div>
          <div>Genes / UPD</div>
        </div>
        {(regions?.rows ?? []).map((r) => (
          <div
            key={`${r.chrom}:${r.start}-${r.end}`}
            style={{
              display: 'grid',
              gridTemplateColumns: '0.6fr 1.4fr 0.7fr 1fr 0.7fr 0.6fr 1.6fr',
              padding: '10px 20px',
              alignItems: 'center',
              fontSize: 12.5,
              fontFamily: mono,
              borderBottom: `1px solid ${colors.borderSoftest}`,
            }}
          >
            <div style={{ fontWeight: 600 }}>{r.chrom}</div>
            <div>
              {r.start.toLocaleString()}–{r.end.toLocaleString()}
            </div>
            <div style={{ fontWeight: 600 }}>{r.length_mb} Mb</div>
            <div style={{ color: colors.textSecondary }}>{r.position ?? '—'}</div>
            <div>{r.n_markers.toLocaleString()}</div>
            <div>{r.quality.toFixed(1)}</div>
            <div style={{ overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }}>
              {r.upd_flag && <span style={{ color: colors.statusIntermediate, fontWeight: 600, marginRight: 8 }}>{r.upd_flag}</span>}
              {r.genes.map((g) => (
                <span key={g.symbol} style={{ color: g.inheritance === 'AR' ? colors.accentButtonBg : colors.textSecondary, marginRight: 8 }}>
                  {g.symbol} ({g.inheritance})
                </span>
              ))}
              {!r.upd_flag && r.genes.length === 0 && <span style={{ color: colors.textTertiary }}>—</span>}
            </div>
          </div>
        ))}
        {regions && regions.rows.length === 0 && (
          <div style={{ padding: '40px 20px', textAlign: 'center', color: colors.textSecondary, fontSize: 13.5 }}>No ROH regions match these filters.</div>
        )}
      </div>

      {regions && (
        <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
          <div style={{ fontSize: 12.5, color: colors.textSecondary }}>
            {regions.total_matching === 0
              ? 'No results'
              : `Showing ${(page * PAGE_SIZE + 1).toLocaleString()}–${Math.min(regions.total_matching, page * PAGE_SIZE + PAGE_SIZE).toLocaleString()} of ${regions.total_matching.toLocaleString()}`}
          </div>
          <div style={{ display: 'flex', gap: 8 }}>
            <PageButton label="Prev" disabled={page === 0} onClick={() => onChangePage(Math.max(0, page - 1))} />
            <PageButton
              label="Next"
              disabled={page >= Math.max(0, Math.ceil(regions.total_matching / PAGE_SIZE) - 1)}
              onClick={() => onChangePage(page + 1)}
            />
          </div>
        </div>
      )}
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
