import { useEffect, useState } from 'react'
import { api, type CatalogResponse, type CatalogRow, type RexprtPathogenicity } from '../api'
import { colors, mono } from '../theme'
import { StatCard, StatGrid } from '../components/StatCard'

const PAGE_SIZE = 50

// Overlapping gene(s) take priority for display; otherwise the nearest gene
// with its real distance/side. "—" only when no gene reference data was
// available server-side at all -- never a fabricated "no gene nearby" claim.
function formatGene(r: CatalogRow): string {
  if (r.genes.length > 0) return r.genes.join(' / ')
  if (!r.nearest_gene) return '—'
  if (r.nearest_distance_bp === 0) return `adjacent to ${r.nearest_gene}`
  const kb = (r.nearest_distance_bp! / 1000).toFixed(1)
  return `${kb}kb ${r.nearest_direction} of ${r.nearest_gene}`
}

// Real evidence behind a RExPRT score -- shown as a hover tooltip rather than
// its own column, since it's ~8 real feature values per region across 225k
// evaluated regions.
function formatPathogenicityEvidence(p: RexprtPathogenicity): string {
  const flags = [
    p.promoter && 'promoter',
    p.utr_5 && "5'UTR",
    p.utr_3 && "3'UTR",
    p.tad && 'TAD boundary',
    p.estr && 'eSTR',
    p.open_reg_region && 'open regulatory region',
    (p.rad21 || p.smc3 || p.ctcf) && 'TFBS overlap',
  ].filter(Boolean)
  const tissue = p.tissue.replace(/_/g, ' ').toLowerCase()
  return [
    `RExPRT ensemble score ${p.ensemble_score.toFixed(2)} (SVM ${p.svm.toFixed(2)}, XGB ${p.xgb.toFixed(2)})`,
    `${p.region}${p.location ? ` · ${p.location}` : ''} of ${p.gene}`,
    `GERP ${p.gerp.toFixed(2)} · LOEUF ${p.loeuf.toFixed(2)} · pLI ${p.pli.toFixed(2)}`,
    flags.length > 0 ? flags.join(', ') : null,
    tissue !== 'no_expression' ? `expressed: ${tissue}` : null,
  ]
    .filter(Boolean)
    .join(' — ')
}

function PathogenicityBadge({ pathogenicity }: { pathogenicity: CatalogRow['pathogenicity'] }) {
  if (!pathogenicity) {
    return <span style={{ color: colors.textSecondary, fontSize: 12 }} title="Not evaluated by RExPRT (e.g. motif length outside its 3-8bp trained range)">—</span>
  }
  const pathogenic = pathogenicity.predicted_pathogenic
  return (
    <span
      title={formatPathogenicityEvidence(pathogenicity)}
      style={{
        display: 'inline-block',
        padding: '3px 9px',
        borderRadius: 999,
        fontSize: 11,
        fontWeight: 700,
        whiteSpace: 'nowrap',
        background: pathogenic ? colors.interpBg.pathogenic : colors.interpBg.normal,
        color: pathogenic ? colors.interpColor.pathogenic : colors.interpColor.normal,
      }}
    >
      {pathogenic ? 'Predicted pathogenic' : 'Predicted benign'}
    </span>
  )
}

export function RepeatCatalogScreen({ mode, sample }: { mode: string; sample: string }) {
  const [data, setData] = useState<CatalogResponse | null>(null)
  const [chrom, setChrom] = useState('all')
  const [motif, setMotif] = useState('')
  const [minLen, setMinLen] = useState('')
  const [rexprtPathogenicOnly, setRexprtPathogenicOnly] = useState(false)
  const [page, setPage] = useState(0)

  useEffect(() => {
    setData(null)
    api
      .catalog({
        mode, sample, chrom, motif, minLen: minLen ? Number(minLen) : undefined,
        rexprtPathogenicOnly, page, pageSize: PAGE_SIZE,
      })
      .then(setData)
  }, [mode, sample, chrom, motif, minLen, rexprtPathogenicOnly, page])

  const clearFilters = () => {
    setChrom('all')
    setMotif('')
    setMinLen('')
    setRexprtPathogenicOnly(false)
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
        <StatCard
          label="RExPRT predicted pathogenic"
          value={data.n_rexprt_predicted_pathogenic.toLocaleString()}
          sub={`Of ${data.n_rexprt_evaluated.toLocaleString()} evaluated (research signal, not diagnostic)`}
        />
      </StatGrid>

      <div style={{ background: colors.interpBg.uncertain, border: `1px solid ${colors.statusUncertain}`, borderRadius: 12, padding: '14px 20px', fontSize: 11.5, color: colors.textBody, lineHeight: 1.5 }}>
        <strong style={{ color: colors.interpColor.uncertain }}>Pathogenicity column: RExPRT (Danzi et al. 2024)</strong> — real, published genome-wide
        predictions scoring each repeat's reference-length genomic context (conservation, gene constraint, regulatory
        overlap, expression), not any individual sample's expanded allele. This is a research-grade candidate
        -prioritization signal for the genome-wide catalog, not a diagnostic call, and does not override this
        product's STRchive disease-panel loci, which retain their own curated clinical classification. Known
        limitation: intronic repeats are prone to false negatives -- e.g. the real C9orf72 GGGGCC expansion scores
        below RExPRT's own threshold ("Predicted benign") despite being a well-established pathogenic locus.
      </div>

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
        <Field label="RExPRT">
          <div
            onClick={() => {
              setRexprtPathogenicOnly((v) => !v)
              setPage(0)
            }}
            style={{
              padding: '8px 14px',
              borderRadius: 7,
              border: `1px solid ${rexprtPathogenicOnly ? colors.interpColor.pathogenic : colors.borderInput}`,
              fontSize: 13,
              fontWeight: 600,
              cursor: 'pointer',
              background: rexprtPathogenicOnly ? colors.interpBg.pathogenic : 'transparent',
              color: rexprtPathogenicOnly ? colors.interpColor.pathogenic : colors.textBody,
              whiteSpace: 'nowrap',
            }}
          >
            Predicted pathogenic only
          </div>
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
            gridTemplateColumns: '0.7fr 0.9fr 0.9fr 0.6fr 0.8fr 1.3fr 1fr 0.5fr 0.5fr 0.5fr 1fr',
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
          <div>Gene</div>
          <div>Pathogenicity</div>
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
              gridTemplateColumns: '0.7fr 0.9fr 0.9fr 0.6fr 0.8fr 1.3fr 1fr 0.5fr 0.5fr 0.5fr 1fr',
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
            <div
              title={formatGene(r)}
              style={{
                color: r.genes.length > 0 ? colors.textPrimary : colors.textSecondary,
                fontWeight: r.genes.length > 0 ? 600 : 400,
                overflow: 'hidden',
                textOverflow: 'ellipsis',
                whiteSpace: 'nowrap',
              }}
            >
              {formatGene(r)}
            </div>
            <div>
              <PathogenicityBadge pathogenicity={r.pathogenicity} />
            </div>
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
