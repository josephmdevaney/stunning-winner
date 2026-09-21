import { useEffect, useState } from 'react'
import { rohApi, type RohGeneHit, type RohOverviewResponse, type RohRegion, type RohRegionsResponse, type RohRelationshipFValue } from '../api'
import { colors, mono, rohBandColor } from '../theme'
import { StatCard, StatGrid } from '../components/StatCard'

const PAGE_SIZE = 50
const LENGTH_BANDS = ['≥5 Mb', '3–5 Mb', '1–3 Mb', '<1 Mb']

export function ROHScreen({ sample }: { sample: string }) {
  const [overview, setOverview] = useState<RohOverviewResponse | null>(null)
  const [notFound, setNotFound] = useState(false)
  const [regions, setRegions] = useState<RohRegionsResponse | null>(null)
  const [chrom, setChrom] = useState('all')
  const [minLenMb, setMinLenMb] = useState('')
  const [page, setPage] = useState(0)
  const [allRegions, setAllRegions] = useState<RohRegion[]>([])
  const [ideogramFilter, setIdeogramFilter] = useState('all')
  const [selChrom, setSelChrom] = useState('')

  useEffect(() => {
    setOverview(null)
    setNotFound(false)
    setAllRegions([])
    rohApi
      .overview(sample)
      .then(setOverview)
      .catch(() => setNotFound(true))
  }, [sample])

  useEffect(() => {
    if (!overview) return
    const withRegions = overview.by_chromosome.filter((c) => c.centromere_mid_bp != null && c.n_regions > 0)
    if (withRegions.length === 0) return
    setSelChrom((prev) => {
      if (prev && withRegions.some((c) => c.chrom === prev)) return prev
      const updChrom = overview.upd_flags[0]?.chrom
      if (updChrom && withRegions.some((c) => c.chrom === updChrom)) return updChrom
      return withRegions[0].chrom
    })
  }, [overview])

  useEffect(() => {
    if (!overview) return
    rohApi
      .regions({ sample, chrom: 'all', minLenMb: 0, page: 0, pageSize: Math.max(1, overview.stats.n_regions_total) })
      .then((r) => setAllRegions(r.rows))
  }, [sample, overview])

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

  return (
    <div style={{ padding: 28, display: 'flex', flexDirection: 'column', gap: 20 }}>
      <SexInferenceBanner sexInference={overview.sex_inference} />

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
        <StatCard label="AR candidate genes" value={String(s.n_ar_candidate_genes)} sub={`Within ROH ≥${s.ar_candidate_min_length_mb} Mb`} />
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

      <div style={{ display: 'grid', gridTemplateColumns: '1.55fr 1fr', gap: 20, alignItems: 'start' }}>
        <IdeogramCard overview={overview} allRegions={allRegions} filter={ideogramFilter} onFilterChange={setIdeogramFilter} />

        <div style={{ display: 'flex', flexDirection: 'column', gap: 16 }}>
          <ConsanguinityCard overview={overview} allRegions={allRegions} />

          <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, padding: '18px 20px' }}>
            <div style={{ fontSize: 13, fontWeight: 600, marginBottom: 4 }}>AR candidate genes</div>
            <div style={{ fontSize: 11.5, color: colors.textSecondary, marginBottom: 12 }}>
              Autosomal-recessive OMIM genes overlapping an ROH ≥{s.ar_candidate_min_length_mb} Mb (autozygosity-mapping candidates).
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
      </div>

      {selChrom && <PerChromosomeTrackCard overview={overview} allRegions={allRegions} selChrom={selChrom} onSelectChrom={setSelChrom} />}

      <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, padding: '18px 20px' }}>
        <div style={{ fontSize: 11, color: colors.textSecondary, lineHeight: 1.5 }}>{overview.non_autosomal_note}</div>
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

// ---- Sex inference banner (drives chrX/chrY hemizygous-vs-autozygous interpretation below) ----

function SexInferenceBanner({ sexInference }: { sexInference: RohOverviewResponse['sex_inference'] }) {
  const { sex } = sexInference
  const label = sex === 'M' ? 'Male (XY)' : sex === 'F' ? 'Female (XX)' : 'Unknown'
  const color = sex ? colors.textBody : colors.statusIntermediate
  const bg = sex ? colors.borderSoftest : colors.interpBg.intermediate
  const border = sex ? colors.border : colors.statusIntermediate

  let evidence: string
  if (sex && sexInference.chrX_mean_depth != null) {
    evidence =
      `Inferred from mosdepth coverage — chrX mean depth ${sexInference.chrX_mean_depth.toFixed(2)}x is ` +
      `${((sexInference.chrX_autosomal_ratio ?? 0) * 100).toFixed(0)}% of the autosomal mean ` +
      `(${sexInference.autosomal_mean_depth?.toFixed(2)}x); chrY mean depth ${sexInference.chrY_mean_depth?.toFixed(2)}x.`
  } else if (sex) {
    evidence = `From ${sexInference.method ?? 'known sample provenance'} (no mosdepth coverage on file).`
  } else {
    evidence = 'Could not be determined — no mosdepth coverage or known provenance on file. chrX/chrY ROH below cannot be resolved to hemizygosity vs. autozygosity.'
  }

  return (
    <div style={{ display: 'flex', alignItems: 'center', gap: 12, padding: '10px 16px', borderRadius: 10, background: bg, border: `1px solid ${border}` }}>
      <div style={{ fontSize: 12, fontWeight: 700, color, whiteSpace: 'nowrap' }}>Sex: {label}</div>
      <div style={{ fontSize: 11.5, color: colors.textSecondary, lineHeight: 1.4 }}>{evidence}</div>
    </div>
  )
}

// ---- Genome-wide autozygosity ideogram (project/Runs of Homozygosity.dc.html, Option 1A: buildIdeogram) ----

const IDEOGRAM_COL_TRACK = 300
const IDEOGRAM_LABEL_W = 34
const IDEOGRAM_COL_GAP = 64
const IDEOGRAM_ROW_H = 30
const IDEOGRAM_TOP_PAD = 10

function IdeogramCard({
  overview,
  allRegions,
  filter,
  onFilterChange,
}: {
  overview: RohOverviewResponse
  allRegions: RohRegion[]
  filter: string
  onFilterChange: (band: string) => void
}) {
  const imprintedChroms = overview.by_chromosome
    .filter((c) => c.imprinted)
    .map((c) => c.chrom.replace('chr', ''))
    .join('·')

  return (
    <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, padding: '20px 22px' }}>
      <div style={{ fontSize: 15, fontWeight: 600, color: colors.textPrimary }}>Genome-wide autozygosity map</div>
      <div style={{ fontSize: 12, color: colors.textTertiary, marginTop: 2 }}>ROH painted on {overview.genome_build} chromosomes, drawn to genomic scale</div>

      <div style={{ display: 'flex', gap: 8, margin: '14px 0 6px', flexWrap: 'wrap' }}>
        {['all', ...LENGTH_BANDS].map((band) => {
          const active = filter === band
          return (
            <div
              key={band}
              onClick={() => onFilterChange(band)}
              style={{
                display: 'flex',
                alignItems: 'center',
                gap: 7,
                padding: '5px 11px',
                borderRadius: 20,
                fontSize: 12,
                fontWeight: 600,
                cursor: 'pointer',
                background: active ? colors.sidebarBg : colors.chipBg,
                color: active ? 'white' : colors.textLabel,
              }}
            >
              <div
                style={{
                  width: 10,
                  height: 10,
                  borderRadius: 3,
                  background: band === 'all' ? `linear-gradient(90deg, ${colors.rohBand5mb}, ${colors.rohBand1to3mb})` : rohBandColor(band),
                }}
              />
              {band === 'all' ? 'All' : band}
            </div>
          )
        })}
      </div>

      <Ideogram byChromosome={overview.by_chromosome} allRegions={allRegions} filter={filter} />

      <div style={{ display: 'flex', alignItems: 'center', gap: 18, marginTop: 10, fontSize: 11, color: colors.textSecondary, flexWrap: 'wrap' }}>
        <div style={{ display: 'flex', alignItems: 'center', gap: 6 }}>
          <div style={{ width: 8, height: 8, borderRadius: '50%', background: colors.rohCentromere }} />
          Centromere
        </div>
        <div style={{ display: 'flex', alignItems: 'center', gap: 6 }}>
          <div style={{ width: 8, height: 8, borderRadius: 2, background: colors.rohImprinted }} />
          Imprinted chromosome ({imprintedChroms})
        </div>
        <div style={{ display: 'flex', alignItems: 'center', gap: 6 }}>
          <div style={{ width: 9, height: 9, border: `2px solid ${colors.rohUpdOutline}`, borderRadius: 2 }} />
          UPD-flagged ROH
        </div>
        {overview.by_chromosome.some((c) => c.hemizygous) && (
          <div style={{ display: 'flex', alignItems: 'center', gap: 6 }}>
            <svg width={10} height={10}>
              <defs>
                <pattern id="hemi-hatch-legend" width={4} height={4} patternUnits="userSpaceOnUse" patternTransform="rotate(45)">
                  <rect width={4} height={4} fill={colors.rohHemizygous} />
                  <line x1={0} y1={0} x2={0} y2={4} stroke={colors.rohHemizygousStroke} strokeWidth={1.4} />
                </pattern>
              </defs>
              <rect width={10} height={10} rx={2} fill="url(#hemi-hatch-legend)" stroke={colors.rohHemizygousStroke} />
            </svg>
            Hemizygous chrX/chrY (not autozygous)
          </div>
        )}
      </div>
    </div>
  )
}

function Ideogram({
  byChromosome,
  allRegions,
  filter,
}: {
  byChromosome: RohOverviewResponse['by_chromosome']
  allRegions: RohRegion[]
  filter: string
}) {
  const chroms = byChromosome.filter((c) => c.centromere_mid_bp != null && c.chrom_length_bp != null)
  if (chroms.length === 0) return null
  const maxLen = Math.max(...chroms.map((c) => c.chrom_length_bp as number))
  const half = Math.ceil(chroms.length / 2)
  const columns = [
    { list: chroms.slice(0, half), trackX: IDEOGRAM_LABEL_W },
    { list: chroms.slice(half), trackX: IDEOGRAM_LABEL_W + IDEOGRAM_COL_TRACK + IDEOGRAM_COL_GAP + IDEOGRAM_LABEL_W },
  ]
  const width = columns[1].trackX + IDEOGRAM_COL_TRACK
  const height = IDEOGRAM_TOP_PAD + half * IDEOGRAM_ROW_H + 6

  const regionsByChrom = new Map<string, RohRegion[]>()
  for (const r of allRegions) {
    if (!regionsByChrom.has(r.chrom)) regionsByChrom.set(r.chrom, [])
    regionsByChrom.get(r.chrom)!.push(r)
  }

  return (
    <svg width="100%" viewBox={`0 0 ${width} ${height}`} style={{ display: 'block' }}>
      <defs>
        <pattern id="hemi-hatch" width={4} height={4} patternUnits="userSpaceOnUse" patternTransform="rotate(45)">
          <rect width={4} height={4} fill={colors.rohHemizygous} />
          <line x1={0} y1={0} x2={0} y2={4} stroke={colors.rohHemizygousStroke} strokeWidth={1.4} />
        </pattern>
      </defs>
      {columns.map(({ list, trackX }, colIdx) =>
        list.map((c, idx) => {
          const cy = IDEOGRAM_TOP_PAD + idx * IDEOGRAM_ROW_H + IDEOGRAM_ROW_H / 2
          const barY = cy - 6
          const barH = 12
          const trackW = ((c.chrom_length_bp as number) / maxLen) * IDEOGRAM_COL_TRACK
          const cenX = trackX + ((c.centromere_mid_bp as number) / maxLen) * IDEOGRAM_COL_TRACK
          const labelColor = c.imprinted ? colors.rohImprinted : colors.textLabel
          const segs = regionsByChrom.get(c.chrom) ?? []
          return (
            <g key={`${colIdx}-${c.chrom}`}>
              <text x={trackX - 7} y={cy + 4} fontSize={11} fontWeight={600} fill={labelColor} textAnchor="end" fontFamily={mono}>
                {c.chrom.replace('chr', '')}
              </text>
              <rect x={trackX} y={barY} width={trackW} height={barH} rx={6} fill={colors.trackBg} stroke={colors.borderInputAlt} />
              <circle cx={cenX} cy={cy} r={3} fill={colors.rohCentromere} />
              {segs.map((r) => {
                if (filter !== 'all' && r.length_band !== filter) return null
                const x = trackX + (r.start / maxLen) * IDEOGRAM_COL_TRACK
                const w = Math.max((r.length_bp / maxLen) * IDEOGRAM_COL_TRACK, 1.6)
                return (
                  <g key={`${r.chrom}:${r.start}`}>
                    <rect x={x} y={barY} width={w} height={barH} rx={2} fill={c.hemizygous ? 'url(#hemi-hatch)' : rohBandColor(r.length_band)} />
                    {r.upd_flag && (
                      <rect x={x - 1.5} y={barY - 1.5} width={w + 3} height={barH + 3} rx={3} fill="none" stroke={colors.rohUpdOutline} strokeWidth={1.6} />
                    )}
                  </g>
                )
              })}
            </g>
          )
        }),
      )}
    </svg>
  )
}

// ---- Per-chromosome track (project/Runs of Homozygosity.dc.html, Option 1B: buildTrack) ----

const TRACK_W = 1000
const TRACK_LEFT_PAD = 46
const TRACK_RIGHT_PAD = 22
const TRACK_BAR_Y = 44
const TRACK_BAR_H = 22
const TRACK_RULER_Y = TRACK_BAR_Y + TRACK_BAR_H + 16
const TRACK_GENE_Y = TRACK_BAR_Y + TRACK_BAR_H + 50
const TRACK_H = 150

function PerChromosomeTrackCard({
  overview,
  allRegions,
  selChrom,
  onSelectChrom,
}: {
  overview: RohOverviewResponse
  allRegions: RohRegion[]
  selChrom: string
  onSelectChrom: (c: string) => void
}) {
  const chromsWithRegions = overview.by_chromosome.filter((c) => c.centromere_mid_bp != null && c.n_regions > 0)
  const updChroms = new Set(overview.upd_flags.map((f) => f.chrom))
  const selMeta = overview.by_chromosome.find((c) => c.chrom === selChrom)
  const selSegs = allRegions.filter((r) => r.chrom === selChrom)
  const totalMb = selSegs.reduce((a, r) => a + r.length_mb, 0)

  if (!selMeta) return null

  return (
    <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, padding: '18px 22px' }}>
      <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', marginBottom: 12, flexWrap: 'wrap', gap: 12 }}>
        <div>
          <div style={{ fontSize: 15, fontWeight: 600, color: colors.textPrimary }}>Per-chromosome track</div>
          <div style={{ fontSize: 12, color: colors.textTertiary, marginTop: 2 }}>
            {selChrom} · {selSegs.length} {selMeta.hemizygous ? 'hemizygous region' : 'ROH'} · {totalMb.toFixed(1)} Mb
            {selMeta.imprinted ? ' · imprinted chromosome' : ''}
            {selMeta.hemizygous ? ' · not autozygosity' : ''}
          </div>
        </div>
        <div style={{ display: 'flex', gap: 6, flexWrap: 'wrap', maxWidth: 640, justifyContent: 'flex-end' }}>
          {chromsWithRegions.map((c) => {
            const active = c.chrom === selChrom
            const flagged = updChroms.has(c.chrom)
            return (
              <div
                key={c.chrom}
                onClick={() => onSelectChrom(c.chrom)}
                style={{
                  padding: '5px 11px',
                  borderRadius: 7,
                  fontSize: 12,
                  fontWeight: 600,
                  fontFamily: mono,
                  cursor: 'pointer',
                  background: active ? colors.sidebarBg : flagged ? colors.interpBg.intermediate : c.hemizygous ? colors.rohHemizygous : colors.chipBg,
                  color: active ? 'white' : flagged ? colors.interpColor.intermediate : colors.textLabel,
                }}
              >
                {c.chrom.replace('chr', '')}
              </div>
            )
          })}
        </div>
      </div>
      {selMeta.hemizygous && (
        <div style={{ fontSize: 11.5, color: colors.textSecondary, background: colors.borderSoftest, borderRadius: 8, padding: '8px 12px', marginBottom: 10, lineHeight: 1.4 }}>
          Hemizygous chromosome (one copy) — these regions reflect bcftools calling a single copy as
          "homozygous," not real autozygosity.
        </div>
      )}
      <ChromosomeTrack chromMeta={selMeta} segments={selSegs} />
    </div>
  )
}

function ChromosomeTrack({
  chromMeta,
  segments,
}: {
  chromMeta: RohOverviewResponse['by_chromosome'][number]
  segments: RohRegion[]
}) {
  const L = chromMeta.chrom_length_bp
  const cen = chromMeta.centromere_mid_bp
  if (!L || cen == null) return null
  const trackW = TRACK_W - TRACK_LEFT_PAD - TRACK_RIGHT_PAD
  const scale = trackW / L
  const cenX = TRACK_LEFT_PAD + cen * scale

  const genesBySymbol = new Map<string, RohGeneHit>()
  segments.forEach((s) => s.genes.forEach((g) => {
    if (!genesBySymbol.has(g.symbol)) genesBySymbol.set(g.symbol, g)
  }))
  const genes = [...genesBySymbol.values()].sort((a, b) => a.start - b.start)

  const ticks: number[] = []
  for (let mb = 0; mb * 1_000_000 <= L; mb += 25) ticks.push(mb)

  return (
    <svg width="100%" viewBox={`0 0 ${TRACK_W} ${TRACK_H}`} style={{ display: 'block' }}>
      <defs>
        <pattern id="hemi-hatch-track" width={5} height={5} patternUnits="userSpaceOnUse" patternTransform="rotate(45)">
          <rect width={5} height={5} fill={colors.rohHemizygous} />
          <line x1={0} y1={0} x2={0} y2={5} stroke={colors.rohHemizygousStroke} strokeWidth={1.6} />
        </pattern>
      </defs>
      <rect x={TRACK_LEFT_PAD} y={TRACK_BAR_Y} width={trackW} height={TRACK_BAR_H} rx={11} fill={colors.trackBg} stroke={colors.borderInputAlt} />
      <rect x={TRACK_LEFT_PAD} y={TRACK_BAR_Y} width={4} height={TRACK_BAR_H} rx={2} fill={colors.borderInput} />
      <rect x={TRACK_LEFT_PAD + trackW - 4} y={TRACK_BAR_Y} width={4} height={TRACK_BAR_H} rx={2} fill={colors.borderInput} />
      <circle cx={cenX} cy={TRACK_BAR_Y + TRACK_BAR_H / 2} r={5} fill={colors.rohCentromere} />

      {segments.map((s) => {
        const x = TRACK_LEFT_PAD + s.start * scale
        const w = Math.max(s.length_bp * scale, 2.5)
        return (
          <g key={`${s.chrom}:${s.start}`}>
            <rect x={x} y={TRACK_BAR_Y} width={w} height={TRACK_BAR_H} rx={3} fill={s.hemizygous ? 'url(#hemi-hatch-track)' : rohBandColor(s.length_band)} />
            {s.upd_flag && (
              <rect x={x - 2} y={TRACK_BAR_Y - 2} width={w + 4} height={TRACK_BAR_H + 4} rx={4} fill="none" stroke={colors.rohUpdOutline} strokeWidth={2} />
            )}
            {w > 50 && (
              <text x={x + w / 2} y={TRACK_BAR_Y - 8} fontSize={10.5} fontWeight={600} fill={colors.textLabel} textAnchor="middle" fontFamily={mono}>
                {s.length_mb.toFixed(1)} Mb
              </text>
            )}
          </g>
        )
      })}

      <line x1={TRACK_LEFT_PAD} y1={TRACK_RULER_Y} x2={TRACK_LEFT_PAD + trackW} y2={TRACK_RULER_Y} stroke={colors.axisLine} />
      {ticks.map((mb) => {
        const x = TRACK_LEFT_PAD + mb * 1_000_000 * scale
        return (
          <g key={mb}>
            <line x1={x} y1={TRACK_RULER_Y} x2={x} y2={TRACK_RULER_Y + 4} stroke={colors.textTertiary} />
            <text x={x} y={TRACK_RULER_Y + 15} fontSize={9} fill={colors.textSecondary} textAnchor="middle" fontFamily={mono}>
              {mb}
            </text>
          </g>
        )
      })}
      <text x={TRACK_LEFT_PAD + trackW} y={TRACK_RULER_Y + 15} fontSize={9} fill={colors.textSecondary} textAnchor="end" fontFamily={mono}>
        Mb
      </text>

      {genes.length > 0 ? (
        <>
          <text x={TRACK_LEFT_PAD} y={TRACK_GENE_Y - 4} fontSize={10} fontWeight={600} fill={colors.textSecondary}>
            OMIM disease genes within ROH
          </text>
          {genes.map((g) => {
            const x = TRACK_LEFT_PAD + g.start * scale
            const col = g.inheritance === 'AR' ? colors.cluster2 : colors.cluster1
            return (
              <g key={g.symbol}>
                <line x1={x} y1={TRACK_GENE_Y + 2} x2={x} y2={TRACK_GENE_Y + 14} stroke={col} strokeWidth={1.6} />
                <circle cx={x} cy={TRACK_GENE_Y + 2} r={3} fill={col} />
                <text x={x + 6} y={TRACK_GENE_Y + 12} fontSize={10.5} fontWeight={600} fill={col} fontFamily={mono}>
                  {g.symbol}
                </text>
                <text x={x + 6} y={TRACK_GENE_Y + 23} fontSize={8.5} fill={colors.textTertiary}>
                  {g.inheritance}
                </text>
              </g>
            )
          })}
        </>
      ) : (
        <text x={TRACK_LEFT_PAD} y={TRACK_GENE_Y + 6} fontSize={11} fill={colors.textSecondary}>
          {chromMeta.hemizygous
            ? `Hemizygous chromosome — regions above reflect one copy, not autozygosity; no AR/AD gene search is run on ${chromMeta.chrom}`
            : chromMeta.imprinted
              ? 'Imprinted chromosome · no cataloged OMIM disease genes within ROH'
              : `No cataloged OMIM disease genes within ROH on ${chromMeta.chrom}`}
        </text>
      )}
    </svg>
  )
}

// ---- Consanguinity assessment: F-value rows + ladder (project/Runs of Homozygosity.dc.html: fRowsA, buildLadder) ----

function consanguinityVerdict(f3Pct: number, secondCousinPct: number) {
  if (f3Pct >= 10) {
    return {
      headline: 'Consanguinity report: indicated',
      text: `Autozygosity (${f3Pct.toFixed(2)}% at ≥3 Mb) meets or exceeds first-cousin-level autozygosity. Consider a consanguinity/close-relation history discussion and expanded autosomal-recessive candidate-gene review.`,
      bg: colors.interpBg.pathogenic,
      border: colors.statusPathogenic,
      color: colors.interpColor.pathogenic,
    }
  }
  if (f3Pct >= 2) {
    return {
      headline: 'Consanguinity report: possible',
      text: `Autozygosity (${f3Pct.toFixed(2)}% at ≥3 Mb) falls within the ACMG 2–3% consanguinity reporting range. Consider discussing consanguinity/close-relation history with the family.`,
      bg: colors.interpBg.intermediate,
      border: colors.statusIntermediate,
      color: colors.interpColor.intermediate,
    }
  }
  return {
    headline: 'Consanguinity report: not indicated',
    text: `Autozygosity (${f3Pct.toFixed(2)}% at ≥3 Mb) is below the ACMG 2–3% reporting threshold and below the second-cousin value (${secondCousinPct.toFixed(2)}%). A consanguinity statement is not indicated; ROH ≥3 Mb are provided to support autozygosity mapping against the proband's phenotype.`,
    bg: colors.interpBg.normal,
    border: colors.statusNormal,
    color: colors.interpColor.normal,
  }
}

function ConsanguinityCard({ overview, allRegions }: { overview: RohOverviewResponse; allRegions: RohRegion[] }) {
  const s = overview.stats
  const f3Pct = s.f_value_pct
  const autosomalGenomeBp = s.autosomal_genome_length_bp
  const f1Bp = allRegions.filter((r) => r.autosomal && r.length_mb >= 1).reduce((sum, r) => sum + r.length_bp, 0)
  const f1Pct = autosomalGenomeBp ? (f1Bp / autosomalGenomeBp) * 100 : 0
  const secondCousinPct = (overview.relationship_f_values.find((r) => r.relationship === 'Second cousins')?.f_value ?? 0) * 100
  const verdict = consanguinityVerdict(f3Pct, secondCousinPct)

  return (
    <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, padding: '20px 22px' }}>
      <div style={{ fontSize: 15, fontWeight: 600, color: colors.textPrimary }}>Consanguinity assessment</div>
      <div style={{ fontSize: 12, color: colors.textTertiary, marginTop: 2, marginBottom: 14 }}>Per ACMG 2021 technical standard</div>

      <div style={{ display: 'flex', flexDirection: 'column', gap: 10 }}>
        <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', padding: '11px 14px', borderRadius: 9, background: colors.borderSoftest }}>
          <div>
            <div style={{ fontSize: 12.5, fontWeight: 600, color: colors.textBody }}>All ROH ≥1 Mb</div>
            <div style={{ fontSize: 11, color: colors.textTertiary, marginTop: 1 }}>includes short LD / IBS segments</div>
          </div>
          <div style={{ fontFamily: mono, fontSize: 18, fontWeight: 600, color: colors.textQuaternary }}>{f1Pct.toFixed(2)}%</div>
        </div>
        <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', padding: '11px 14px', borderRadius: 9, background: colors.interpBg.uncertain }}>
          <div>
            <div style={{ fontSize: 12.5, fontWeight: 600, color: colors.textBody }}>IBD-informative ≥3 Mb</div>
            <div style={{ fontSize: 11, color: colors.textTertiary, marginTop: 1 }}>ACMG size threshold applied</div>
          </div>
          <div style={{ fontFamily: mono, fontSize: 18, fontWeight: 600, color: colors.accentButtonBg }}>{f3Pct.toFixed(2)}%</div>
        </div>
      </div>

      <div style={{ marginTop: 16 }}>
        <ConsanguinityLadder relationshipFValues={overview.relationship_f_values} f1Pct={f1Pct} f3Pct={f3Pct} />
      </div>

      <div style={{ marginTop: 14, padding: '12px 14px', borderRadius: 9, background: verdict.bg, border: `1px solid ${verdict.border}` }}>
        <div style={{ fontSize: 12, fontWeight: 600, color: verdict.color }}>{verdict.headline}</div>
        <div style={{ fontSize: 11.5, color: verdict.color, marginTop: 3, lineHeight: 1.5 }}>{verdict.text}</div>
      </div>
    </div>
  )
}

function ConsanguinityLadder({
  relationshipFValues,
  f1Pct,
  f3Pct,
}: {
  relationshipFValues: RohRelationshipFValue[]
  f1Pct: number
  f3Pct: number
}) {
  const W = 560
  const H = 94
  const leftPad = 10
  const rightPad = 10
  const trackW = W - leftPad - rightPad
  const max = Math.max(0.26, (f3Pct / 100) * 1.05)
  const sc = trackW / max
  const baseY = 58
  const px = (v: number) => leftPad + v * sc

  const findRel = (name: string) => relationshipFValues.find((r) => r.relationship === name)?.f_value ?? 0
  const marks: [number, string][] = [
    [findRel('Second cousins'), '2nd cousin'],
    [findRel('First cousins once removed'), '1c 1R'],
    [findRel('First cousins'), '1st cousin'],
    [findRel('Uncle-niece / half-siblings / double first cousins'), '2° / dbl-1c'],
    [findRel('Parent-child / full siblings'), '1st degree'],
  ]

  const f1 = Math.min(f1Pct / 100, max)
  const f3 = Math.min(f3Pct / 100, max)

  return (
    <svg width="100%" viewBox={`0 0 ${W} ${H}`} style={{ display: 'block' }}>
      <rect x={px(0.02)} y={baseY - 9} width={(0.03 - 0.02) * sc} height={18} fill={colors.interpBg.intermediate} opacity={0.6} />
      <rect x={px(0.1)} y={baseY - 9} width={(max - 0.1) * sc} height={18} fill={colors.interpBg.pathogenic} opacity={0.5} />
      <line x1={leftPad} y1={baseY} x2={leftPad + trackW} y2={baseY} stroke={colors.axisLine} strokeWidth={2} />
      {marks.map(([v, label], i) => {
        const x = px(v)
        const row = i % 2
        const ly = row ? baseY + 30 : baseY + 20
        return (
          <g key={label}>
            <line x1={x} y1={baseY - 9} x2={x} y2={baseY + 9} stroke={colors.textTertiary} strokeWidth={1} />
            <text x={x} y={ly} fontSize={9} fill={colors.textSecondary} textAnchor="middle" fontFamily={mono}>
              {(v * 100).toFixed(v < 0.05 ? 2 : 1)}%
            </text>
            <text x={x} y={ly + 11} fontSize={8.5} fill={colors.textTertiary} textAnchor="middle">
              {label}
            </text>
          </g>
        )
      })}
      <polygon points={`${px(f1) - 4},${baseY - 11} ${px(f1) + 4},${baseY - 11} ${px(f1)},${baseY - 4}`} fill={colors.textTertiary} />
      <polygon points={`${px(f3) - 6},${baseY - 14} ${px(f3) + 6},${baseY - 14} ${px(f3)},${baseY - 4}`} fill={colors.accentButtonBg} />
      <text x={px(f3)} y={baseY - 18} fontSize={10.5} fontWeight={600} fill={colors.accentButtonBg} textAnchor="middle" fontFamily={mono}>
        F {f3Pct.toFixed(2)}%
      </text>
    </svg>
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
              {r.hemizygous && (
                <span style={{ color: colors.textSecondary, fontWeight: 600, marginRight: 8, fontStyle: 'italic' }}>hemizygous, not autozygous</span>
              )}
              {r.upd_flag && <span style={{ color: colors.statusIntermediate, fontWeight: 600, marginRight: 8 }}>{r.upd_flag}</span>}
              {r.genes.map((g) => (
                <span key={g.symbol} style={{ color: g.inheritance === 'AR' ? colors.accentButtonBg : colors.textSecondary, marginRight: 8 }}>
                  {g.symbol} ({g.inheritance})
                </span>
              ))}
              {!r.hemizygous && !r.upd_flag && r.genes.length === 0 && <span style={{ color: colors.textTertiary }}>—</span>}
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
