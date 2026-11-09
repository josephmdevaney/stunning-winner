import { colors, fmtInt, fmtPct, mono } from '../theme'
import { smallvarApi, type FilterState, type Overview } from '../api'
import { Button, LoadState, Pill, SectionTitle, Stat } from '../components/ui'
import { cardStyle, pageStyle } from '../components/styles'
import { useLoad } from '../useLoad'

type Open = (f: Partial<FilterState>) => void

export function OverviewScreen({ sample, onOpenVariants }: { sample: string; onOpenVariants: Open }) {
  const { data, error } = useLoad(() => smallvarApi.overview(sample), sample)
  if (!data) return <LoadState data={data} error={error} />

  return (
    <div style={pageStyle}>
      <CallSetCard o={data} />
      <FunnelCard o={data} onOpen={onOpenVariants} />
      <ClinvarCard o={data} onOpen={onOpenVariants} />
      <EvidenceCard o={data} />
    </div>
  )
}

function CallSetCard({ o }: { o: Overview }) {
  const q = o.qc
  const indels = q.insertions + q.deletions
  const maxIndel = Math.max(1, ...q.indel_sizes.map((b) => b.count))
  return (
    <div style={cardStyle}>
      <SectionTitle title="Call set" sub="Every call in the VCF that passes its FILTER, counted from the file." />
      <div style={{ display: 'flex', gap: 36, flexWrap: 'wrap', marginBottom: 18 }}>
        <Stat label="Calls" value={fmtInt(q.calls_passing)} sub={q.calls_filtered ? `${fmtInt(q.calls_filtered)} more fail the FILTER` : 'none fail the FILTER'} />
        <Stat label="SNVs" value={fmtInt(q.snv)} sub={q.mnv + q.complex ? `${fmtInt(q.mnv + q.complex)} MNV/complex` : undefined} />
        <Stat label="Indels" value={fmtInt(indels)} sub={`${fmtInt(q.insertions)} ins · ${fmtInt(q.deletions)} del`} />
        <Stat label="Het / hom" value={q.het_hom_ratio ?? '—'} sub={`${fmtInt(q.het)} het · ${fmtInt(q.hom)} hom`} />
        <Stat label="Ti / Tv" value={q.ti_tv ?? '—'} sub="SNV transitions to transversions" />
        <Stat
          label="Phased het"
          value={fmtPct(q.phased_het_fraction, 0)}
          sub={q.phased_het_fraction ? 'heterozygous calls with a phased genotype' : 'unphased: no compound-het calls'}
        />
        <Stat label="Sex" value={o.sex.sex ?? 'undetermined'} sub={o.sex.source ?? undefined} />
      </div>
      <div style={{ fontSize: 11, fontWeight: 600, color: colors.textSecondary, textTransform: 'uppercase', letterSpacing: '0.04em', marginBottom: 6 }}>
        Indel size (bp)
      </div>
      <div style={{ display: 'flex', alignItems: 'flex-end', gap: 10, height: 74 }}>
        {q.indel_sizes.map((b) => (
          <div key={b.bin} style={{ flex: 1, maxWidth: 70, textAlign: 'center' }} title={`${fmtInt(b.count)} indels of ${b.bin} bp`}>
            <div style={{ height: Math.max(2, (b.count / maxIndel) * 52), background: colors.cluster1, borderRadius: '3px 3px 0 0', opacity: 0.85 }} />
            <div style={{ fontSize: 10.5, color: colors.textSecondary, marginTop: 3, fontFamily: mono }}>{b.bin}</div>
          </div>
        ))}
      </div>
    </div>
  )
}

// Each funnel step is a real filter on the Variants screen, so the count can be reproduced there.
function presetFor(label: string, hasFrequency: boolean): Partial<FilterState> | null {
  const rare: Partial<FilterState> = hasFrequency ? { maxAf: '0.01' } : {}
  if (label === 'Pass the VCF filter') return {}
  if (label.startsWith('Rare')) return rare
  if (label === 'In a gene body') return { ...rare, inGene: true }
  if (label.startsWith('ClinVar')) return { ...rare, inGene: true, clinvar: 'plp' }
  return null
}

function FunnelCard({ o, onOpen }: { o: Overview; onOpen: Open }) {
  const top = o.funnel[0]?.count || 1
  const hasFrequency = o.funnel.some((f) => f.label.startsWith('Rare'))
  return (
    <div style={cardStyle}>
      <SectionTitle
        title="From genome to review list"
        sub="Each step adds one filter to the one before. Open a step to see exactly those variants."
      />
      {o.funnel.map((step, i) => {
        const preset = presetFor(step.label, hasFrequency)
        const prev = i > 0 ? o.funnel[i - 1].count : null
        return (
          <div key={step.label} style={{ display: 'grid', gridTemplateColumns: '260px 1fr 150px 90px', gap: 14, alignItems: 'center', padding: '7px 0', borderBottom: `1px solid ${colors.borderSoftest}` }}>
            <div style={{ fontSize: 13, fontWeight: 500 }}>{step.label}</div>
            <div style={{ background: colors.trackBg, borderRadius: 4, height: 14, overflow: 'hidden' }}>
              <div style={{ width: `${Math.max(0.6, (step.count / top) * 100)}%`, height: '100%', background: i === o.funnel.length - 1 ? colors.statusPathogenic : colors.cluster1, opacity: 0.85 }} />
            </div>
            <div style={{ fontFamily: mono, fontSize: 12.5, textAlign: 'right' }}>
              {fmtInt(step.count)}
              {prev ? <span style={{ color: colors.textTertiary }}> · {((step.count / Math.max(prev, 1)) * 100).toFixed(step.count / Math.max(prev, 1) < 0.1 ? 2 : 0)}%</span> : null}
            </div>
            <div style={{ textAlign: 'right' }}>{preset && <Button onClick={() => onOpen(preset)}>View</Button>}</div>
          </div>
        )
      })}
      {!hasFrequency && (
        <div style={{ fontSize: 12, color: colors.interpColor.intermediate, marginTop: 10 }}>
          No population frequency is installed, so there is no &ldquo;rare&rdquo; step; see the resources below.
        </div>
      )}
    </div>
  )
}

function ClinvarCard({ o, onOpen }: { o: Overview; onOpen: Open }) {
  const cv = o.clinvar
  return (
    <div style={cardStyle}>
      <SectionTitle
        title="ClinVar matches"
        sub="ClinVar's own classification and review status for variants matched by exact allele. This tool does not classify variants itself."
      />
      {!cv ? (
        <div style={{ fontSize: 13, color: colors.textSecondary }}>ClinVar is not installed for this sample, so no matches were looked for.</div>
      ) : (
        <>
          <div style={{ display: 'flex', gap: 36, flexWrap: 'wrap', marginBottom: 14 }}>
            <Stat label="Pathogenic / likely pathogenic" value={fmtInt(cv.plp)} sub="calls passing the FILTER" />
            <Stat label="Conflicting" value={fmtInt(cv.conflicting)} sub="conflicting classifications of pathogenicity" />
            {o.clinvar_plp_rare != null && <Stat label="P/LP and rare" value={fmtInt(o.clinvar_plp_rare)} sub="gnomAD AF ≤ 1% or absent" />}
          </div>
          <div style={{ display: 'flex', gap: 8, flexWrap: 'wrap', alignItems: 'center' }}>
            <Button onClick={() => onOpen({ clinvar: 'plp' })} primary>
              View all P/LP
            </Button>
            <Button onClick={() => onOpen({ clinvar: 'plp_conflicting' })}>P/LP and conflicting</Button>
            {o.top_plp_genes.length > 0 && <span style={{ fontSize: 12, color: colors.textSecondary, marginLeft: 8 }}>Genes with a P/LP call:</span>}
            {o.top_plp_genes.map((g) => (
              <span
                key={g.gene}
                onClick={() => onOpen({ clinvar: 'plp', genes: g.gene })}
                style={{ cursor: 'pointer', fontSize: 12, fontFamily: mono, padding: '3px 9px', borderRadius: 999, background: colors.chipBg, border: `1px solid ${colors.borderSoft}` }}
              >
                {g.gene}
                {g.count > 1 ? ` ×${g.count}` : ''}
              </span>
            ))}
          </div>
        </>
      )}
      {o.gnomad_coverage && (
        <div style={{ fontSize: 12.5, color: colors.textSecondary, marginTop: 14, lineHeight: 1.6 }}>
          gnomAD has a frequency for {fmtInt(o.gnomad_coverage.found)} of {fmtInt(o.gnomad_coverage.of)} calls ({fmtPct(o.gnomad_coverage.found / Math.max(o.gnomad_coverage.of, 1))}); {fmtInt(o.gnomad_coverage.common_over_1pct)} of those are
          above 1%.
        </div>
      )}
    </div>
  )
}

function EvidenceCard({ o }: { o: Overview }) {
  return (
    <div style={cardStyle}>
      <SectionTitle
        title="What this tool can and cannot say"
        sub="ACMG/AMP criteria and the evidence behind each. A criterion marked not evaluated has no data here; it is never treated as met or as ruled out."
      />
      <div style={{ display: 'grid', gridTemplateColumns: '150px 220px 90px 1fr', gap: 10, padding: '6px 0', fontSize: 10.5, fontWeight: 600, color: colors.textSecondary, letterSpacing: '0.04em', textTransform: 'uppercase', borderBottom: `1px solid ${colors.borderSofter}` }}>
        <div>Criteria</div>
        <div>Evidence</div>
        <div>Status</div>
        <div>Note</div>
      </div>
      {o.evidence.map((e) => (
        <div key={e.criteria} style={{ display: 'grid', gridTemplateColumns: '150px 220px 90px 1fr', gap: 10, padding: '8px 0', fontSize: 12.5, alignItems: 'center', borderBottom: `1px solid ${colors.borderSoftest}` }}>
          <div style={{ fontFamily: mono }}>{e.criteria}</div>
          <div>{e.evidence}</div>
          <div>
            {e.available ? (
              <Pill color={colors.interpColor.normal} bg={colors.interpBg.normal}>available</Pill>
            ) : (
              <Pill color={colors.textSecondary} bg={colors.chipBg}>not evaluated</Pill>
            )}
          </div>
          <div style={{ color: colors.textSecondary }}>{e.note}</div>
        </div>
      ))}

      <div style={{ fontSize: 12.5, fontWeight: 600, margin: '18px 0 6px' }}>Annotation resources used for this sample</div>
      {o.resources.map((r) => (
        <div key={r.key} style={{ display: 'grid', gridTemplateColumns: '210px 90px 1fr', gap: 10, padding: '6px 0', fontSize: 12.5, alignItems: 'start', borderBottom: `1px solid ${colors.borderSoftest}` }}>
          <div style={{ fontWeight: 600 }}>{r.label}</div>
          <div>
            {r.installed ? (
              <Pill color={colors.interpColor.normal} bg={colors.interpBg.normal}>installed</Pill>
            ) : (
              <Pill color={colors.textSecondary} bg={colors.chipBg}>missing</Pill>
            )}
          </div>
          <div style={{ color: colors.textSecondary, lineHeight: 1.5 }}>{r.detail}</div>
        </div>
      ))}
    </div>
  )
}
