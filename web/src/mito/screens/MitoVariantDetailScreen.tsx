import { useEffect, useState } from 'react'
import { colors, mono, mito, interpBg, interpColor } from '../theme'
import { mitoApi, type VariantDetailResponse, type Variant } from '../api'
import { HeteroplasmyBar } from '../charts/HeteroplasmyBar'
import { CircularMap } from '../charts/CircularMap'
import { Pileup } from '../charts/Pileup'
import { ConsensusTrack } from '../charts/ConsensusTrack'

const cardStyle: React.CSSProperties = {
  background: 'white',
  border: `1px solid ${colors.border}`,
  borderRadius: 12,
  padding: '20px 24px',
}

export function MitoVariantDetailScreen({ sample, variantId }: { sample: string; variantId: string }) {
  const [data, setData] = useState<VariantDetailResponse | null>(null)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    if (!variantId) return
    setData(null)
    setError(null)
    mitoApi.variant(variantId, sample).then(setData).catch((e) => setError(String(e)))
  }, [sample, variantId])

  if (!variantId)
    return <div style={{ padding: 28, color: colors.textSecondary }}>Select a variant from the Overview table to see its detail.</div>
  if (error) return <div style={{ padding: 28, color: colors.statusPathogenic }}>Failed to load: {error}</div>
  if (!data) return <div style={{ padding: 28, color: colors.textSecondary }}>Loading…</div>

  const v = data.variant

  return (
    <div style={{ padding: 28, display: 'grid', gridTemplateColumns: '1fr 320px', gap: 20, alignItems: 'start' }}>
      <div style={{ display: 'flex', flexDirection: 'column', gap: 20, minWidth: 0 }}>
        <div style={cardStyle}>
          <div style={{ display: 'flex', alignItems: 'baseline', gap: 12 }}>
            <div style={{ fontSize: 22, fontWeight: 700, fontFamily: mono }}>{v.position}</div>
            <div style={{ fontSize: 13, color: colors.textSecondary }}>{v.locus}</div>
          </div>
          {v.disease && <div style={{ fontSize: 13.5, color: colors.textQuaternary, marginTop: 4 }}>{v.disease}</div>}
          <div style={{ display: 'flex', gap: 28, marginTop: 14, flexWrap: 'wrap' }}>
            <HeaderStat label="Type" value={v.type} />
            <HeaderStat label="Heteroplasmy" value={`${v.heteroplasmy}%`} mono />
            <HeaderStat label="Inheritance" value={v.inheritance} />
            <HeaderStat label="Read depth" value={`${v.depth}x`} />
          </div>
        </div>

        <div style={cardStyle}>
          <div style={{ fontSize: 14.5, fontWeight: 600, marginBottom: 2 }}>Heteroplasmy fraction</div>
          <div style={{ fontSize: 12, color: colors.textTertiary, marginBottom: 14 }}>
            Proportion of mutant vs wild-type mtDNA molecules across mapped reads
          </div>
          <div style={{ width: '100%', overflowX: 'auto' }}>
            <HeteroplasmyBar het={v.heteroplasmy} />
          </div>
          <div style={{ display: 'flex', gap: 18, marginTop: 10, fontSize: 12, color: colors.textSecondary, flexWrap: 'wrap' }}>
            <Legend swatch={mito.mutant} label="Mutant / variant reads" />
            <Legend swatch={mito.wildType} label="Wild-type reads" />
            <Legend line label="Reporting threshold" />
          </div>
        </div>

        <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 20 }}>
          <div style={cardStyle}>
            <div style={{ fontSize: 14.5, fontWeight: 600, marginBottom: 2 }}>Circular mtDNA map</div>
            <div style={{ fontSize: 12, color: colors.textTertiary, marginBottom: 14 }}>
              rCRS {data.genome_length.toLocaleString()} bp · variant position highlighted
            </div>
            <CircularMap geneMap={data.gene_map} genomeLength={data.genome_length} variant={v} />
          </div>
          <div style={cardStyle}>
            <div style={{ fontSize: 14.5, fontWeight: 600, marginBottom: 2 }}>Read pileup at position</div>
            <div style={{ fontSize: 12, color: colors.textTertiary, marginBottom: 14 }}>IGV-style · one row per HiFi read (subsampled)</div>
            <Pileup reads={data.pileup_reads} />
          </div>
        </div>

        <div style={cardStyle}>
          <div style={{ fontSize: 14.5, fontWeight: 600, marginBottom: 2 }}>Base-level consensus</div>
          <div style={{ fontSize: 12, color: colors.textTertiary, marginBottom: 14 }}>
            Consensus sequence flanking the variant position, reference vs sample
          </div>
          <ConsensusTrack altBase={data.alt_base} />
        </div>
      </div>

      <div style={{ display: 'flex', flexDirection: 'column', gap: 16 }}>
        <div style={{ ...cardStyle, padding: 20 }}>
          <RailHeading>Interpretation</RailHeading>
          <div style={{ marginTop: 12, padding: '12px 14px', borderRadius: 8, background: interpBg(v.status) }}>
            <div style={{ fontSize: 14.5, fontWeight: 700, color: interpColor(v.status) }}>{v.status}</div>
            <div style={{ fontSize: 13, color: colors.textBody, marginTop: 4, lineHeight: 1.5 }}>{v.note}</div>
          </div>
          <div style={{ marginTop: 16, display: 'flex', flexDirection: 'column', gap: 10, fontSize: 13 }}>
            <KeyVal k="Classification" v={v.status} />
            <KeyVal k="Evidence level" v={v.evidence} />
            <KeyVal k="Inheritance" v={v.inheritance} />
            <KeyVal k="Caller" v={data.caller} />
            {v.ad_alt != null && v.ad_ref != null && (
              <KeyVal k="Read support" v={`${v.ad_alt} alt / ${v.ad_ref} ref`} />
            )}
          </div>
        </div>

        <div style={{ ...cardStyle, padding: 20 }}>
          <RailHeading>Population &amp; constraint</RailHeading>
          <PopulationConstraint v={v} />
        </div>

        <div style={{ ...cardStyle, padding: 20 }}>
          <RailHeading>Reference ranges</RailHeading>
          <div style={{ marginTop: 12, display: 'flex', flexDirection: 'column', gap: 8, fontSize: 13 }}>
            {data.reference_ranges.map((r) => (
              <div key={r.label} style={{ display: 'flex', justifyContent: 'space-between' }}>
                <span style={{ color: colors.textSecondary }}>{r.label}</span>
                <span style={{ fontFamily: mono }}>{r.range}</span>
              </div>
            ))}
          </div>
        </div>

        <div
          style={{
            padding: '14px 16px',
            borderRadius: 8,
            background: colors.accentButtonBg,
            color: 'white',
            fontSize: 13,
            fontWeight: 600,
            textAlign: 'center',
            cursor: 'pointer',
          }}
        >
          Add to clinical report
        </div>
      </div>
    </div>
  )
}

function HeaderStat({ label, value, mono: isMono }: { label: string; value: string; mono?: boolean }) {
  return (
    <div>
      <div style={{ fontSize: 11, fontWeight: 600, color: colors.textSecondary, textTransform: 'uppercase', letterSpacing: '0.04em' }}>{label}</div>
      <div style={{ fontWeight: isMono ? 700 : 600, fontSize: 15, marginTop: 2, fontFamily: isMono ? mono : undefined }}>{value}</div>
    </div>
  )
}

function Legend({ swatch, line, label }: { swatch?: string; line?: boolean; label: string }) {
  return (
    <div style={{ display: 'flex', alignItems: 'center', gap: 6 }}>
      {line ? (
        <div style={{ width: 14, height: 2, background: 'oklch(0.7 0.01 90)' }} />
      ) : (
        <div style={{ width: 10, height: 10, borderRadius: 2, background: swatch }} />
      )}
      {label}
    </div>
  )
}

function RailHeading({ children }: { children: React.ReactNode }) {
  return (
    <div style={{ fontSize: 13, fontWeight: 600, color: colors.textSecondary, textTransform: 'uppercase', letterSpacing: '0.04em' }}>{children}</div>
  )
}

function PopulationConstraint({ v }: { v: Variant }) {
  const gnomadChecked = 'gnomad' in v
  const constraintChecked = 'constraint' in v

  if (!gnomadChecked && !constraintChecked) {
    return (
      <div style={{ marginTop: 12, fontSize: 12.5, color: colors.textTertiary }}>
        Synthetic demo variant — not checked against real population data.
      </div>
    )
  }

  return (
    <div style={{ marginTop: 12, display: 'flex', flexDirection: 'column', gap: 14, fontSize: 13 }}>
      <div>
        <div style={{ fontSize: 11, fontWeight: 600, color: colors.textSecondary, textTransform: 'uppercase', letterSpacing: '0.04em', marginBottom: 6 }}>
          gnomAD v3.1 frequency
        </div>
        {v.gnomad ? (
          <div style={{ display: 'flex', flexDirection: 'column', gap: 6 }}>
            <KeyVal k="AF (homoplasmic)" v={`${(v.gnomad.af_hom * 100).toFixed(4)}%`} />
            <KeyVal k="AF (heteroplasmic)" v={`${(v.gnomad.af_het * 100).toFixed(4)}%`} />
            <KeyVal k="Allele number" v={v.gnomad.an.toLocaleString()} />
            {!v.gnomad.pass_qc && (
              <div style={{ fontSize: 12, color: colors.statusIntermediate, marginTop: 2 }}>
                Site flagged by gnomAD: {v.gnomad.filters.join(', ')}
              </div>
            )}
          </div>
        ) : (
          <div style={{ fontSize: 12.5, color: colors.textTertiary }}>
            Not observed in gnomAD's mtDNA callset (n=56,434 genomes).
          </div>
        )}
      </div>

      <div>
        <div style={{ fontSize: 11, fontWeight: 600, color: colors.textSecondary, textTransform: 'uppercase', letterSpacing: '0.04em', marginBottom: 6 }}>
          gnomAD v4.1 gene constraint
        </div>
        {v.constraint ? (
          <div style={{ display: 'flex', flexDirection: 'column', gap: 4 }}>
            {v.constraint.gene_rows.map((r) => (
              <div key={r.consequence} style={{ display: 'flex', justifyContent: 'space-between', gap: 8 }}>
                <span style={{ color: colors.textSecondary }}>{r.consequence}</span>
                <span style={{ fontFamily: mono }}>
                  {r.obs_exp.toFixed(3)} obs:exp
                </span>
              </div>
            ))}
            {v.constraint.regional && (
              <div style={{ fontSize: 12, color: colors.statusUncertain, marginTop: 4 }}>
                Falls within a flagged regional constraint hotspot ({v.constraint.regional.start}–{v.constraint.regional.end}, obs:exp{' '}
                {v.constraint.regional.obs_exp.toFixed(3)}).
              </div>
            )}
          </div>
        ) : (
          <div style={{ fontSize: 12.5, color: colors.textTertiary }}>
            Not covered — outside the 37 annotated mtDNA genes (e.g. the D-loop).
          </div>
        )}
      </div>
    </div>
  )
}

function KeyVal({ k, v }: { k: string; v: string }) {
  return (
    <div style={{ display: 'flex', justifyContent: 'space-between' }}>
      <span style={{ color: colors.textSecondary }}>{k}</span>
      <span style={{ fontWeight: 600 }}>{v}</span>
    </div>
  )
}
