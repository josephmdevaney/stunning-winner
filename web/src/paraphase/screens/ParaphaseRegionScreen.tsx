import { useEffect, useState } from 'react'
import { colors, mono } from '../../theme'
import { paraphaseApi, type RegionResponse, type VariantRecord } from '../api'

const BOUNDARY_LABEL: Record<string, string> = {
  fully_contained: 'Fully contained',
  partial: 'Partial',
  truncated: 'Truncated (pseudogene-like)',
  no_overlap: 'No overlap',
}

const BOUNDARY_COLOR: Record<string, string> = {
  fully_contained: colors.statusNormal,
  partial: colors.statusIntermediate,
  truncated: colors.statusPathogenic,
  no_overlap: colors.textTertiary,
}

export function ParaphaseRegionScreen({ sample, region }: { sample: string; region: string }) {
  const [data, setData] = useState<RegionResponse | null>(null)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    setData(null)
    setError(null)
    paraphaseApi.region(region, sample).then(setData).catch((e) => setError(String(e)))
  }, [sample, region])

  if (error) return <div style={{ padding: 28, color: colors.statusPathogenic }}>Failed to load: {error}</div>
  if (!data) return <div style={{ padding: 28, color: colors.textSecondary }}>Loading…</div>

  return (
    <div style={{ padding: 28, display: 'flex', flexDirection: 'column', gap: 16 }}>
      <div style={{ fontSize: 15, fontWeight: 600 }}>{region}</div>
      <div style={{ fontSize: 12, color: colors.textSecondary, marginTop: -10 }}>
        Boundary check explains whether each haplotype covers the whole gene of interest (fully_contained),
        only part of it (partial), is truncated mid-transcript (truncated — most likely a nonfunctional
        pseudogene copy), or doesn't overlap it at all (no_overlap).
      </div>

      {data.haplotypes.length === 0 ? (
        <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, padding: '40px 20px', textAlign: 'center', color: colors.textSecondary, fontSize: 13.5 }}>
          No haplotype rows for this region.
        </div>
      ) : (
        <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, overflow: 'hidden' }}>
          <div
            style={{
              display: 'grid',
              gridTemplateColumns: '1fr 1.4fr 1.4fr 1.2fr 3fr',
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
            <div>Gene</div>
            <div>Haplotype</div>
            <div>Haplotype boundary</div>
            <div>Boundary check</div>
            <div>Pathogenic variants</div>
          </div>

          {data.haplotypes.map((h, i) => (
            <div
              key={`${h.haplotype}-${i}`}
              style={{
                display: 'grid',
                gridTemplateColumns: '1fr 1.4fr 1.4fr 1.2fr 3fr',
                gap: 10,
                padding: '12px 20px',
                alignItems: 'start',
                fontSize: 12.5,
                borderBottom: `1px solid ${colors.borderSoftest}`,
              }}
            >
              <div style={{ fontWeight: 600 }}>{h.gene}</div>
              <div style={{ fontFamily: mono }}>{h.haplotype}</div>
              <div style={{ fontFamily: mono, color: colors.textSecondary }}>{h.haplotype_boundary}</div>
              <div style={{ fontWeight: 600, color: BOUNDARY_COLOR[h.boundary_check] ?? colors.textPrimary }}>
                {BOUNDARY_LABEL[h.boundary_check] ?? h.boundary_check}
              </div>
              <div style={{ display: 'flex', flexDirection: 'column', gap: 6 }}>
                {h.known_nonfunctional && (
                  <span style={{ color: colors.statusIntermediate, fontWeight: 600 }}>
                    Known nonfunctional haplotype (from Paraphase JSON)
                  </span>
                )}
                {h.pathogenic_variants.length === 0 && !h.known_nonfunctional && (
                  <span style={{ color: colors.textTertiary }}>None annotated</span>
                )}
                {h.pathogenic_variants.map((v, vi) => (
                  <VariantCard key={vi} v={v} />
                ))}
              </div>
            </div>
          ))}
        </div>
      )}
    </div>
  )
}

// Classification color: red only for ACMG's own definitive terms
// (Pathogenic / Likely pathogenic), backed by a real, freshly-checked ClinVar
// record -- never assumed just because a variant made it into this list.
// Amber for a real ClinVar entry that isn't a clean pathogenic call
// (Conflicting classifications, Uncertain significance, ...). Gray for a
// variant flagged purely on VEP's own HIGH-impact prediction, with no
// ClinVar entry at all.
function VariantCard({ v }: { v: VariantRecord }) {
  const hasClinvar = !!v.clinvar_id
  const color = v.is_definitive_pathogenic ? colors.statusPathogenic : hasClinvar ? colors.statusIntermediate : colors.textTertiary
  const bg = v.is_definitive_pathogenic ? colors.interpBg.pathogenic : hasClinvar ? colors.interpBg.intermediate : colors.borderSoftest

  return (
    <div style={{ padding: '8px 10px', borderRadius: 8, background: bg, border: `1px solid ${color}` }}>
      <div style={{ fontFamily: mono, fontWeight: 700, fontSize: 12 }}>
        {v.hgvsc ?? v.variant}
        {v.hgvsp && <span style={{ marginLeft: 8, fontWeight: 600 }}>{v.hgvsp}</span>}
      </div>
      <div style={{ fontSize: 11, color: colors.textSecondary, marginTop: 2 }}>
        {v.consequence_terms.join(', ')}
        {v.impact && ` · ${v.impact}`}
      </div>
      <div style={{ display: 'flex', alignItems: 'baseline', gap: 6, marginTop: 4, flexWrap: 'wrap' }}>
        <span style={{ fontWeight: 700, fontSize: 11.5, color }}>
          {hasClinvar ? v.clinvar_classification ?? 'ClinVar lookup failed' : 'No ClinVar entry (VEP-only)'}
        </span>
        {v.clinvar_id && (
          <span style={{ fontSize: 10.5, color: colors.textTertiary, fontFamily: mono }}>ClinVar {v.clinvar_id}</span>
        )}
      </div>
      {v.clinvar_review_status && (
        <div style={{ fontSize: 10.5, color: colors.textTertiary, marginTop: 2 }}>{v.clinvar_review_status}</div>
      )}
      {v.clinvar_condition && (
        <div style={{ fontSize: 10.5, color: colors.textSecondary, marginTop: 2 }}>{v.clinvar_condition}</div>
      )}
    </div>
  )
}
