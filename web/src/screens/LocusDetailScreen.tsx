import { useEffect, useState } from 'react'
import { api, type LocusDetail } from '../api'
import { colors, interpBg, interpColor, mono } from '../theme'
import { Card, CardHeading } from '../components/Card'
import { ExpressionEffectsCard } from '../components/ExpressionEffectsCard'
import { Pileup } from '../charts/Pileup'
import { AlleleHistogram } from '../charts/AlleleHistogram'
import { MethylationPanel } from '../charts/MethylationPanel'

function fmtRange(min: number | null, max: number | null, prefix = ''): string {
  if (min == null && max == null) return '—'
  if (min != null && max == null) return `≥${min}`
  if (min == null && max != null) return `≤${max}`
  return `${prefix}${min}–${max}`
}

// Some loci (RFC1, the FAME genes, ...) have several candidate motifs at the same
// tract, where which motif expanded matters as much as the length -- call it out
// whenever the actual dominant motif differs from the locus's reference motif.
function fmtAllele(count: number | null, motif: string | null, referenceMotif: string): string {
  if (count == null) return '—'
  if (motif && motif !== referenceMotif) return `${count} repeats (${motif})`
  return `${count} repeats`
}

export function LocusDetailScreen({
  mode, sample, locusId,
}: {
  mode: string; sample: string; locusId: string
}) {
  const [locus, setLocus] = useState<LocusDetail | null>(null)
  const [notFound, setNotFound] = useState(false)

  useEffect(() => {
    setLocus(null)
    setNotFound(false)
    if (!locusId) return
    api.locus(locusId, mode, sample).then(setLocus).catch(() => setNotFound(true))
  }, [locusId, mode, sample])

  if (!locusId || notFound) {
    return (
      <div style={{ padding: 28, color: colors.textSecondary }}>
        Select a locus from the Overview table to see its detail view.
      </div>
    )
  }
  if (!locus) return <div style={{ padding: 28 }}>Loading…</div>

  return (
    <div style={{ padding: 28, display: 'grid', gridTemplateColumns: '1fr 320px', gap: 20, alignItems: 'start' }}>
      <div style={{ display: 'flex', flexDirection: 'column', gap: 20, minWidth: 0 }}>
        <Card>
          <div style={{ display: 'flex', alignItems: 'baseline', gap: 12 }}>
            <div style={{ fontSize: 22, fontWeight: 700 }}>{locus.gene}</div>
            <div style={{ fontFamily: mono, fontSize: 13, color: colors.textSecondary }}>{locus.region}</div>
          </div>
          <div style={{ fontSize: 13.5, color: colors.textQuaternary, marginTop: 4 }}>{locus.disease}</div>
          <div style={{ display: 'flex', gap: 28, marginTop: 14, flexWrap: 'wrap' }}>
            <Field label="Reference motif" value={locus.motif} />
            <Field label="Allele 1" value={fmtAllele(locus.a1, locus.a1_motif, locus.motif)} />
            <Field label="Allele 2" value={fmtAllele(locus.a2, locus.a2_motif, locus.motif)} />
            <Field label="Zygosity" value={locus.zygosity ?? '—'} mono={false} />
            <Field label="Read depth" value={locus.depth != null ? `${locus.depth}x` : '—'} mono={false} />
          </div>
        </Card>

        <Card>
          <CardHeading title="Read pileup" sub="Per-read repeat count, sorted · colored by allele cluster · estimated from SD/ALLR" />
          <div style={{ width: '100%', overflowX: 'auto' }}>
            <Pileup reads={locus.pileup_reads} />
          </div>
          <div style={{ display: 'flex', gap: 18, marginTop: 10, fontSize: 12, color: colors.textSecondary, flexWrap: 'wrap' }}>
            <Legend color={colors.cluster1} label="Allele 1 cluster" />
            <Legend color={colors.cluster2} label="Allele 2 cluster" />
          </div>
        </Card>

        <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 20 }}>
          <Card>
            <CardHeading title="Allele size distribution" sub="Estimated read counts by repeat length" />
            <AlleleHistogram reads={locus.pileup_reads} a1={locus.a1} a2={locus.a2} />
          </Card>
          <Card>
            <CardHeading title="CpG methylation" sub="Real per-allele mean 5mCpG level (TRGT AM field)" />
            <MethylationPanel a1Methylation={locus.a1_methylation} a2Methylation={locus.a2_methylation} genotyped={locus.genotyped} />
          </Card>
        </div>

        {locus.expression_effects && <ExpressionEffectsCard data={locus.expression_effects} motif={locus.motif} />}
      </div>

      <div style={{ display: 'flex', flexDirection: 'column', gap: 16 }}>
        <Card>
          <SectionLabel text="Interpretation" />
          <div style={{ marginTop: 12, padding: '12px 14px', borderRadius: 8, background: interpBg(locus.status) }}>
            <div style={{ fontSize: 14.5, fontWeight: 700, color: interpColor(locus.status) }}>{locus.status}</div>
            <div style={{ fontSize: 13, color: colors.textBody, marginTop: 4, lineHeight: 1.5, wordBreak: 'break-word' }}>{locus.note}</div>
          </div>
          <div style={{ marginTop: 16, display: 'flex', flexDirection: 'column', gap: 10, fontSize: 13 }}>
            <Row label="Classification" value={locus.status} />
            <Row label="Inheritance pattern" value={locus.inheritance || '—'} />
            <Row label="Reads spanning locus" value={locus.depth != null ? String(locus.depth) : '—'} />
            <Row label="Caller" value={locus.caller} />
          </div>
        </Card>

        <Card>
          <SectionLabel text="Reference ranges" />
          <div style={{ marginTop: 12, display: 'flex', flexDirection: 'column', gap: 8, fontSize: 13 }}>
            <Row label="Normal" value={fmtRange(null, locus.benign_max)} mono />
            <Row label="Intermediate" value={fmtRange(locus.intermediate_min, locus.intermediate_max)} mono />
            <Row label="Pathogenic" value={fmtRange(locus.pathogenic_min, null)} mono />
          </div>
        </Card>

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

function Field({ label, value, mono: useMono = true }: { label: string; value: string; mono?: boolean }) {
  return (
    <div>
      <div style={{ fontSize: 11, fontWeight: 600, color: colors.textSecondary, textTransform: 'uppercase', letterSpacing: '0.04em' }}>{label}</div>
      <div style={{ fontFamily: useMono ? mono : undefined, fontWeight: 700, fontSize: 15, marginTop: 2 }}>{value}</div>
    </div>
  )
}

function Legend({ color, label }: { color: string; label: string }) {
  return (
    <div style={{ display: 'flex', alignItems: 'center', gap: 6 }}>
      <div style={{ width: 10, height: 10, borderRadius: 2, background: color }} />
      {label}
    </div>
  )
}

function SectionLabel({ text }: { text: string }) {
  return (
    <div style={{ fontSize: 13, fontWeight: 600, color: colors.textSecondary, textTransform: 'uppercase', letterSpacing: '0.04em' }}>{text}</div>
  )
}

function Row({ label, value, mono: useMono = false }: { label: string; value: string; mono?: boolean }) {
  return (
    <div style={{ display: 'flex', justifyContent: 'space-between' }}>
      <span style={{ color: colors.textSecondary }}>{label}</span>
      <span style={{ fontWeight: useMono ? 400 : 600, fontFamily: useMono ? mono : undefined }}>{value}</span>
    </div>
  )
}
