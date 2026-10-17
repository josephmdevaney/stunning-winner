import { callBg, callColor, colors, fmtPct, mono } from '../theme'
import { kivviApi, type D4Z4Allele } from '../api'
import { AlleleStrip } from '../components/AlleleStrip'
import { ReadPlot } from '../components/ReadPlot'
import { VariantTable } from '../components/VariantTable'
import { LoadState, Notes, SectionTitle, StatusPill } from '../components/ui'
import { cardStyle, pageStyle } from '../components/styles'
import { useLoad } from '../useLoad'

export function D4Z4Screen({ sample }: { sample: string }) {
  const { data, error } = useLoad(() => kivviApi.d4z4(sample), sample)
  if (!data) return <LoadState data={data} error={error} />

  const { summary } = data
  return (
    <div style={pageStyle}>
      <div style={{ ...cardStyle, background: callBg(summary.status), borderColor: callColor(summary.status) }}>
        <div style={{ fontSize: 17, fontWeight: 600, color: callColor(summary.status) }}>{summary.label}</div>
        <div style={{ fontSize: 13, color: colors.textBody, marginTop: 6 }}>
          {summary.complete_alleles} complete · {summary.partial_alleles} partial allele{summary.partial_alleles === 1 ? '' : 's'} · overall methylation{' '}
          <span style={{ fontFamily: mono }}>{fmtPct(data.overall_methylation)}</span>
        </div>
      </div>
      <Notes notes={summary.notes} />

      {data.alleles.map((a, i) => (
        <AlleleCard key={i} index={i + 1} allele={a} />
      ))}

      <div style={cardStyle}>
        <SectionTitle title="In-cis D4Z4 arrays" />
        {data.in_cis_arrays.length === 0 ? (
          <div style={{ fontSize: 13, color: colors.textSecondary }}>No in-cis duplication of D4Z4 arrays was called.</div>
        ) : (
          <pre style={{ margin: 0, fontFamily: mono, fontSize: 12, whiteSpace: 'pre-wrap' }}>{JSON.stringify(data.in_cis_arrays, null, 2)}</pre>
        )}
      </div>

      <div style={{ ...cardStyle, fontSize: 12.5, color: colors.textBody, lineHeight: 1.6 }}>
        <SectionTitle title="How alleles are called" />
        FSHD1 (about 95% of FSHD) results from contraction of D4Z4 to 1–10 units on a chr4 qA haplotype with an intact PolyA signal. Contracted
        alleles on a <b>qB</b> or <b>qADisruptedPolyA</b> haplotype, or on the chr4 <b>Group2.1</b> (4A166) haplotype, are not pathogenic. Partially
        assembled alleles report a lower bound (≥N). Methylation is shown per allele (last three units) and per unit; Kivvi defines no
        hypomethylation cutoff, so none is applied — low overall methylation is the signature of FSHD2.
      </div>

      <ReadPlot target="d4z4" sample={sample} available={data.svg_available} />
      <VariantTable variants={data.variants} />
    </div>
  )
}

function AlleleCard({ index, allele }: { index: number; allele: D4Z4Allele }) {
  const chrom = allele.upstream_group ? `${allele.chromosome} · ${allele.upstream_group}` : allele.chromosome
  const facts: [string, string][] = [
    ['Chromosome', chrom],
    ['Distal haplotype', allele.distal_haplotype],
    ['Methylation', fmtPct(allele.methylation)],
    ['Assembly', allele.allele_type],
  ]
  return (
    <div style={cardStyle}>
      <div style={{ display: 'flex', alignItems: 'baseline', gap: 12, flexWrap: 'wrap' }}>
        <div style={{ fontSize: 12, fontWeight: 600, color: colors.textSecondary }}>Allele {index}</div>
        <div style={{ fontSize: 26, fontWeight: 700 }}>{allele.size_label}</div>
        <div style={{ fontSize: 12.5, color: colors.textTertiary }}>units</div>
        <StatusPill status={allele.interpretation.status} label={allele.interpretation.label} />
        {allele.ending_in_qal && <span style={{ fontSize: 11.5, color: colors.textSecondary, background: colors.chipBg, padding: '2px 8px', borderRadius: 999 }}>A-L allele</span>}
      </div>
      <div style={{ fontSize: 12.5, color: colors.textBody, margin: '6px 0 14px' }}>{allele.interpretation.reason}</div>
      <div style={{ display: 'flex', gap: 28, flexWrap: 'wrap', marginBottom: 14 }}>
        {facts.map(([k, v]) => (
          <div key={k}>
            <div style={{ fontSize: 10.5, fontWeight: 600, color: colors.textSecondary, textTransform: 'uppercase', letterSpacing: '0.04em' }}>{k}</div>
            <div style={{ fontFamily: mono, fontSize: 13, marginTop: 2 }}>{v}</div>
          </div>
        ))}
      </div>
      <AlleleStrip units={allele.units} methylation={allele.per_unit_methylation} />
    </div>
  )
}
