import { colors, mono } from '../theme'
import { kivviApi } from '../api'
import { AlleleStrip } from '../components/AlleleStrip'
import { ReadPlot } from '../components/ReadPlot'
import { VariantTable } from '../components/VariantTable'
import { LoadState, Notes, SectionTitle } from '../components/ui'
import { cardStyle, pageStyle } from '../components/styles'
import { useLoad } from '../useLoad'

export function Kiv2Screen({ sample }: { sample: string }) {
  const { data, error } = useLoad(() => kivviApi.kiv2(sample), sample)
  if (!data) return <LoadState data={data} error={error} />

  return (
    <div style={pageStyle}>
      <div style={cardStyle}>
        <SectionTitle title="KIV2 copy number" sub="LPA kringle IV type 2 repeat, 5.5 kb unit. A short KIV2 allele is associated with higher cardiovascular risk." />
        <div style={{ fontSize: 30, fontWeight: 700, fontFamily: mono }}>{data.allele_cn || '—'}</div>
        {data.alleles.length > 1 && data.shortest_copy_number != null && (
          <div style={{ fontSize: 12.5, color: colors.textSecondary, marginTop: 4 }}>
            Shorter allele: {data.shortest_copy_number} copies. Kivvi defines no clinical cutoff, so none is applied.
          </div>
        )}
      </div>
      <Notes notes={data.notes} />

      {data.alleles.map((a) => (
        <div key={a.index} style={cardStyle}>
          <div style={{ display: 'flex', alignItems: 'baseline', gap: 12, flexWrap: 'wrap', marginBottom: 12 }}>
            <div style={{ fontSize: 12, fontWeight: 600, color: colors.textSecondary }}>Allele {a.index}</div>
            <div style={{ fontSize: 26, fontWeight: 700 }}>{a.copy_number}</div>
            <div style={{ fontSize: 12.5, color: colors.textTertiary }}>KIV2 copies</div>
            <div style={{ fontSize: 12.5, color: colors.textSecondary, marginLeft: 'auto' }}>
              <span style={{ fontFamily: mono }}>{a.supporting_reads}</span> supporting reads
            </div>
          </div>
          <AlleleStrip units={a.units} />
        </div>
      ))}

      {data.partial_alleles.length > 0 && (
        <div style={cardStyle}>
          <SectionTitle title="Partially assembled alleles" sub="Lower bound only — these did not span the whole repeat." />
          {data.partial_alleles.map((p, i) => (
            <div key={i} style={{ marginBottom: 10 }}>
              <div style={{ fontSize: 12.5, marginBottom: 6 }}>
                ≥<span style={{ fontFamily: mono }}>{p.lower_bound}</span> copies
              </div>
              <AlleleStrip units={p.units} />
            </div>
          ))}
        </div>
      )}

      <ReadPlot target="kiv2" sample={sample} available={data.svg_available} />
      <VariantTable variants={data.variants} />
    </div>
  )
}
