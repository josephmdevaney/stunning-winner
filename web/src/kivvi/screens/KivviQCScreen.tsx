import { colors, fmtNum } from '../theme'
import { kivviApi, type TargetQC } from '../api'
import { LoadState, Notes, SectionTitle } from '../components/ui'
import { cardStyle, pageStyle } from '../components/styles'
import { useLoad } from '../useLoad'

export function KivviQCScreen({ sample }: { sample: string }) {
  const { data, error } = useLoad(() => kivviApi.qc(sample), sample)
  if (!data) return <LoadState data={data} error={error} />

  return (
    <div style={pageStyle}>
      {data.input_note && <Notes notes={[`${data.input_note}. Kivvi reports the genome depth only from a whole-genome BAM.`]} />}
      {data.target_qc.map((t) => (
        <TargetCard key={t.target} qc={t} />
      ))}
    </div>
  )
}

function Tile({ label, value, sub }: { label: string; value: string; sub?: string }) {
  return (
    <div style={{ border: `1px solid ${colors.border}`, borderRadius: 10, padding: '14px 16px' }}>
      <div style={{ fontSize: 12, fontWeight: 600, color: colors.textSecondary }}>{label}</div>
      <div style={{ fontSize: 24, fontWeight: 700, marginTop: 6 }}>{value}</div>
      {sub && <div style={{ fontSize: 12, color: colors.textTertiary, marginTop: 4 }}>{sub}</div>}
    </div>
  )
}

function TargetCard({ qc }: { qc: TargetQC }) {
  const resolvedAll = qc.alleles_resolved >= qc.alleles_expected
  return (
    <div style={cardStyle}>
      <SectionTitle title={qc.label} />
      <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))', gap: 14 }}>
        <Tile label="Repeat depth" value={fmtNum(qc.repeat_depth, '×')} sub="After realignment to one repeat unit" />
        <Tile label="Genome depth" value={qc.genome_depth == null ? 'Unavailable' : fmtNum(qc.genome_depth, '×')} sub={qc.genome_depth == null ? 'Needs a whole-genome BAM' : undefined} />
        <Tile
          label="Alleles resolved"
          value={`${qc.alleles_resolved} / ${qc.alleles_expected}`}
          sub={resolvedAll ? (qc.partial_alleles ? `${qc.partial_alleles} partial` : 'All fully assembled') : 'Fewer than expected'}
        />
        <Tile label="Reads linking units" value={fmtNum(qc.reads_linking_units)} />
        {qc.median_read_length != null && <Tile label="Median read length" value={fmtNum(qc.median_read_length, ' bp')} />}
        {qc.per_allele_depth != null && <Tile label="Depth per allele" value={fmtNum(qc.per_allele_depth, '×')} sub="D4Z4 is GC-rich; can differ from genome depth" />}
        <Tile label="Small variants" value={fmtNum(qc.variants_total)} sub={`${qc.variants_on_assembled} on assembled alleles`} />
      </div>
    </div>
  )
}
