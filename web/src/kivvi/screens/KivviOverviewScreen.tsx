import { colors, fmtPct, mono } from '../theme'
import { kivviApi, type OverviewResponse } from '../api'
import type { KivviScreen } from '../KivviSidebar'
import { LoadState, Notes, SectionTitle, StatusPill } from '../components/ui'
import { cardStyle, pageStyle } from '../components/styles'
import { useLoad } from '../useLoad'

export function KivviOverviewScreen({ sample, onNavigate }: { sample: string; onNavigate: (s: KivviScreen) => void }) {
  const { data, error } = useLoad(() => kivviApi.overview(sample), sample)
  if (!data) return <LoadState data={data} error={error} />

  return (
    <div style={pageStyle}>
      <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(420px, 1fr))', gap: 20, alignItems: 'start' }}>
        <D4Z4Card data={data} onOpen={() => onNavigate('d4z4')} />
        <Kiv2Card data={data} onOpen={() => onNavigate('kiv2')} />
      </div>
    </div>
  )
}

function OpenLink({ onClick }: { onClick: () => void }) {
  return (
    <button
      onClick={onClick}
      style={{ border: `1px solid ${colors.borderInput}`, background: 'white', borderRadius: 6, padding: '6px 12px', fontSize: 12.5, fontWeight: 600, cursor: 'pointer', flexShrink: 0 }}
    >
      View detail →
    </button>
  )
}

function Missing({ label }: { label: string }) {
  return (
    <div style={cardStyle}>
      <SectionTitle title={label} />
      <div style={{ color: colors.textSecondary, fontSize: 13 }}>No Kivvi output for this target.</div>
    </div>
  )
}

function D4Z4Card({ data, onOpen }: { data: OverviewResponse; onOpen: () => void }) {
  const d = data.d4z4
  if (!d) return <Missing label="D4Z4 (FSHD)" />
  const cols = '0.7fr 1.1fr 1.3fr 0.8fr 1.4fr'
  return (
    <div style={cardStyle}>
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', gap: 12 }}>
        <SectionTitle title="D4Z4 (FSHD)" sub="Contraction to 1–10 units on a chr4 qA haplotype with an intact PolyA signal is FSHD1-range." />
        <OpenLink onClick={onOpen} />
      </div>
      <div style={{ margin: '4px 0 14px' }}>
        <StatusPill status={d.summary.status} label={d.summary.label} />
      </div>

      <div style={{ display: 'grid', gridTemplateColumns: cols, gap: 10, padding: '8px 0', fontSize: 10.5, fontWeight: 600, color: colors.textSecondary, letterSpacing: '0.04em', textTransform: 'uppercase', borderBottom: `1px solid ${colors.borderSofter}` }}>
        <div>Units</div>
        <div>Chromosome</div>
        <div>Distal haplotype</div>
        <div>Methyl.</div>
        <div>Call</div>
      </div>
      {d.alleles.map((a, i) => (
        <div key={i} style={{ display: 'grid', gridTemplateColumns: cols, gap: 10, padding: '9px 0', fontSize: 12.5, alignItems: 'center', borderBottom: `1px solid ${colors.borderSoftest}` }}>
          <div style={{ fontFamily: mono, fontWeight: 600 }}>{a.size_label}</div>
          <div style={{ fontFamily: mono }}>
            {a.chromosome}
            {a.upstream_group ? ` · ${a.upstream_group}` : ''}
          </div>
          <div style={{ fontFamily: mono }}>{a.distal_haplotype}</div>
          <div style={{ fontFamily: mono }}>{fmtPct(a.methylation)}</div>
          <div>
            <StatusPill status={a.interpretation.status} label={a.interpretation.label} />
          </div>
        </div>
      ))}
      <div style={{ fontSize: 12.5, color: colors.textSecondary, marginTop: 12 }}>
        Overall methylation across D4Z4 units: <span style={{ fontFamily: mono }}>{fmtPct(d.overall_methylation)}</span>
      </div>
      <div style={{ marginTop: 12 }}>
        <Notes notes={d.summary.notes} />
      </div>
    </div>
  )
}

function Kiv2Card({ data, onOpen }: { data: OverviewResponse; onOpen: () => void }) {
  const k = data.kiv2
  if (!k) return <Missing label="LPA KIV2" />
  return (
    <div style={cardStyle}>
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', gap: 12 }}>
        <SectionTitle title="LPA KIV2" sub="Kringle IV type 2 copy number per allele. A short allele is associated with higher cardiovascular risk." />
        <OpenLink onClick={onOpen} />
      </div>
      {k.alleles.length > 0 ? (
        <div style={{ display: 'flex', gap: 16, margin: '8px 0 6px' }}>
          {k.alleles.map((a) => (
            <div key={a.index} style={{ flex: 1, border: `1px solid ${colors.border}`, borderRadius: 10, padding: '14px 16px', background: a.copy_number === k.shortest_copy_number && k.alleles.length > 1 ? colors.interpBg.intermediate : undefined }}>
              <div style={{ fontSize: 11.5, fontWeight: 600, color: colors.textSecondary }}>Allele {a.index}</div>
              <div style={{ fontSize: 28, fontWeight: 700, marginTop: 4 }}>{a.copy_number}</div>
              <div style={{ fontSize: 12, color: colors.textTertiary }}>KIV2 copies</div>
            </div>
          ))}
        </div>
      ) : (
        <div style={{ color: colors.textSecondary, fontSize: 13 }}>No KIV2 allele was assembled.</div>
      )}
      {k.alleles.length > 1 && k.shortest_copy_number != null && (
        <div style={{ fontSize: 12.5, color: colors.textSecondary, marginTop: 8 }}>
          Shorter allele: <span style={{ fontFamily: mono }}>{k.shortest_copy_number} copies</span>
        </div>
      )}
      <div style={{ marginTop: 12 }}>
        <Notes notes={k.notes} />
      </div>
    </div>
  )
}
