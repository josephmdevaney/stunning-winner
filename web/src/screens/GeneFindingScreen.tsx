import { useEffect, useState } from 'react'
import { api, type GeneFindingDetail } from '../api'
import { colors, interpBg, interpColor, mono } from '../theme'
import { Card, CardHeading } from '../components/Card'

export function GeneFindingScreen({
  mode, sample, gene,
}: {
  mode: string; sample: string; gene: string
}) {
  const [data, setData] = useState<GeneFindingDetail | null>(null)
  const [notFound, setNotFound] = useState(false)

  useEffect(() => {
    setData(null)
    setNotFound(false)
    if (!gene) return
    api.geneFinding(gene, mode, sample).then(setData).catch(() => setNotFound(true))
  }, [gene, mode, sample])

  if (!gene || notFound) {
    return (
      <div style={{ padding: 28, color: colors.textSecondary }}>
        Select a gene finding from the Overview table to see its detail view.
      </div>
    )
  }
  if (!data) return <div style={{ padding: 28 }}>Loading…</div>

  const rawFieldEntries = Object.entries(data.raw_fields)

  return (
    <div style={{ padding: 28, display: 'grid', gridTemplateColumns: '1fr 320px', gap: 20, alignItems: 'start' }}>
      <div style={{ display: 'flex', flexDirection: 'column', gap: 20, minWidth: 0 }}>
        <Card>
          <div style={{ display: 'flex', alignItems: 'baseline', gap: 12 }}>
            <div style={{ fontSize: 22, fontWeight: 700 }}>{data.gene}</div>
            <div style={{ fontFamily: mono, fontSize: 13, color: colors.textSecondary }}>{data.region}</div>
          </div>
          <div style={{ fontSize: 13.5, color: colors.textQuaternary, marginTop: 4 }}>{data.disease}</div>
          <div style={{ fontSize: 12.5, color: colors.textTertiary, marginTop: 10, lineHeight: 1.5 }}>
            Detected via <strong>{data.tool}</strong> -- {data.mechanism}. Not a tandem repeat, so this finding has
            no read-pileup or allele-length view; everything shown here is read directly off a field the tool
            itself asserted.
          </div>
        </Card>

        <Card>
          <CardHeading
            title={`Named variants (${data.variants.length})`}
            sub={
              data.variants.length
                ? 'Every variant this gene’s pathogenic-variant annotation flagged, deduplicated across redundant haplotype-phase calls'
                : 'No named pathogenic variants flagged for this gene'
            }
          />
          {data.variants.length > 0 ? (
            <div style={{ display: 'flex', flexDirection: 'column', gap: 6, marginTop: 4 }}>
              {data.variants.map((v) => (
                <div
                  key={v}
                  style={{
                    fontFamily: mono,
                    fontSize: 12.5,
                    padding: '8px 12px',
                    borderRadius: 6,
                    background: 'oklch(0.975 0.005 220)',
                    border: `1px solid ${colors.borderSoftest}`,
                    wordBreak: 'break-all',
                  }}
                >
                  {v}
                </div>
              ))}
            </div>
          ) : (
            <div style={{ fontSize: 13, color: colors.textTertiary }}>—</div>
          )}
        </Card>

        <Card>
          <CardHeading
            title="Raw tool fields"
            sub="Verbatim fields the caller reported for this gene -- nothing here is re-derived or summarized"
          />
          {rawFieldEntries.length > 0 ? (
            <div style={{ display: 'flex', flexDirection: 'column', gap: 8, fontSize: 13 }}>
              {rawFieldEntries.map(([k, v]) => (
                <div key={k} style={{ display: 'flex', gap: 12 }}>
                  <div style={{ fontFamily: mono, color: colors.textSecondary, minWidth: 140, flexShrink: 0 }}>{k}</div>
                  <div style={{ fontFamily: mono, wordBreak: 'break-all' }}>{v}</div>
                </div>
              ))}
            </div>
          ) : (
            <div style={{ fontSize: 13, color: colors.textTertiary }}>—</div>
          )}
        </Card>
      </div>

      <div style={{ display: 'flex', flexDirection: 'column', gap: 16 }}>
        <Card>
          <SectionLabel text="Interpretation" />
          <div style={{ marginTop: 12, padding: '12px 14px', borderRadius: 8, background: interpBg(data.status) }}>
            <div style={{ fontSize: 14.5, fontWeight: 700, color: interpColor(data.status) }}>{data.status}</div>
            <div style={{ fontSize: 13, color: colors.textBody, marginTop: 4, lineHeight: 1.5, wordBreak: 'break-word' }}>{data.note}</div>
          </div>
          <div style={{ marginTop: 16, display: 'flex', flexDirection: 'column', gap: 10, fontSize: 13 }}>
            <Row label="Classification" value={data.status} />
            <Row label="Zygosity" value={data.zygosity ?? '—'} />
            <Row label="Caller" value={data.caller} />
          </div>
        </Card>
      </div>
    </div>
  )
}

function SectionLabel({ text }: { text: string }) {
  return (
    <div style={{ fontSize: 13, fontWeight: 600, color: colors.textSecondary, textTransform: 'uppercase', letterSpacing: '0.04em' }}>{text}</div>
  )
}

function Row({ label, value }: { label: string; value: string }) {
  return (
    <div style={{ display: 'flex', justifyContent: 'space-between' }}>
      <span style={{ color: colors.textSecondary }}>{label}</span>
      <span style={{ fontWeight: 600 }}>{value}</span>
    </div>
  )
}
