import { useEffect, useState } from 'react'
import { colors, mono } from '../../theme'
import { StatCard } from '../../components/StatCard'
import { paraphaseApi, type JsonRegionDetail } from '../api'
import { FLAG_COLOR, FLAG_LABEL } from '../flags'
import { headerRowStyle } from '../styles'
import { Card, Chip } from '../ui'

const GRID = '1.8fr 1.2fr 1.2fr 0.9fr 0.9fr 0.7fr 0.7fr'

const fmtBoundary = (b: (number | null)[] | null) =>
  b && b.every((x) => x != null) ? `${b[0]?.toLocaleString()}–${b[1]?.toLocaleString()}` : '—'

export function ParaphaseCallsScreen({ sample, region }: { sample: string; region: string }) {
  const [data, setData] = useState<JsonRegionDetail | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [openHap, setOpenHap] = useState<string | null>(null)

  useEffect(() => {
    setData(null)
    setError(null)
    setOpenHap(null)
    paraphaseApi.jsonRegion(region, sample).then(setData).catch((e) => setError(String(e)))
  }, [sample, region])

  if (error) return <div style={{ padding: 28, color: colors.statusPathogenic }}>Failed to load: {error}</div>
  if (!data) return <div style={{ padding: 28, color: colors.textSecondary }}>Loading…</div>

  const s = data.summary
  const fields = Object.entries(data.region_fields)

  return (
    <div style={{ padding: 28, display: 'flex', flexDirection: 'column', gap: 20 }}>
      <div>
        <div style={{ display: 'flex', alignItems: 'baseline', gap: 12 }}>
          <div style={{ fontSize: 18, fontWeight: 700, fontFamily: mono }}>{s.region}</div>
          <div style={{ fontSize: 13, color: colors.textSecondary }}>{s.genes.join(', ')}</div>
        </div>
        {s.phase_region && (
          <div style={{ fontSize: 11.5, color: colors.textTertiary, fontFamily: mono, marginTop: 3 }}>
            Phased over {s.phase_region}
          </div>
        )}
      </div>

      {s.flags.length > 0 && (
        <div style={{ display: 'flex', flexDirection: 'column', gap: 8 }}>
          {s.flags.map((f) => (
            <div
              key={f.kind}
              style={{
                display: 'flex',
                gap: 10,
                alignItems: 'center',
                padding: '10px 14px',
                borderRadius: 10,
                background: 'white',
                border: `1px solid ${FLAG_COLOR[f.kind]}`,
                fontSize: 12.5,
              }}
            >
              <Chip color={FLAG_COLOR[f.kind]}>{FLAG_LABEL[f.kind]}</Chip>
              <span>{f.text}</span>
            </div>
          ))}
        </div>
      )}

      <div style={{ display: 'grid', gridTemplateColumns: 'repeat(4,1fr)', gap: 16 }}>
        <StatCard label="Total copies" value={s.total_cn != null ? String(s.total_cn) : '—'} sub="Gene plus paralog or pseudogene" />
        <StatCard label="Gene copies" value={s.gene_cn != null ? String(s.gene_cn) : '—'} sub={s.gene_cn == null ? 'Not reported for this region' : 'Excluding paralogs'} />
        <StatCard
          label="Region depth"
          value={s.depth_median != null ? `${s.depth_median}x` : '—'}
          sub={
            s.depth_median == null
              ? 'No depth reported'
              : `80th percentile ${s.depth_p80 ?? '—'}x${s.genome_depth ? ` · genome ${s.genome_depth}x` : ''}`
          }
        />
        <StatCard
          label="Haplotypes"
          value={String(s.n_haplotypes)}
          sub={s.two_copy_haplotypes ? `${s.two_copy_haplotypes} inferred at two copies` : 'As phased'}
        />
      </div>

      {s.headline && (
        <div style={{ fontSize: 13, color: colors.textBody }}>
          <span style={{ fontWeight: 600 }}>Headline: </span>
          {s.headline}
        </div>
      )}

      <Card title="Haplotypes" sub="Each phased copy of the gene family. Variants are shown as Paraphase recorded them; use its VCF for accurate small-variant calls.">
        {data.haplotypes.length === 0 ? (
          <div style={{ padding: '40px 20px', textAlign: 'center', color: colors.textSecondary, fontSize: 13.5 }}>
            Paraphase phased no haplotypes in this region.
          </div>
        ) : (
          <>
            <div style={{ display: 'grid', gridTemplateColumns: GRID, gap: 10, ...headerRowStyle }}>
              <div>Haplotype</div>
              <div>Assigned to</div>
              <div>Boundary</div>
              <div>Truncated</div>
              <div>Haplogroup</div>
              <div>Reads</div>
              <div>Variants</div>
            </div>
            {data.haplotypes.map((h) => (
              <div key={h.name} style={{ borderBottom: `1px solid ${colors.borderSoftest}` }}>
                <div
                  onClick={h.n_variants ? () => setOpenHap(openHap === h.name ? null : h.name) : undefined}
                  style={{
                    display: 'grid',
                    gridTemplateColumns: GRID,
                    gap: 10,
                    padding: '12px 20px',
                    alignItems: 'center',
                    fontSize: 12.5,
                    cursor: h.n_variants ? 'pointer' : 'default',
                  }}
                >
                  <div style={{ fontFamily: mono, fontWeight: 600 }}>
                    {h.name}
                    {h.two_copy && (
                      <span style={{ marginLeft: 6, color: colors.textTertiary, fontWeight: 400 }} title="Inferred to exist in two identical copies">
                        ×2
                      </span>
                    )}
                  </div>
                  <div style={{ color: colors.textBody }}>{h.kind ?? h.group ?? '—'}</div>
                  <div style={{ fontFamily: mono, color: colors.textSecondary }}>{fmtBoundary(h.boundary)}</div>
                  <div style={{ fontFamily: mono }}>{h.truncated.length ? h.truncated.join(', ') : '—'}</div>
                  <div style={{ fontFamily: mono }}>{h.haplogroup ?? '—'}</div>
                  <div style={{ fontFamily: mono }}>{h.unique_reads ?? '—'}</div>
                  <div style={{ fontFamily: mono }}>{h.n_variants}</div>
                </div>
                {openHap === h.name && (
                  <div style={{ padding: '0 20px 14px', display: 'flex', flexWrap: 'wrap', gap: 6 }}>
                    {h.variants.map((v) => (
                      <span key={v} style={{ fontFamily: mono, fontSize: 11, background: colors.chipBg, borderRadius: 6, padding: '2px 7px' }}>
                        {v}
                      </span>
                    ))}
                    {h.n_variants > h.variants.length && (
                      <span style={{ fontSize: 11, color: colors.textTertiary, alignSelf: 'center' }}>
                        +{h.n_variants - h.variants.length} more in the Paraphase VCF
                      </span>
                    )}
                  </div>
                )}
              </div>
            ))}
          </>
        )}
      </Card>

      {data.fusions.length > 0 && (
        <Card title="Fusion haplotypes" sub="Deletions or duplications from unequal crossing over between paralogs, called from the flanking sequence of phased haplotypes.">
          {data.fusions.map((f) => (
            <div key={f.name} style={{ display: 'grid', gridTemplateColumns: '2fr 1fr 2fr', gap: 10, padding: '12px 20px', fontSize: 12.5, borderBottom: `1px solid ${colors.borderSoftest}` }}>
              <div style={{ fontFamily: mono, fontWeight: 600 }}>{f.name}</div>
              <div>{f.type ?? '—'}</div>
              <div style={{ fontFamily: mono, color: colors.textSecondary }}>
                {f.breakpoint != null ? JSON.stringify(f.breakpoint) : '—'}
              </div>
            </div>
          ))}
        </Card>
      )}

      {data.alleles && data.alleles.length > 0 && (
        <Card title="Alleles" sub="Haplotypes linked into alleles by read-level evidence (tandem segmental duplications only).">
          {data.alleles.map((a, i) => (
            <div key={i} style={{ padding: '12px 20px', fontSize: 12.5, borderBottom: `1px solid ${colors.borderSoftest}`, fontFamily: mono }}>
              <span style={{ color: colors.textTertiary }}>Allele {i + 1}: </span>
              {a.join(' + ')}
            </div>
          ))}
        </Card>
      )}

      {fields.length > 0 && (
        <Card title="Region-specific fields" sub="Scalar values Paraphase reports only for this kind of region, as written.">
          <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 0 }}>
            {fields.map(([k, v]) => (
              <div key={k} style={{ display: 'flex', justifyContent: 'space-between', gap: 12, padding: '10px 20px', fontSize: 12.5, borderBottom: `1px solid ${colors.borderSoftest}` }}>
                <span style={{ color: colors.textSecondary, fontFamily: mono }}>{k}</span>
                <span style={{ fontFamily: mono, fontWeight: 600 }}>{v == null ? 'null' : String(v)}</span>
              </div>
            ))}
          </div>
        </Card>
      )}
    </div>
  )
}
