import type { ReactNode } from 'react'
import { colors, fmtPct, mono } from '../theme'
import { kivviApi, type Fshd2Allele, type Fshd2Evidence, type SizeRange } from '../api'
import { LoadState, SectionTitle } from '../components/ui'
import { cardStyle, pageStyle } from '../components/styles'
import { useLoad } from '../useLoad'

const SIZE_RANGE_LABEL: Record<SizeRange, string> = {
  contracted: '1–10 units',
  fshd2_range: '11–20 units',
  above: 'Above 20',
  unresolved: 'Not resolved',
}

const SIZE_RANGE_HELP: Record<SizeRange, string> = {
  contracted: 'The FSHD1 contraction range (1–10 units).',
  fshd2_range: 'The range the Kivvi preprint calls commonly indicated for FSHD2.',
  above: 'Above both ranges.',
  unresolved: 'A partially assembled allele: its lower bound does not say which range it is in.',
}

// Permissive alleles in the accent colour, the rest in grey; neither means "normal" or "abnormal".
const PERMISSIVE = colors.accentButtonBg
const OTHER = colors.statusUnknown

export function FSHD2Screen({ sample }: { sample: string }) {
  const { data, error } = useLoad(() => kivviApi.d4z4(sample), sample)
  if (!data) return <LoadState data={data} error={error} />
  const { fshd2 } = data
  const { summary } = fshd2

  return (
    <div style={pageStyle}>
      <div style={cardStyle}>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', gap: 16 }}>
          <div>
            <div style={{ fontSize: 17, fontWeight: 600 }}>FSHD2 evidence</div>
            <div style={{ fontSize: 13, color: colors.textBody, marginTop: 6, maxWidth: 760, lineHeight: 1.55 }}>
              What Kivvi reports that bears on FSHD2, side by side. No call is made: Kivvi defines no FSHD2 rule and no methylation cut-off, and its preprint validated it on FSHD1 and control samples, not on confirmed FSHD2 cases.
            </div>
          </div>
          <span style={{ padding: '3px 10px', borderRadius: 999, fontSize: 11.5, fontWeight: 600, color: colors.textSecondary, background: colors.chipBg, whiteSpace: 'nowrap' }}>No call made</span>
        </div>
      </div>

      <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(230px, 1fr))', gap: 16 }}>
        <Stat
          label="Sample methylation"
          value={fmtPct(fshd2.sample_methylation)}
          sub="Methylated reads over all reads across every D4Z4 unit (Kivvi's methylation_all_sites). Kivvi notes it is low in FSHD2 samples."
        />
        <Stat
          label="Permissive alleles"
          value={String(summary.permissive_alleles)}
          sub="chr4, qAIntactPolyA and not the 4A166 haplotype: the kind that can express DUX4."
        />
        <Stat
          label="Permissive, 11–20 units"
          value={String(summary.permissive_fshd2_range)}
          sub={`The range the preprint calls commonly indicated for FSHD2.${summary.permissive_contracted ? ` ${summary.permissive_contracted} permissive allele${summary.permissive_contracted === 1 ? ' is' : 's are'} contracted (1–10 units).` : ''}`}
        />
      </div>

      <div style={{ ...cardStyle, padding: 0, overflow: 'hidden' }}>
        <div style={{ padding: '18px 24px 8px' }}>
          <SectionTitle title="Alleles" sub="Methylation is Kivvi's per-allele value: methylated reads over all reads at its CpG sites in the allele's last three repeat units." />
        </div>
        <AlleleTable alleles={fshd2.alleles} />
      </div>

      <div style={cardStyle}>
        <SectionTitle title="Methylation by allele" sub="FSHD2 is expected to lower methylation on every D4Z4 allele, whatever its haplotype or size, so the permissive allele is not the only one to look at." />
        <MethylationPlot evidence={fshd2} />
      </div>

      <div style={cardStyle}>
        <SectionTitle title="What an FSHD2 finding rests on" />
        <div style={{ display: 'flex', flexDirection: 'column', gap: 14 }}>
          <Item
            title="A permissive D4Z4 allele"
            state={`${summary.permissive_alleles} found`}
            text="Only chr4 qA alleles with an intact PolyA signal can express DUX4. They are marked in the table above."
          />
          <Item
            title="Hypomethylation across D4Z4, including alleles that are not permissive"
            state="Shown above, no cut-off applied"
            text="In the preprint's six samples with candidate SMCHD1 variants, most D4Z4 alleles were hypomethylated whatever their allele type and size, and their sample-level methylation was lower than the FSHD1 samples'. The two groups overlapped: one FSHD1 sample fell within the SMCHD1 samples' range."
          />
          <Item
            title="A pathogenic variant in SMCHD1, DNMT3B or LRIF1"
            state="Not assessed"
            text="Kivvi reports D4Z4 and nothing else, and this view does not look for variants in these genes. They have to be checked in the genome's small-variant calls."
          />
        </div>
      </div>

      <div style={{ ...cardStyle, fontSize: 12.5, color: colors.textBody, lineHeight: 1.65 }}>
        <SectionTitle title="Read with care" />
        <ul style={{ margin: 0, paddingLeft: 18, display: 'flex', flexDirection: 'column', gap: 6 }}>
          <li>
            <b>Not validated for FSHD2.</b> The preprint tested Kivvi on FSHD1 and control samples. Its six candidate SMCHD1 samples had unknown FSHD status, and it states that Kivvi's accuracy in the 11–20 unit range was not measured.
          </li>
          <li>
            <b>Size is not decisive.</b> About 10% of FSHD2 cases carry a permissive allele of 8–10 units, and a contracted permissive allele with low methylation can mean FSHD1 and FSHD2 together (the preprint flags one such candidate).
          </li>
          <li>
            <b>Sample type matters.</b> The preprint left cell lines out of its methylation analysis because methylation in lymphoblastoid cell lines can be unreliable.
          </li>
          <li>
            <b>Mosaic contractions and other D4Z4 variants</b> are not evaluated.
          </li>
        </ul>
        <div style={{ marginTop: 12, color: colors.textSecondary }}>
          Source: Chen et al., "Resolution of the D4Z4 repeat responsible for facioscapulohumeral muscular dystrophy with HiFi sequencing", bioRxiv 2026 (doi 10.64898/2026.04.10.717730), a preprint that has not been peer reviewed.
        </div>
      </div>
    </div>
  )
}

function Stat({ label, value, sub }: { label: string; value: string; sub: string }) {
  return (
    <div style={{ ...cardStyle, padding: '18px 20px' }}>
      <div style={{ fontSize: 12, fontWeight: 600, color: colors.textSecondary }}>{label}</div>
      <div style={{ fontSize: 26, fontWeight: 700, marginTop: 8 }}>{value}</div>
      <div style={{ fontSize: 12, color: colors.textTertiary, marginTop: 4, lineHeight: 1.5 }}>{sub}</div>
    </div>
  )
}

function Item({ title, state, text }: { title: string; state: string; text: ReactNode }) {
  return (
    <div style={{ display: 'grid', gridTemplateColumns: '1fr 190px', gap: 16, alignItems: 'start' }}>
      <div>
        <div style={{ fontSize: 13.5, fontWeight: 600 }}>{title}</div>
        <div style={{ fontSize: 12.5, color: colors.textBody, marginTop: 3, lineHeight: 1.55 }}>{text}</div>
      </div>
      <div style={{ fontSize: 12, fontWeight: 600, color: colors.textSecondary, textAlign: 'right' }}>{state}</div>
    </div>
  )
}

const th = { padding: '8px 14px', fontSize: 10.5, fontWeight: 600, color: colors.textSecondary, letterSpacing: '0.04em', textTransform: 'uppercase' as const, textAlign: 'left' as const, borderBottom: `1px solid ${colors.borderSofter}` }
const td = { padding: '10px 14px', fontSize: 13, borderTop: `1px solid ${colors.borderSoftest}`, verticalAlign: 'middle' as const }

function AlleleTable({ alleles }: { alleles: Fshd2Allele[] }) {
  if (!alleles.length) return <div style={{ padding: '0 24px 20px', fontSize: 13, color: colors.textSecondary }}>No D4Z4 allele was assembled.</div>
  return (
    <div style={{ overflowX: 'auto' }}>
      <table style={{ width: '100%', borderCollapse: 'collapse' }}>
        <thead>
          <tr>
            <th style={th}>Allele</th>
            <th style={th}>Units</th>
            <th style={th}>Chromosome</th>
            <th style={th}>Distal haplotype</th>
            <th style={th}>Permissive</th>
            <th style={th}>Methylation</th>
          </tr>
        </thead>
        <tbody>
          {alleles.map((a) => (
            <tr key={a.index}>
              <td style={{ ...td, fontWeight: 600 }}>{a.index}</td>
              <td style={td}>
                <span style={{ fontWeight: 700 }}>{a.size_label}</span>{' '}
                <span title={SIZE_RANGE_HELP[a.size_range]} style={{ fontSize: 11.5, color: colors.textSecondary, marginLeft: 6 }}>
                  {SIZE_RANGE_LABEL[a.size_range]}
                </span>
              </td>
              <td style={{ ...td, fontFamily: mono, fontSize: 12.5 }}>{a.upstream_group ? `${a.chromosome} · ${a.upstream_group}` : a.chromosome}</td>
              <td style={{ ...td, fontFamily: mono, fontSize: 12.5 }}>{a.distal_haplotype}</td>
              <td style={{ ...td, fontWeight: a.permissive ? 700 : 400, color: a.permissive ? PERMISSIVE : colors.textSecondary }}>{a.permissive ? 'Yes' : 'No'}</td>
              <td style={td}>
                <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
                  <div style={{ width: 110, height: 8, borderRadius: 4, background: colors.trackBg, overflow: 'hidden' }}>
                    <div style={{ width: `${(a.methylation ?? 0) * 100}%`, height: '100%', background: a.permissive ? PERMISSIVE : OTHER }} />
                  </div>
                  <span style={{ fontFamily: mono, fontSize: 12.5 }}>{fmtPct(a.methylation)}</span>
                </div>
              </td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  )
}

const W = 720
const H = 230
const PAD = { left: 44, right: 118, top: 14, bottom: 46 }

// Each allele's methylation as a point on a 0-100% axis, the permissive ones in the accent colour, with the sample's
// overall methylation as a dashed line. No threshold is drawn: none is defined.
function MethylationPlot({ evidence }: { evidence: Fshd2Evidence }) {
  const { alleles, sample_methylation } = evidence
  if (!alleles.length) return null
  const x = (i: number) => PAD.left + ((i + 0.5) / alleles.length) * (W - PAD.left - PAD.right)
  const y = (v: number) => PAD.top + (1 - v) * (H - PAD.top - PAD.bottom)
  return (
    <div>
      <svg viewBox={`0 0 ${W} ${H}`} style={{ width: '100%', maxWidth: W, height: 'auto', display: 'block' }} role="img" aria-label="Methylation of each D4Z4 allele against the sample's overall methylation">
        {[0, 0.25, 0.5, 0.75, 1].map((t) => (
          <g key={t}>
            <line x1={PAD.left} x2={W - PAD.right} y1={y(t)} y2={y(t)} stroke={colors.gridLine} />
            <text x={PAD.left - 8} y={y(t) + 4} textAnchor="end" fontSize={11} fill={colors.textTertiary}>
              {Math.round(t * 100)}%
            </text>
          </g>
        ))}
        {sample_methylation != null && (
          <g>
            <line x1={PAD.left} x2={W - PAD.right} y1={y(sample_methylation)} y2={y(sample_methylation)} stroke={colors.textQuaternary} strokeDasharray="4 3" />
            <text x={W - PAD.right + 6} y={y(sample_methylation) + 4} fontSize={10.5} fill={colors.textSecondary}>
              sample overall {fmtPct(sample_methylation)}
            </text>
          </g>
        )}
        {alleles.map((a, i) => (
          <g key={a.index}>
            {a.methylation != null && (
              <circle cx={x(i)} cy={y(a.methylation)} r={6} fill={a.permissive ? PERMISSIVE : OTHER}>
                <title>{`Allele ${a.index} · ${a.size_label} units · ${fmtPct(a.methylation)}`}</title>
              </circle>
            )}
            <text x={x(i)} y={H - 26} textAnchor="middle" fontSize={11.5} fontWeight={600} fill={colors.textPrimary}>
              {a.index}
            </text>
            <text x={x(i)} y={H - 11} textAnchor="middle" fontSize={10.5} fill={colors.textTertiary}>
              {a.size_label} units
            </text>
          </g>
        ))}
      </svg>
      <div style={{ display: 'flex', gap: 18, fontSize: 12, color: colors.textSecondary, marginTop: 6 }}>
        <Key color={PERMISSIVE} label="Permissive allele" />
        <Key color={OTHER} label="Other allele" />
      </div>
    </div>
  )
}

function Key({ color, label }: { color: string; label: string }) {
  return (
    <span style={{ display: 'inline-flex', alignItems: 'center', gap: 5 }}>
      <span style={{ width: 9, height: 9, borderRadius: '50%', background: color }} />
      {label}
    </span>
  )
}
