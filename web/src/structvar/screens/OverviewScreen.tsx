import { colors, fmtInt, fmtPct, mono, TYPE_COLOR } from '../theme'
import { structvarApi, type FilterState, type Overview, type SvType } from '../api'
import { Button, LoadState, Pill, SectionTitle, Stat } from '../components/ui'
import { cardStyle, pageStyle } from '../components/styles'
import { useLoad } from '../useLoad'

type Open = (f: Partial<FilterState>) => void

const TYPES: SvType[] = ['DEL', 'DUP', 'INS', 'INV', 'CNV', 'BND']
const SIZE_TYPES: SvType[] = ['DEL', 'DUP', 'INS', 'INV', 'CNV']

export function OverviewScreen({ sample, onOpenSvs }: { sample: string; onOpenSvs: Open }) {
  const { data, error } = useLoad(() => structvarApi.overview(sample), sample)
  if (!data) return <LoadState data={data} error={error} />

  return (
    <div style={pageStyle}>
      <CallSetCard o={data} />
      <FamilyCard o={data} />
      <FunnelCard o={data} onOpen={onOpenSvs} />
      <GenesCard o={data} onOpen={onOpenSvs} />
      <ResourcesCard o={data} />
    </div>
  )
}

function CallSetCard({ o }: { o: Overview }) {
  const q = o.qc
  const proband = q.per_sample.find((s) => s.role === 'proband')
  const maxBin = Math.max(1, ...q.size_bins.map((b) => SIZE_TYPES.reduce((n, t) => n + (b[t] ?? 0), 0)))
  return (
    <div style={cardStyle}>
      <SectionTitle
        title="Call set"
        sub={`${q.source ?? o.caller ?? 'Unknown caller'}. Counts are read from the VCF${q.samples.length > 1 ? ", which holds every family member's calls" : ''}.`}
      />
      {(q.invalid_records > 0 || q.skipped > 0) && (
        <div style={{ background: colors.interpBg.intermediate, color: colors.interpColor.intermediate, borderRadius: 8, padding: '8px 12px', fontSize: 12.5, marginBottom: 14, lineHeight: 1.5 }}>
          {q.invalid_records > 0 && (
            <div>
              {fmtInt(q.invalid_records)} record{q.invalid_records === 1 ? '' : 's'} with an invalid allele {q.invalid_records === 1 ? 'was' : 'were'} skipped (a tool that stops at the
              first malformed record, such as bcftools, would have dropped every record after it).
            </div>
          )}
          {q.skipped > 0 && <div>{fmtInt(q.skipped)} records of a type this suite does not handle were skipped.</div>}
        </div>
      )}
      <div style={{ display: 'flex', gap: 36, flexWrap: 'wrap', marginBottom: 20 }}>
        <Stat label="Records" value={fmtInt(q.records)} sub="every call in the VCF" />
        <Stat label="PASS, not breakend" value={fmtInt(q.pass_non_bnd)} sub="passes the caller's FILTER" />
        {proband && <Stat label="Carried by proband" value={fmtInt(proband.carrier)} sub={`of those, ${proband.sample}`} />}
        <Stat label="Breakends" value={fmtInt(q.by_type_all.BND)} sub="most are artifacts" />
        <Stat label="Imprecise" value={fmtInt(q.imprecise)} sub="breakpoints not resolved to the base" />
        {q.de_novo_candidates != null && <Stat label="De novo candidates" value={fmtInt(q.de_novo_candidates)} sub="mostly parental missed calls" />}
      </div>

      <div style={{ display: 'grid', gridTemplateColumns: 'minmax(260px, 340px) 1fr', gap: 36, alignItems: 'start' }}>
        <div>
          <SubHead>By type</SubHead>
          <div style={{ display: 'grid', gridTemplateColumns: '60px 1fr 1fr', gap: 8, fontSize: 10.5, fontWeight: 600, color: colors.textSecondary, textTransform: 'uppercase', letterSpacing: '0.04em', padding: '4px 0', borderBottom: `1px solid ${colors.borderSofter}` }}>
            <div />
            <div style={{ textAlign: 'right' }}>All</div>
            <div style={{ textAlign: 'right' }}>PASS</div>
          </div>
          {TYPES.map((t) => (
            <div key={t} style={{ display: 'grid', gridTemplateColumns: '60px 1fr 1fr', gap: 8, fontSize: 12.5, padding: '5px 0', borderBottom: `1px solid ${colors.borderSoftest}` }}>
              <div style={{ fontWeight: 600, color: TYPE_COLOR[t] }}>{t}</div>
              <div style={{ textAlign: 'right', fontFamily: mono }}>{fmtInt(q.by_type_all[t])}</div>
              <div style={{ textAlign: 'right', fontFamily: mono }}>{fmtInt(q.by_type_pass[t])}</div>
            </div>
          ))}
          <SubHead style={{ marginTop: 18 }}>FILTER values</SubHead>
          {q.filters.slice(0, 8).map((f) => (
            <div key={f.filter} style={{ display: 'flex', justifyContent: 'space-between', fontSize: 12.5, padding: '4px 0', borderBottom: `1px solid ${colors.borderSoftest}` }}>
              <span style={{ fontFamily: mono }}>
                {f.filter === '.' ? '. (not set, counted as PASS)' : f.filter}
              </span>
              <span style={{ fontFamily: mono, color: colors.textSecondary }}>{fmtInt(f.count)}</span>
            </div>
          ))}
        </div>

        <div>
          <SubHead>Size of PASS calls (excluding breakends)</SubHead>
          <div style={{ display: 'flex', alignItems: 'flex-end', gap: 14, height: 150, marginTop: 8 }}>
            {q.size_bins.map((b) => {
              const total = SIZE_TYPES.reduce((n, t) => n + (b[t] ?? 0), 0)
              return (
                <div key={b.bin} style={{ flex: 1, maxWidth: 90, textAlign: 'center' }} title={`${b.bin}: ${fmtInt(total)} calls`}>
                  <div style={{ fontSize: 11, fontFamily: mono, color: colors.textSecondary, marginBottom: 3 }}>{fmtInt(total)}</div>
                  <div style={{ height: Math.max(2, (total / maxBin) * 100), display: 'flex', flexDirection: 'column-reverse', borderRadius: '3px 3px 0 0', overflow: 'hidden' }}>
                    {SIZE_TYPES.map((t) => {
                      const n = b[t] ?? 0
                      return n ? <div key={t} title={`${t} ${fmtInt(n)}`} style={{ flexGrow: n, background: TYPE_COLOR[t], opacity: 0.9 }} /> : null
                    })}
                  </div>
                  <div style={{ fontSize: 10.5, color: colors.textSecondary, marginTop: 4 }}>{b.bin}</div>
                </div>
              )
            })}
          </div>
          <div style={{ display: 'flex', gap: 14, marginTop: 12, fontSize: 12, color: colors.textSecondary }}>
            {SIZE_TYPES.map((t) => (
              <span key={t} style={{ display: 'flex', alignItems: 'center', gap: 5 }}>
                <span style={{ width: 10, height: 10, borderRadius: 2, background: TYPE_COLOR[t] }} />
                {t}
              </span>
            ))}
          </div>
        </div>
      </div>
    </div>
  )
}

function SubHead({ children, style }: { children: React.ReactNode; style?: React.CSSProperties }) {
  return <div style={{ fontSize: 11, fontWeight: 600, color: colors.textSecondary, textTransform: 'uppercase', letterSpacing: '0.04em', ...style }}>{children}</div>
}

function FamilyCard({ o }: { o: Overview }) {
  const q = o.qc
  const r = q.roles
  return (
    <div style={cardStyle}>
      <SectionTitle title="Samples and family structure" sub="Who is the proband and who are the parents, and what that is based on." />
      <div style={{ display: 'flex', gap: 10, alignItems: 'baseline', flexWrap: 'wrap', marginBottom: 12 }}>
        {r.source === 'inferred' && <Pill color={colors.interpColor.intermediate} bg={colors.interpBg.intermediate}>inferred from the genotypes</Pill>}
        {r.source === 'sample.json' && <Pill color={colors.interpColor.normal} bg={colors.interpBg.normal}>stated in sample.json</Pill>}
        {r.source === 'none' && <Pill color={colors.textSecondary} bg={colors.chipBg}>no family structure assumed</Pill>}
        <span style={{ fontSize: 12.5, color: colors.textBody, lineHeight: 1.55 }}>{r.basis}</span>
      </div>
      <div style={{ display: 'grid', gridTemplateColumns: '1.3fr 1fr 1fr 1fr 1fr', gap: 10, padding: '6px 0', fontSize: 10.5, fontWeight: 600, color: colors.textSecondary, textTransform: 'uppercase', letterSpacing: '0.04em', borderBottom: `1px solid ${colors.borderSofter}` }}>
        <div>Sample</div>
        <div>Role</div>
        <div style={{ textAlign: 'right' }}>Carries</div>
        <div style={{ textAlign: 'right' }}>Reference</div>
        <div style={{ textAlign: 'right' }}>No call</div>
      </div>
      {q.per_sample.map((s) => {
        const total = s.carrier + s.ref + s.nocall
        return (
          <div key={s.sample} style={{ display: 'grid', gridTemplateColumns: '1.3fr 1fr 1fr 1fr 1fr', gap: 10, padding: '7px 0', fontSize: 12.5, alignItems: 'center', borderBottom: `1px solid ${colors.borderSoftest}` }}>
            <div style={{ fontFamily: mono, fontWeight: 600 }}>{s.sample}</div>
            <div>{s.role}</div>
            <div style={{ textAlign: 'right', fontFamily: mono }}>{fmtInt(s.carrier)}</div>
            <div style={{ textAlign: 'right', fontFamily: mono }}>{fmtInt(s.ref)}</div>
            <div style={{ textAlign: 'right', fontFamily: mono, color: s.nocall / Math.max(total, 1) > 0.05 ? colors.interpColor.intermediate : colors.textSecondary }}>
              {fmtInt(s.nocall)} <span style={{ color: colors.textTertiary }}>({fmtPct(s.nocall / Math.max(total, 1), 1)})</span>
            </div>
          </div>
        )
      })}
      <div style={{ fontSize: 11.5, color: colors.textTertiary, marginTop: 8 }}>Counts are over PASS calls that are not breakends.</div>
    </div>
  )
}

// Each funnel step is a real filter on the SV screen, so the count can be reproduced there. A step's
// filters are those of every step up to and including it.
function presetFor(i: number, steps: { label: string }[]): Partial<FilterState> | null {
  const at = (prefix: string) => steps.findIndex((s) => s.label.startsWith(prefix))
  const reached = (prefix: string) => at(prefix) >= 0 && i >= at(prefix)
  if (i === 0) return null
  const f: Partial<FilterState> = { carried: reached('Carried'), includeBnd: !reached('Not a breakend') }
  if (reached('Rare')) f.maxPopAf = '0.01'
  if (reached('Hits a gene')) f.impact = 'exon'
  if (reached('Outside segmental')) {
    f.maxSegdup = '50'
    f.noRepeat = true
  }
  if (reached('De novo')) f.inheritance = 'de_novo'
  return f
}

function FunnelCard({ o, onOpen }: { o: Overview; onOpen: Open }) {
  const top = o.funnel[0]?.count || 1
  const hasRare = o.funnel.some((f) => f.label.startsWith('Rare'))
  const hasGene = o.funnel.some((f) => f.label.startsWith('Hits a gene'))
  return (
    <div style={cardStyle}>
      <SectionTitle title="From call set to review list" sub="Each step adds one filter to the ones before. Open a step to see exactly those calls." />
      {o.funnel.map((step, i) => {
        const preset = presetFor(i, o.funnel)
        const prev = i > 0 ? o.funnel[i - 1].count : null
        return (
          <div key={step.label} style={{ display: 'grid', gridTemplateColumns: '290px 1fr 150px 90px', gap: 14, alignItems: 'center', padding: '7px 0', borderBottom: `1px solid ${colors.borderSoftest}` }}>
            <div style={{ fontSize: 13, fontWeight: 500 }}>{step.label}</div>
            <div style={{ background: colors.trackBg, borderRadius: 4, height: 14, overflow: 'hidden' }}>
              <div style={{ width: `${Math.max(0.6, (step.count / top) * 100)}%`, height: '100%', background: i === o.funnel.length - 1 ? colors.statusUncertain : colors.cluster1, opacity: 0.85 }} />
            </div>
            <div style={{ fontFamily: mono, fontSize: 12.5, textAlign: 'right' }}>
              {fmtInt(step.count)}
              {prev ? <span style={{ color: colors.textTertiary }}> · {((step.count / Math.max(prev, 1)) * 100).toFixed(step.count / Math.max(prev, 1) < 0.1 ? 1 : 0)}%</span> : null}
            </div>
            <div style={{ textAlign: 'right' }}>{preset && <Button onClick={() => onOpen(preset)}>View</Button>}</div>
          </div>
        )
      })}
      {!hasRare && <Note>No population catalog is installed, so there is no &ldquo;rare&rdquo; step.</Note>}
      {!hasGene && <Note>No gene models are installed, so there is no gene step.</Note>}
      {o.trio && (
        <Note>
          &ldquo;De novo&rdquo; here means the proband carries a call that neither parent does. Most such calls are a parent&rsquo;s missed call, so confirm them: look at the
          parents&rsquo; read depth, the proband&rsquo;s supporting reads, and whether the call sits in a repeat or segmental duplication.
        </Note>
      )}
    </div>
  )
}

function Note({ children }: { children: React.ReactNode }) {
  return <div style={{ fontSize: 12, color: colors.interpColor.intermediate, marginTop: 10, lineHeight: 1.55 }}>{children}</div>
}

function GenesCard({ o, onOpen }: { o: Overview; onOpen: Open }) {
  if (!o.top_genes.length) return null
  return (
    <div style={cardStyle}>
      <SectionTitle title="Protein-coding genes hit by rare calls" sub="Rare calls the proband carries that hit an exon of, or delete or duplicate, a protein-coding gene. The number is how many calls. Pseudogenes and non-coding RNAs are in the table but not listed here." />
      <div style={{ display: 'flex', gap: 8, flexWrap: 'wrap' }}>
        {o.top_genes.map((g) => (
          <span
            key={g.gene}
            onClick={() => onOpen({ maxPopAf: '0.01', impact: 'exon', genes: g.gene })}
            style={{ cursor: 'pointer', fontSize: 12.5, fontFamily: mono, padding: '4px 10px', borderRadius: 999, background: colors.chipBg, border: `1px solid ${colors.borderSoft}` }}
          >
            {g.gene}
            {g.count > 1 ? ` ×${g.count}` : ''}
          </span>
        ))}
      </div>
    </div>
  )
}

function ResourcesCard({ o }: { o: Overview }) {
  return (
    <div style={cardStyle}>
      <SectionTitle
        title="What this tool can and cannot say"
        sub="Gene impact is overlap with the MANE transcript, not a predicted consequence. Population matches are by breakpoints and size, not sequence. An evidence source marked not evaluated has no data here."
      />
      {o.resources.map((r) => (
        <div key={r.key} style={{ display: 'grid', gridTemplateColumns: '240px 90px 1fr', gap: 10, padding: '7px 0', fontSize: 12.5, alignItems: 'start', borderBottom: `1px solid ${colors.borderSoftest}` }}>
          <div style={{ fontWeight: 600 }}>{r.label}</div>
          <div>
            {r.installed ? <Pill color={colors.interpColor.normal} bg={colors.interpBg.normal}>installed</Pill> : <Pill color={colors.textSecondary} bg={colors.chipBg}>missing</Pill>}
          </div>
          <div style={{ color: colors.textSecondary, lineHeight: 1.5 }}>{r.detail}</div>
        </div>
      ))}
    </div>
  )
}
