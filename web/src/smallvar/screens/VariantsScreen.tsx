import { useEffect, useState } from 'react'
import { DEFAULT_FILTERS, smallvarApi, type FilterState, type FreqInfo, type GenePreset, type Variant, type VariantsResponse } from '../api'
import { clinvarBg, clinvarColor, colors, fmtAf, fmtInt, mono } from '../theme'
import { Button, Pill } from '../components/ui'
import { cardStyle, pageStyle } from '../components/styles'

const PAGE_SIZE = 50
const CHROMS = [...Array.from({ length: 22 }, (_, i) => `chr${i + 1}`), 'chrX', 'chrY', 'chrM']
const KINDS: { key: string; label: string }[] = [
  { key: 'snv', label: 'SNV' },
  { key: 'ins', label: 'Insertion' },
  { key: 'del', label: 'Deletion' },
  { key: 'mnv,complex', label: 'MNV / complex' },
]
const GRID = '1.25fr 1.1fr 1.5fr 0.9fr 2.3fr 0.95fr 0.8fr 0.8fr'
const BA1 = 0.05 // ACMG/AMP stand-alone benign frequency threshold

export function VariantsScreen({ sample, initial }: { sample: string; initial: FilterState }) {
  const [form, setForm] = useState<FilterState>(initial)
  const [applied, setApplied] = useState<FilterState>(initial)
  const [page, setPage] = useState(1)
  const [resp, setResp] = useState<VariantsResponse | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [loading, setLoading] = useState(false)
  const [open, setOpen] = useState<number | null>(null)
  const [presets, setPresets] = useState<GenePreset[]>([])

  useEffect(() => {
    smallvarApi.genePresets().then((r) => setPresets(r.presets)).catch(() => undefined)
  }, [])

  useEffect(() => {
    let current = true
    setLoading(true)
    setError(null)
    smallvarApi
      .variants(sample, applied, page, PAGE_SIZE)
      .then((r) => {
        if (!current) return
        setResp(r)
        setLoading(false)
      })
      .catch((e) => {
        if (!current) return
        setError(String(e))
        setLoading(false)
      })
    return () => {
      current = false
    }
  }, [sample, applied, page])

  const set = <K extends keyof FilterState>(key: K, value: FilterState[K]) => setForm((f) => ({ ...f, [key]: value }))
  const apply = () => {
    setApplied(form)
    setPage(1)
    setOpen(null)
  }
  const reset = () => {
    setForm(DEFAULT_FILTERS)
    setApplied(DEFAULT_FILTERS)
    setPage(1)
    setOpen(null)
  }
  const resort = (sort: FilterState['sort'], order: FilterState['order']) => {
    setForm((f) => ({ ...f, sort, order }))
    setApplied((f) => ({ ...f, sort, order }))
    setPage(1)
  }
  const toggleKind = (key: string) => {
    const parts = key.split(',')
    const on = parts.every((k) => form.kinds.includes(k))
    set('kinds', on ? form.kinds.filter((k) => !parts.includes(k)) : [...form.kinds, ...parts.filter((k) => !form.kinds.includes(k))])
  }

  const total = resp?.total_matching ?? 0
  const first = total === 0 ? 0 : (page - 1) * PAGE_SIZE + 1
  const last = Math.min(page * PAGE_SIZE, total)
  const acmg = presets.find((p) => p.key === 'acmg_sf_v3_3')

  return (
    <div style={pageStyle}>
      <div style={cardStyle}>
        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fill, minmax(190px, 1fr))', gap: 14 }}>
          <Field label="ClinVar">
            <Select value={form.clinvar} onChange={(v) => set('clinvar', v as FilterState['clinvar'])} options={[
              ['any', 'Any'], ['plp', 'Pathogenic / likely pathogenic'], ['plp_conflicting', 'P/LP or conflicting'], ['has_match', 'Any ClinVar match'], ['no_match', 'No ClinVar match'],
            ]} />
          </Field>
          <Field label="ClinVar review">
            <Select value={String(form.minStars)} onChange={(v) => set('minStars', Number(v))} options={[
              ['0', 'Any'], ['1', '≥ 1 star'], ['2', '≥ 2 stars'], ['3', '≥ 3 stars (expert panel)'], ['4', '4 stars (guideline)'],
            ]} />
          </Field>
          <Field label="gnomAD frequency">
            <Select value={form.maxAf} onChange={(v) => set('maxAf', v)} options={[
              ['', 'Any'], ['0.05', '≤ 5%'], ['0.01', '≤ 1%'], ['0.001', '≤ 0.1%'], ['0.0001', '≤ 0.01%'], ['0.00001', '≤ 0.001%'],
            ]} />
          </Field>
          <Field label="Absent from gnomAD">
            <Select value={form.afMissing} onChange={(v) => set('afMissing', v as FilterState['afMissing'])} options={[
              ['include', 'Keep (treat as rare)'], ['exclude', 'Drop'],
            ]} />
          </Field>
          <Field label="HPRC frequency">
            <Select value={form.maxHprcAf} onChange={(v) => set('maxHprcAf', v)} options={[
              ['', 'Any'], ['0.1', '≤ 10%'], ['0.01', '≤ 1%'], ['0.001', '≤ 0.1%'],
            ]} />
          </Field>
          <Field label="Zygosity">
            <Select value={form.zygosity} onChange={(v) => set('zygosity', v as FilterState['zygosity'])} options={[
              ['any', 'Any'], ['het', 'Heterozygous'], ['hom', 'Homozygous'],
            ]} />
          </Field>
          <Field label="Chromosome">
            <Select value={form.chrom} onChange={(v) => set('chrom', v)} options={[['', 'All (genome-wide)'], ...CHROMS.map((c): [string, string] => [c, c])]} />
          </Field>
          <Field label="Min GQ">
            <input
              type="number" min={0} value={form.minGq} placeholder="any"
              onChange={(e) => set('minGq', e.target.value)} onKeyDown={(e) => e.key === 'Enter' && apply()}
              style={inputStyle}
            />
          </Field>
        </div>

        <div style={{ display: 'grid', gridTemplateColumns: '1fr auto', gap: 14, marginTop: 14, alignItems: 'end' }}>
          <Field label="Genes (comma or space separated; empty = every gene)">
            <input
              type="text" value={form.genes} placeholder="e.g. BRCA1, MYBPC3, F11" title={form.genes}
              onChange={(e) => set('genes', e.target.value)} onKeyDown={(e) => e.key === 'Enter' && apply()}
              style={{ ...inputStyle, fontFamily: mono }}
            />
          </Field>
          {acmg && (
            <Button onClick={() => set('genes', acmg.genes.join(', '))} title={acmg.note}>
              Fill: {acmg.label} ({acmg.genes.length})
            </Button>
          )}
        </div>

        <div style={{ display: 'flex', gap: 22, flexWrap: 'wrap', alignItems: 'center', marginTop: 14 }}>
          {KINDS.map((k) => (
            <Check key={k.key} label={k.label} checked={k.key.split(',').every((x) => form.kinds.includes(x))} onChange={() => toggleKind(k.key)} />
          ))}
          <span style={{ width: 1, height: 18, background: colors.border }} />
          <Check label="Pass the VCF FILTER only" checked={form.passOnly} onChange={(v) => set('passOnly', v)} />
          <Check label="In a gene body" checked={form.inGene} onChange={(v) => set('inGene', v)} />
          <div style={{ marginLeft: 'auto', display: 'flex', gap: 8 }}>
            <Button onClick={reset}>Reset</Button>
            <Button onClick={apply} primary>
              Apply filters
            </Button>
          </div>
        </div>
      </div>

      {resp && resp.ignored_filters.length > 0 && (
        <div style={{ ...cardStyle, padding: '12px 18px', background: colors.interpBg.intermediate, color: colors.interpColor.intermediate, fontSize: 12.5 }}>
          Not applied, because the resource is not installed for this sample: {resp.ignored_filters.map(ignoredLabel).join(', ')}. The results are not filtered on
          {resp.ignored_filters.length > 1 ? ' these' : ' this'}.
        </div>
      )}

      <div style={{ ...cardStyle, padding: 0, overflow: 'hidden', opacity: loading ? 0.6 : 1 }}>
        <div style={{ padding: '14px 20px', borderBottom: `1px solid ${colors.borderSoft}`, display: 'flex', justifyContent: 'space-between', alignItems: 'center', gap: 16, flexWrap: 'wrap' }}>
          <div style={{ fontSize: 13 }}>
            {error ? (
              <span style={{ color: colors.statusPathogenic }}>Failed to load: {error}</span>
            ) : resp ? (
              <>
                <b>{fmtInt(total)}</b> of {fmtInt(resp.total_all)} calls
                {total > 0 && <span style={{ color: colors.textSecondary }}> · showing {fmtInt(first)}–{fmtInt(last)}</span>}
              </>
            ) : (
              'Loading…'
            )}
          </div>
          <div style={{ display: 'flex', gap: 8, alignItems: 'center', fontSize: 12.5 }}>
            <span style={{ color: colors.textSecondary }}>Sort</span>
            <select value={applied.sort} onChange={(e) => resort(e.target.value as FilterState['sort'], applied.order)} style={{ ...inputStyle, width: 'auto' }}>
              <option value="position">Position</option>
              <option value="af">gnomAD frequency</option>
              <option value="stars">ClinVar review</option>
              <option value="gq">GQ</option>
            </select>
            <Button onClick={() => resort(applied.sort, applied.order === 'asc' ? 'desc' : 'asc')} title="Reverse the order">
              {applied.order === 'asc' ? '↑ ascending' : '↓ descending'}
            </Button>
          </div>
        </div>

        <div style={{ display: 'grid', gridTemplateColumns: GRID, gap: 10, padding: '9px 20px', fontSize: 10.5, fontWeight: 600, color: colors.textSecondary, letterSpacing: '0.04em', textTransform: 'uppercase', borderBottom: `1px solid ${colors.borderSofter}` }}>
          <div>Position</div>
          <div>Gene</div>
          <div>Change</div>
          <div>Genotype</div>
          <div>ClinVar</div>
          <div>gnomAD AF</div>
          <div>HPRC AF</div>
          <div>GQ · DP</div>
        </div>

        {resp?.rows.map((v) => (
          <Row key={v.id} v={v} open={open === v.id} onToggle={() => setOpen(open === v.id ? null : v.id)} gnomadInstalled={resp.gnomad_installed} />
        ))}

        {resp && total === 0 && <Empty resp={resp} />}

        {total > PAGE_SIZE && (
          <div style={{ padding: '12px 20px', display: 'flex', gap: 8, alignItems: 'center' }}>
            <Button onClick={() => setPage((p) => Math.max(1, p - 1))} disabled={page === 1}>
              ← Previous
            </Button>
            <Button onClick={() => setPage((p) => p + 1)} disabled={last >= total}>
              Next →
            </Button>
            <span style={{ fontSize: 12, color: colors.textSecondary }}>
              page {fmtInt(page)} of {fmtInt(Math.ceil(total / PAGE_SIZE))}
            </span>
          </div>
        )}
      </div>
    </div>
  )
}

function ignoredLabel(key: string): string {
  return { clinvar: 'ClinVar filter', min_stars: 'ClinVar review', max_af: 'gnomAD frequency', max_hprc_af: 'HPRC frequency' }[key] ?? key
}

const inputStyle: React.CSSProperties = {
  width: '100%', boxSizing: 'border-box', padding: '6px 8px', fontSize: 13, border: `1px solid ${colors.borderInput}`, borderRadius: 6, background: 'white', color: colors.textPrimary,
}

function Field({ label, children }: { label: string; children: React.ReactNode }) {
  return (
    <label style={{ display: 'block' }}>
      <div style={{ fontSize: 11, fontWeight: 600, color: colors.textSecondary, textTransform: 'uppercase', letterSpacing: '0.04em', marginBottom: 4 }}>{label}</div>
      {children}
    </label>
  )
}

function Select({ value, onChange, options }: { value: string; onChange: (v: string) => void; options: [string, string][] }) {
  return (
    <select value={value} onChange={(e) => onChange(e.target.value)} style={inputStyle}>
      {options.map(([v, label]) => (
        <option key={v} value={v}>
          {label}
        </option>
      ))}
    </select>
  )
}

function Check({ label, checked, onChange }: { label: string; checked: boolean; onChange: (v: boolean) => void }) {
  return (
    <label style={{ display: 'flex', alignItems: 'center', gap: 6, fontSize: 12.5, cursor: 'pointer' }}>
      <input type="checkbox" checked={checked} onChange={(e) => onChange(e.target.checked)} />
      {label}
    </label>
  )
}

const KIND_LABEL: Record<Variant['kind'], string> = { snv: 'SNV', ins: 'INS', del: 'DEL', mnv: 'MNV', complex: 'CPX' }

function trunc(s: string, n = 12): string {
  return s.length > n ? `${s.slice(0, n)}…` : s
}

function stars(n: number): string {
  return '★'.repeat(n) + '☆'.repeat(Math.max(0, 4 - n))
}

function Row({ v, open, onToggle, gnomadInstalled }: { v: Variant; open: boolean; onToggle: () => void; gnomadInstalled: boolean }) {
  const cv = v.clinvar
  const common = v.gnomad != null && v.gnomad.af > BA1
  return (
    <>
      <div
        onClick={onToggle}
        style={{ display: 'grid', gridTemplateColumns: GRID, gap: 10, padding: '8px 20px', fontSize: 12.5, alignItems: 'center', borderBottom: `1px solid ${colors.borderSoftest}`, cursor: 'pointer', background: open ? colors.chipBg : undefined }}
      >
        <div style={{ fontFamily: mono }}>
          {v.chrom}:{fmtInt(v.pos)}
        </div>
        <div style={{ fontWeight: v.genes.length ? 600 : 400, color: v.genes.length ? colors.textPrimary : colors.textTertiary, overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }} title={v.genes.join(', ')}>
          {v.genes.length ? v.genes.join(', ') : 'intergenic'}
        </div>
        <div style={{ fontFamily: mono, overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }} title={`${v.ref} → ${v.alt}`}>
          <span style={{ color: colors.textTertiary, fontSize: 10.5, marginRight: 6 }}>{KIND_LABEL[v.kind]}</span>
          {trunc(v.ref)}→{trunc(v.alt)}
        </div>
        <div>
          {v.genotype.zygosity}
          {v.genotype.phased && <span title={`phased${v.genotype.phase_set ? `, phase set ${v.genotype.phase_set}` : ''}`} style={{ color: colors.textTertiary }}> ⛓</span>}
        </div>
        <div style={{ display: 'flex', gap: 8, alignItems: 'center', minWidth: 0 }}>
          {cv ? (
            <>
              <Pill color={clinvarColor(cv.group)} bg={clinvarBg(cv.group)} title={cv.classification}>
                {cv.group === 'plp' ? shortClass(cv.classification) : 'Conflicting'}
              </Pill>
              <span style={{ color: colors.textTertiary, fontSize: 11.5, letterSpacing: 1 }} title={cv.review}>{stars(cv.stars)}</span>
              <span style={{ color: colors.textSecondary, overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }} title={cv.conditions}>{cv.conditions}</span>
            </>
          ) : (
            <span style={{ color: colors.textTertiary }}>—</span>
          )}
        </div>
        <div style={{ fontFamily: mono }}>
          <Freq f={v.gnomad} installed={gnomadInstalled} />
          {common && <span title="Above the 5% stand-alone benign (BA1) frequency threshold" style={{ marginLeft: 4, color: colors.interpColor.intermediate, fontSize: 10.5 }}>common</span>}
        </div>
        <div style={{ fontFamily: mono }}>
          <Freq f={v.hprc} installed />
        </div>
        <div style={{ fontFamily: mono, color: colors.textSecondary }}>
          {v.quality.gq ?? '—'} · {v.quality.dp ?? '—'}
        </div>
      </div>
      {open && <Detail v={v} />}
    </>
  )
}

function shortClass(c: string): string {
  if (c.startsWith('Pathogenic/Likely')) return 'P/LP'
  if (c.startsWith('Likely pathogenic')) return 'Likely pathogenic'
  return 'Pathogenic'
}

function Freq({ f, installed }: { f: FreqInfo | null; installed: boolean }) {
  if (f) return <>{fmtAf(f.af)}</>
  if (!installed) return <span style={{ color: colors.textTertiary }} title="Not evaluated: resource not installed">n/e</span>
  return <span style={{ color: colors.textTertiary }} title="Not seen at a PASS site in this resource">absent</span>
}

function Detail({ v }: { v: Variant }) {
  const cv = v.clinvar
  const q = v.quality
  const g = v.genotype
  const plain = `${v.chrom.replace('chr', '')}-${v.pos}-${v.ref}-${v.alt}`
  const showFlag = cv?.group === 'plp' && v.gnomad != null && v.gnomad.af > BA1
  return (
    <div style={{ padding: '14px 20px 18px 36px', background: colors.chipBg, borderBottom: `1px solid ${colors.borderSoft}`, display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(260px, 1fr))', gap: 22, fontSize: 12.5 }}>
      <Block title="Variant">
        <KV k="Position" v={`${v.chrom}:${fmtInt(v.pos)}`} />
        <KV k="VCF alleles" v={`${trunc(v.ref, 30)} → ${trunc(v.alt, 30)}`} />
        <KV k="Genotype" v={`${g.gt} · ${g.zygosity}`} />
        <KV k="Phasing" v={g.phased ? `phased${g.phase_set ? `, phase set ${g.phase_set}` : ''}${g.haplotypes ? `, haplotype ${g.haplotypes}` : ''}` : 'unphased'} />
        <KV k="QUAL · FILTER" v={`${q.qual ?? '—'} · ${q.filter}`} />
        <KV k="GQ · DP · VAF" v={`${q.gq ?? '—'} · ${q.dp ?? '—'} · ${q.vaf ?? '—'}`} />
      </Block>

      <Block title="ClinVar">
        {cv ? (
          <>
            <KV k="Classification" v={cv.classification} />
            <KV k="Review" v={`${stars(cv.stars)} ${cv.review}`} />
            <KV k="Name" v={cv.name} />
            <KV k="Conditions" v={cv.conditions || '—'} />
            <KV k="Last evaluated" v={`${cv.last_evaluated || '—'} · ${cv.n_submitters} submitter${cv.n_submitters === 1 ? '' : 's'}`} />
            <div style={{ marginTop: 6 }}>
              <a href={`https://www.ncbi.nlm.nih.gov/clinvar/variation/${cv.variation_id}/`} target="_blank" rel="noreferrer" style={{ color: colors.accentButtonBg }}>
                ClinVar variation {cv.variation_id} ↗
              </a>
            </div>
          </>
        ) : (
          <div style={{ color: colors.textSecondary, lineHeight: 1.5 }}>
            No match in the ClinVar subset (pathogenic, likely pathogenic and conflicting). The variant may still be in ClinVar as benign or uncertain.
          </div>
        )}
      </Block>

      <Block title="Population frequency">
        <KV k="gnomAD v3.0" v={v.gnomad ? `AF ${fmtAf(v.gnomad.af)} · AC ${fmtInt(v.gnomad.ac)} · ${fmtInt(v.gnomad.hom)} hom` : 'absent: not seen at a PASS site'} />
        <KV k="HPRC" v={v.hprc ? `AF ${fmtAf(v.hprc.af)} · AC ${fmtInt(v.hprc.ac)} · ${fmtInt(v.hprc.hom)} hom` : 'absent'} />
        {showFlag && v.gnomad && (
          <div style={{ marginTop: 8, padding: '8px 10px', borderRadius: 6, background: colors.interpBg.intermediate, color: colors.interpColor.intermediate, lineHeight: 1.5 }}>
            ClinVar lists this as {cv?.classification.toLowerCase()}, but the gnomAD frequency ({fmtAf(v.gnomad.af)}) is above the 5% stand-alone benign
            (BA1) threshold. The ClinVar assertion may be outdated or specific to a context; check it before relying on it.
          </div>
        )}
        <div style={{ marginTop: 6 }}>
          <a href={`https://gnomad.broadinstitute.org/variant/${plain}`} target="_blank" rel="noreferrer" style={{ color: colors.accentButtonBg }}>
            Open in the gnomAD browser ↗
          </a>
        </div>
      </Block>
    </div>
  )
}

function Block({ title, children }: { title: string; children: React.ReactNode }) {
  return (
    <div>
      <div style={{ fontSize: 11, fontWeight: 600, color: colors.textSecondary, textTransform: 'uppercase', letterSpacing: '0.04em', marginBottom: 6 }}>{title}</div>
      {children}
    </div>
  )
}

function KV({ k, v }: { k: string; v: string }) {
  return (
    <div style={{ display: 'flex', gap: 10, padding: '2px 0' }}>
      <div style={{ width: 104, flexShrink: 0, color: colors.textSecondary }}>{k}</div>
      <div style={{ wordBreak: 'break-word' }}>{v}</div>
    </div>
  )
}

// A filter combination that returns nothing says where the calls went.
function Empty({ resp }: { resp: VariantsResponse }) {
  return (
    <div style={{ padding: '18px 20px' }}>
      <div style={{ fontSize: 13.5, fontWeight: 600 }}>No calls match these filters</div>
      {resp.funnel && (
        <div style={{ marginTop: 10, fontSize: 12.5 }}>
          <div style={{ color: colors.textSecondary, marginBottom: 6 }}>Calls left after each filter, in order:</div>
          {resp.funnel.map((f, i) => (
            <div key={f.label} style={{ display: 'flex', gap: 14, padding: '3px 0', color: f.count === 0 && resp.funnel![i - 1]?.count !== 0 ? colors.statusPathogenic : colors.textBody, fontWeight: f.count === 0 && resp.funnel![i - 1]?.count !== 0 ? 600 : 400 }}>
              <div style={{ width: 70, textAlign: 'right', fontFamily: mono }}>{fmtInt(f.count)}</div>
              <div>{f.label}</div>
            </div>
          ))}
        </div>
      )}
    </div>
  )
}
