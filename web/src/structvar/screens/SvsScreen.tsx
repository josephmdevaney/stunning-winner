import { useEffect, useState } from 'react'
import { DEFAULT_FILTERS, structvarApi, type FamilyMember, type FilterState, type GeneHit, type GeneInfo, type Sv, type SvsResponse } from '../api'
import { colors, fmtAf, fmtInt, fmtSize, IMPACT_LABEL, mono, TYPE_COLOR } from '../theme'
import { Button, Pill } from '../components/ui'
import { cardStyle, pageStyle } from '../components/styles'
import { RegionView } from '../components/RegionView'

const PAGE_SIZE = 50
const CHROMS = [...Array.from({ length: 22 }, (_, i) => `chr${i + 1}`), 'chrX', 'chrY', 'chrM']
const TYPES = ['DEL', 'DUP', 'INS', 'INV', 'CNV']
const GRID = '1.2fr 0.9fr 1.7fr 1.05fr 1.05fr 1.05fr 1.35fr 0.55fr'
const SOURCE_LABEL: Record<string, string> = {
  gnomad_sv: 'gnomAD-SV (short-read)',
  hprc: 'HPRC (long-read, pbsv)',
  decode: 'DECODE (long-read, ONT)',
  eee: 'EEE (long-read discovery)',
}

export function SvsScreen({ sample, initial }: { sample: string; initial: FilterState }) {
  const [form, setForm] = useState<FilterState>(initial)
  const [applied, setApplied] = useState<FilterState>(initial)
  const [page, setPage] = useState(1)
  const [resp, setResp] = useState<SvsResponse | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [loading, setLoading] = useState(false)
  const [open, setOpen] = useState<number | null>(null)

  useEffect(() => {
    let current = true
    setLoading(true)
    setError(null)
    structvarApi
      .svs(sample, applied, page, PAGE_SIZE)
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
  const toggleType = (t: string) => set('types', form.types.includes(t) ? form.types.filter((x) => x !== t) : [...form.types, t])

  const total = resp?.total_matching ?? 0
  const first = total === 0 ? 0 : (page - 1) * PAGE_SIZE + 1
  const last = Math.min(page * PAGE_SIZE, total)
  const trio = resp?.trio ?? true

  return (
    <div style={pageStyle}>
      <div style={cardStyle}>
        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fill, minmax(180px, 1fr))', gap: 14 }}>
          <Field label="Size (bp)">
            <div style={{ display: 'flex', gap: 6 }}>
              <input type="number" min={0} value={form.minSize} placeholder="min" onChange={(e) => set('minSize', e.target.value)} onKeyDown={(e) => e.key === 'Enter' && apply()} style={inputStyle} />
              <input type="number" min={0} value={form.maxSize} placeholder="max" onChange={(e) => set('maxSize', e.target.value)} onKeyDown={(e) => e.key === 'Enter' && apply()} style={inputStyle} />
            </div>
          </Field>
          <Field label="Genotype">
            <Select value={form.zygosity} onChange={(v) => set('zygosity', v as FilterState['zygosity'])} options={[['any', 'Any'], ['het', 'Heterozygous'], ['hom', 'Homozygous or hemizygous']]} />
          </Field>
          <Field label="Gene impact">
            <Select value={form.impact} onChange={(v) => set('impact', v as FilterState['impact'])} options={[
              ['any', 'Any'], ['gene', 'In a gene body'], ['exon', 'Exon or whole gene'], ['coding', 'Coding exon or whole gene'], ['whole_gene', 'Whole gene'],
            ]} />
          </Field>
          <Field label="Population frequency">
            <Select value={form.maxPopAf} onChange={(v) => set('maxPopAf', v)} options={[
              ['', 'Any'], ['0.05', '≤ 5%'], ['0.01', '≤ 1%'], ['0.001', '≤ 0.1%'], ['0.0001', '≤ 0.01%'],
            ]} />
          </Field>
          <Field label="Inheritance (trio)">
            <Select value={trio ? form.inheritance : 'any'} onChange={(v) => set('inheritance', v as FilterState['inheritance'])} disabled={!trio} options={[
              ['any', trio ? 'Any' : 'No parents'], ['de_novo', 'De novo candidate'], ['inherited', 'Inherited'], ['unknown', 'Unknown (parent not called)'],
            ]} />
          </Field>
          <Field label="Parents' depth ≥">
            <input type="number" min={0} value={form.minParentDp} placeholder="any" disabled={!trio} onChange={(e) => set('minParentDp', e.target.value)} onKeyDown={(e) => e.key === 'Enter' && apply()} style={inputStyle} />
          </Field>
          <Field label={`Supporting ${resp?.support_unit ?? 'reads'} ≥`}>
            <input type="number" min={0} value={form.minAlt} placeholder="any" onChange={(e) => set('minAlt', e.target.value)} onKeyDown={(e) => e.key === 'Enter' && apply()} style={inputStyle} />
          </Field>
          <Field label="Min GQ">
            <input type="number" min={0} value={form.minGq} placeholder="any" onChange={(e) => set('minGq', e.target.value)} onKeyDown={(e) => e.key === 'Enter' && apply()} style={inputStyle} />
          </Field>
          <Field label="Segmental dup. below (%)">
            <Select value={form.maxSegdup} onChange={(v) => set('maxSegdup', v)} options={[['', 'Any'], ['10', '10%'], ['50', '50%']]} />
          </Field>
          <Field label="Chromosome">
            <Select value={form.chrom} onChange={(v) => set('chrom', v)} options={[['', 'All (genome-wide)'], ...CHROMS.map((c): [string, string] => [c, c])]} />
          </Field>
        </div>

        <div style={{ marginTop: 14 }}>
          <Field label="Genes (comma or space separated; empty = every gene)">
            <input type="text" value={form.genes} placeholder="e.g. DMD, SMN1, PMP22" onChange={(e) => set('genes', e.target.value)} onKeyDown={(e) => e.key === 'Enter' && apply()} style={{ ...inputStyle, fontFamily: mono }} />
          </Field>
        </div>

        <div style={{ display: 'flex', gap: 20, flexWrap: 'wrap', alignItems: 'center', marginTop: 14 }}>
          {TYPES.map((t) => (
            <Check key={t} label={t} checked={form.types.includes(t)} onChange={() => toggleType(t)} color={TYPE_COLOR[t]} />
          ))}
          <span style={{ width: 1, height: 18, background: colors.border }} />
          <Check label="Carried by the proband" checked={form.carried} onChange={(v) => set('carried', v)} />
          <Check label="Pass the caller's FILTER" checked={form.passOnly} onChange={(v) => set('passOnly', v)} />
          <Check label="Include breakends" checked={form.includeBnd} onChange={(v) => set('includeBnd', v)} />
          <Check label="Precise only" checked={form.preciseOnly} onChange={(v) => set('preciseOnly', v)} />
        </div>
        <div style={{ display: 'flex', gap: 20, flexWrap: 'wrap', alignItems: 'center', marginTop: 10 }}>
          <Check label="In no population catalog" checked={form.novel} onChange={(v) => set('novel', v)} />
          <Check
            label="Not mostly repeat"
            checked={form.noRepeat}
            onChange={(v) => set('noRepeat', v)}
            hint="Most SVs sit in annotated repeats (on the GIAB HG002 truth set about 94% of deletions and 90% of insertions), so this removes most calls. It suits screening for false de novo calls, not browsing."
          />
          <Check label="Not in a gap, telomere, PAR or V(D)J region" checked={form.noOdd} onChange={(v) => set('noOdd', v)} />
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
          Not applied, because the data behind {resp.ignored_filters.length > 1 ? 'them is' : 'it is'} not available for this call set: {resp.ignored_filters.map(ignoredLabel).join(', ')}. The results
          are not filtered on {resp.ignored_filters.length > 1 ? 'these' : 'this'}.
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
              <option value="size">Size</option>
              <option value="impact">Gene impact</option>
              <option value="af">Population frequency</option>
              <option value="gq">GQ</option>
            </select>
            <Button onClick={() => resort(applied.sort, applied.order === 'asc' ? 'desc' : 'asc')} title="Reverse the order">
              {applied.order === 'asc' ? '↑ ascending' : '↓ descending'}
            </Button>
          </div>
        </div>

        <div style={{ display: 'grid', gridTemplateColumns: GRID, gap: 10, padding: '9px 20px', fontSize: 10.5, fontWeight: 600, color: colors.textSecondary, letterSpacing: '0.04em', textTransform: 'uppercase', borderBottom: `1px solid ${colors.borderSofter}` }}>
          <div>Position</div>
          <div>Type · size</div>
          <div>Genes</div>
          <div>Proband</div>
          <div>Inheritance</div>
          <div>Population</div>
          <div>Flags</div>
          <div>GQ</div>
        </div>

        {resp?.rows.map((sv) => (
          <Row key={sv.id} sv={sv} open={open === sv.id} onToggle={() => setOpen(open === sv.id ? null : sv.id)} sample={sample} popInstalled={resp.population_installed} unit={resp.support_unit} />
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
  return { max_pop_af: 'population frequency', novel: 'in no population catalog', impact: 'gene impact', genes: 'gene list', tracks: 'segmental duplication / repeat / region filters', inheritance: 'inheritance', min_parent_dp: "parents' depth" }[key] ?? key
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

function Select({ value, onChange, options, disabled = false }: { value: string; onChange: (v: string) => void; options: [string, string][]; disabled?: boolean }) {
  return (
    <select value={value} onChange={(e) => onChange(e.target.value)} disabled={disabled} style={{ ...inputStyle, opacity: disabled ? 0.55 : 1 }}>
      {options.map(([v, label]) => (
        <option key={v} value={v}>
          {label}
        </option>
      ))}
    </select>
  )
}

function Check({ label, checked, onChange, color, hint }: { label: string; checked: boolean; onChange: (v: boolean) => void; color?: string; hint?: string }) {
  return (
    <label title={hint} style={{ display: 'flex', alignItems: 'center', gap: 6, fontSize: 12.5, cursor: 'pointer', color: color ?? colors.textPrimary, fontWeight: color ? 600 : 400 }}>
      <input type="checkbox" checked={checked} onChange={(e) => onChange(e.target.checked)} />
      {label}
    </label>
  )
}

function exonText(h: GeneHit): string {
  if (h.relation === 'whole_gene') return 'whole gene'
  if (!h.exons.length) return ''
  const runs: string[] = []
  let a = h.exons[0]
  let b = a
  for (const n of h.exons.slice(1)) {
    if (n === b + 1) b = n
    else {
      runs.push(a === b ? `${a}` : `${a}–${b}`)
      a = b = n
    }
  }
  runs.push(a === b ? `${a}` : `${a}–${b}`)
  return `exon${h.exons.length > 1 ? 's' : ''} ${runs.join(', ')} of ${h.n_exons}`
}

const INH_LABEL: Record<string, string> = { de_novo: 'de novo candidate', inherited: 'inherited', unknown: 'unknown' }

function Row({ sv, open, onToggle, sample, popInstalled, unit }: { sv: Sv; open: boolean; onToggle: () => void; sample: string; popInstalled: boolean; unit: string }) {
  const top = sv.genes[0]
  const zyg = sv.focus.copies == null ? '—' : sv.focus.copies >= 2 ? 'hom' : sv.focus.copies === 1 ? (sv.focus.ploidy === 1 ? 'hemizygous' : 'het') : sv.focus.copies === 0 ? 'ref' : 'no call'
  return (
    <>
      <div
        onClick={onToggle}
        style={{ display: 'grid', gridTemplateColumns: GRID, gap: 10, padding: '8px 20px', fontSize: 12.5, alignItems: 'center', borderBottom: `1px solid ${colors.borderSoftest}`, cursor: 'pointer', background: open ? colors.chipBg : undefined }}
      >
        <div style={{ fontFamily: mono }}>
          {sv.chrom}:{fmtInt(sv.pos)}
        </div>
        <div>
          <span style={{ fontFamily: mono, fontWeight: 700, color: TYPE_COLOR[sv.svtype] }}>{sv.svtype}</span>{' '}
          <span style={{ fontFamily: mono, color: colors.textSecondary }}>{sv.svtype === 'BND' ? (sv.mate ? `→ ${sv.mate.chrom}` : '') : fmtSize(sv.size)}</span>
        </div>
        <div style={{ minWidth: 0, overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }} title={sv.genes.map((g) => `${g.symbol} (${IMPACT_LABEL[g.relation]})`).join(', ')}>
          {top ? (
            <>
              <b>{top.symbol}</b> <span style={{ color: colors.textSecondary }}>{exonText(top) || IMPACT_LABEL[top.relation].toLowerCase()}</span>
              {sv.n_genes > 1 && <span style={{ color: colors.textTertiary }}> +{sv.n_genes - 1}</span>}
            </>
          ) : (
            <span style={{ color: colors.textTertiary }}>intergenic</span>
          )}
        </div>
        <div style={{ fontFamily: mono }}>
          {zyg}
          {sv.focus.alt != null && <span style={{ color: colors.textTertiary }}> · {sv.focus.alt} {unit}</span>}
          {sv.focus.cn != null && <span style={{ color: colors.textTertiary }}> · CN {sv.focus.cn}</span>}
        </div>
        <div>
          {sv.inheritance ? (
            <Pill
              color={sv.inheritance === 'de_novo' ? colors.interpColor.uncertain : colors.textSecondary}
              bg={sv.inheritance === 'de_novo' ? colors.interpBg.uncertain : colors.chipBg}
              title={sv.inheritance === 'de_novo' ? 'No call in either parent. Often a parent’s missed call; check the parents’ depth.' : undefined}
            >
              {INH_LABEL[sv.inheritance]}
            </Pill>
          ) : (
            <span style={{ color: colors.textTertiary }}>—</span>
          )}
        </div>
        <div style={{ fontFamily: mono }}>
          {!popInstalled ? (
            <span style={{ color: colors.textTertiary }} title="Not evaluated: no population catalog installed">n/e</span>
          ) : sv.pop_n === 0 ? (
            <span style={{ color: colors.interpColor.uncertain }} title="Matches no call in any population catalog">novel</span>
          ) : (
            <>
              {sv.pop_af_max != null ? fmtAf(sv.pop_af_max) : 'seen'}
              <span style={{ color: colors.textTertiary }}> · {sv.pop_n}/4</span>
            </>
          )}
        </div>
        <div style={{ display: 'flex', gap: 5, flexWrap: 'wrap' }}>
          {sv.imprecise && <Flag>imprecise</Flag>}
          {sv.segdup_frac != null && sv.segdup_frac >= 0.1 && <Flag>segdup {Math.round(sv.segdup_frac * 100)}%</Flag>}
          {sv.repeat && sv.repeat_frac != null && sv.repeat_frac >= 0.5 && <Flag>{sv.repeat}</Flag>}
          {sv.odd.length > 0 && <Flag>{sv.odd.join(', ')}</Flag>}
          {!sv.passed && <Flag warn>{sv.filter}</Flag>}
        </div>
        <div style={{ fontFamily: mono, color: colors.textSecondary }}>{sv.focus.gq ?? '—'}</div>
      </div>
      {open && <Detail sv={sv} sample={sample} unit={unit} />}
    </>
  )
}

function Flag({ children, warn = false }: { children: React.ReactNode; warn?: boolean }) {
  return (
    <span style={{ fontSize: 10.5, padding: '1px 6px', borderRadius: 4, background: warn ? colors.interpBg.intermediate : colors.chipBg, color: warn ? colors.interpColor.intermediate : colors.textSecondary, border: `1px solid ${colors.borderSoft}`, whiteSpace: 'nowrap' }}>
      {children}
    </span>
  )
}

function Detail({ sv, sample, unit }: { sv: Sv; sample: string; unit: string }) {
  const [gene, setGene] = useState<string | null>(sv.genes[0]?.symbol ?? null)
  const loc = sv.svtype === 'BND' && sv.mate ? `${sv.chrom}:${fmtInt(sv.pos)} ↔ ${sv.mate.chrom}:${fmtInt(sv.mate.pos)}` : sv.svtype === 'INS' ? `${sv.chrom}:${fmtInt(sv.pos)}` : `${sv.chrom}:${fmtInt(sv.pos)}–${fmtInt(sv.end)}`
  return (
    <div style={{ padding: '14px 20px 20px 36px', background: colors.chipBg, borderBottom: `1px solid ${colors.borderSoft}`, fontSize: 12.5 }}>
      <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(300px, 1fr))', gap: 24 }}>
        <Block title="Call">
          <KV k="Location" v={loc} />
          <KV k="Type · size" v={`${sv.svtype} · ${fmtSize(sv.size)}${sv.imprecise ? ' · imprecise breakpoints' : ''}`} />
          <KV k="FILTER · QUAL" v={`${sv.filter} · ${sv.qual ?? '—'}`} />
          <KV k="Caller ID" v={sv.vid + (sv.event ? ` (event ${sv.event})` : '')} />
          {sv.svtype === 'BND' && <KV k="Breakend" v="One end of a rearrangement; it is not a call of its own and many are artifacts." />}
        </Block>
        <Block title="Family genotypes">
          <FamilyTable family={sv.family} unit={unit} />
          {sv.inheritance === 'de_novo' && (
            <div style={{ marginTop: 8, padding: '8px 10px', borderRadius: 6, background: colors.interpBg.uncertain, color: colors.interpColor.uncertain, lineHeight: 1.5 }}>
              Neither parent has a call here
              {sv.parent_min_dp != null ? `; the thinner parent was read ${sv.parent_min_dp} times` : ''}. A de novo SV is rare, and most calls like this are a parent&rsquo;s missed call.
              Confirm with the parents&rsquo; reads.
            </div>
          )}
        </Block>
        <Block title="Population catalogs">
          {Object.keys(SOURCE_LABEL).map((s) => {
            const m = sv.population[s]
            return (
              <div key={s} style={{ padding: '3px 0' }}>
                <div style={{ color: colors.textSecondary }}>{SOURCE_LABEL[s]}</div>
                <div style={{ fontFamily: mono }}>
                  {m
                    ? `${m.svtype} match (${Math.round(m.similarity * 100)}% ${sv.svtype === 'INS' ? 'size ratio' : 'overlap'}${m.n > 1 ? `, ${m.n} records` : ''})${m.af != null ? ` · AF ${fmtAf(m.af)}${m.popmax != null ? ` (max ancestry ${fmtAf(m.popmax)})` : ''}` : m.ac != null ? ` · in ${m.ac} genomes` : ' · presence only'}`
                    : 'no match'}
                </div>
              </div>
            )
          })}
          <div style={{ color: colors.textTertiary, marginTop: 6, lineHeight: 1.5 }}>Matched by breakpoints and size, not by sequence.</div>
        </Block>
        <Block title="Genome context">
          <KV k="Segmental dup." v={sv.segdup_frac == null ? 'not evaluated' : `${Math.round(sv.segdup_frac * 100)}% of the call`} />
          <KV k="Repeat" v={sv.repeat ? `${sv.repeat} (${Math.round((sv.repeat_frac ?? 0) * 100)}%)` : 'none found (or not evaluated for calls over 100 kb)'} />
          <KV k="Special region" v={sv.odd.length ? sv.odd.join(', ') : 'none'} />
        </Block>
      </div>

      <div style={{ marginTop: 20 }}>
        <Block title="Where it falls">
          <RegionView sample={sample} sv={sv} />
        </Block>
      </div>

      {sv.genes.length > 0 && (
        <div style={{ marginTop: 20, display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(380px, 1fr))', gap: 24 }}>
          <Block title={`Genes (${sv.n_genes})`}>
            <div style={{ maxHeight: 220, overflowY: 'auto' }}>
              {sv.genes.map((g) => (
                <div key={g.gene_id} onClick={() => setGene(g.symbol)} style={{ display: 'flex', gap: 10, alignItems: 'baseline', padding: '3px 0', cursor: 'pointer', fontWeight: gene === g.symbol ? 700 : 400 }}>
                  <span style={{ fontFamily: mono, width: 90, flexShrink: 0 }}>{g.symbol}</span>
                  <span style={{ color: colors.textSecondary }}>
                    {IMPACT_LABEL[g.relation]}
                    {exonText(g) && g.relation !== 'whole_gene' ? ` · ${exonText(g)}` : ''} · {g.transcript}
                    {g.model_source === 'MANE Select' ? ' (MANE Select)' : ''}
                  </span>
                </div>
              ))}
            </div>
            {sv.n_genes > sv.genes.length && <div style={{ color: colors.textTertiary, marginTop: 4 }}>Showing the {sv.genes.length} most severe of {sv.n_genes}.</div>}
            <div style={{ color: colors.textTertiary, marginTop: 6, lineHeight: 1.5 }}>
              Overlap with the transcript&rsquo;s exons, not a predicted consequence: whether the reading frame is kept is not worked out.
            </div>
          </Block>
          {gene && <GeneBlock symbol={gene} />}
        </div>
      )}
    </div>
  )
}

function FamilyTable({ family, unit }: { family: FamilyMember[]; unit: string }) {
  const cols = '1.2fr 0.8fr 0.7fr 1fr 0.5fr'
  return (
    <div>
      <div style={{ display: 'grid', gridTemplateColumns: cols, gap: 8, fontSize: 10.5, fontWeight: 600, color: colors.textSecondary, textTransform: 'uppercase', letterSpacing: '0.04em', paddingBottom: 3 }}>
        <div>Sample</div>
        <div>Role</div>
        <div>GT</div>
        <div style={{ textTransform: 'none' }}>{unit.charAt(0).toUpperCase() + unit.slice(1)} (alt / ref)</div>
        <div>GQ</div>
      </div>
      {family.map((m) => (
        <div key={m.sample} style={{ display: 'grid', gridTemplateColumns: cols, gap: 8, fontFamily: mono, padding: '2px 0', fontWeight: m.role === 'proband' ? 700 : 400 }}>
          <div>{m.sample}</div>
          <div style={{ fontFamily: 'inherit', color: colors.textSecondary }}>{m.role}</div>
          <div>{m.gt}{m.partial ? ' ½' : ''}</div>
          <div>{m.alt != null || m.ref != null ? `${m.alt ?? '—'} / ${m.ref ?? '—'}` : m.cn != null ? `CN ${m.cn}` : '—'}</div>
          <div>{m.gq ?? '—'}</div>
        </div>
      ))}
    </div>
  )
}

function GeneBlock({ symbol }: { symbol: string }) {
  const [g, setG] = useState<GeneInfo | null>(null)
  const [missing, setMissing] = useState(false)
  useEffect(() => {
    let current = true
    setG(null)
    setMissing(false)
    structvarApi.gene(symbol).then((r) => current && setG(r)).catch(() => current && setMissing(true))
    return () => {
      current = false
    }
  }, [symbol])
  return (
    <Block title={`About ${symbol}`}>
      {missing && <div style={{ color: colors.textSecondary }}>No gene record.</div>}
      {!g && !missing && <div style={{ color: colors.textSecondary }}>Loading…</div>}
      {g && (
        <>
          <KV k="Transcript" v={`${g.transcript} · ${g.n_exons} exons · ${g.model_source}`} />
          <KV k="LoF constraint" v={g.pli != null ? `pLI ${g.pli.toFixed(2)}${g.oe_lof != null ? ` · observed/expected LoF ${g.oe_lof.toFixed(2)}` : ''} (gnomAD v2.1.1)` : 'not available'} />
          <KV k="Inheritance" v={g.inheritance.length ? g.inheritance.join(', ') : 'not recorded'} />
          <KV k="ClinVar conditions" v={g.diseases.length ? g.diseases.join('; ') : 'none recorded'} />
          <KV k="HPO phenotypes" v={g.hpo_n ? `${g.hpo_n} terms, e.g. ${g.hpo_terms.slice(0, 6).join('; ')}` : 'none recorded'} />
          <div style={{ color: colors.textTertiary, marginTop: 6, lineHeight: 1.5 }}>Gene-level context from PacBio&rsquo;s resource bundle (conditions as of 2022-12). Not a dosage-sensitivity curation.</div>
        </>
      )}
    </Block>
  )
}

function Block({ title, children }: { title: string; children: React.ReactNode }) {
  return (
    <div style={{ minWidth: 0 }}>
      <div style={{ fontSize: 11, fontWeight: 600, color: colors.textSecondary, textTransform: 'uppercase', letterSpacing: '0.04em', marginBottom: 6 }}>{title}</div>
      {children}
    </div>
  )
}

function KV({ k, v }: { k: string; v: string }) {
  return (
    <div style={{ display: 'flex', gap: 10, padding: '2px 0' }}>
      <div style={{ width: 112, flexShrink: 0, color: colors.textSecondary }}>{k}</div>
      <div style={{ wordBreak: 'break-word' }}>{v}</div>
    </div>
  )
}

// A filter combination that returns nothing says where the calls went.
function Empty({ resp }: { resp: SvsResponse }) {
  return (
    <div style={{ padding: '18px 20px' }}>
      <div style={{ fontSize: 13.5, fontWeight: 600 }}>No calls match these filters</div>
      {resp.funnel && (
        <div style={{ marginTop: 10, fontSize: 12.5 }}>
          <div style={{ color: colors.textSecondary, marginBottom: 6 }}>Calls left after each filter, in order:</div>
          {resp.funnel.map((f, i) => {
            const cliff = f.count === 0 && resp.funnel![i - 1]?.count !== 0
            return (
              <div key={f.label} style={{ display: 'flex', gap: 14, padding: '3px 0', color: cliff ? colors.statusPathogenic : colors.textBody, fontWeight: cliff ? 600 : 400 }}>
                <div style={{ width: 70, textAlign: 'right', fontFamily: mono }}>{fmtInt(f.count)}</div>
                <div>{f.label}</div>
              </div>
            )
          })}
        </div>
      )}
    </div>
  )
}
