import { useState } from 'react'
import { Card, CardHeading } from '../../components/Card'
import { colors, mono } from '../../theme'
import type { ImprintingOverview, ThresholdOverrides, Window } from '../api'
import { CATEGORY_HELP, CATEGORY_LABEL, STATUS_COLOR, STATUS_LABEL } from '../labels'
import { Chip } from '../ui'

const DEFAULTS = { loss_below: 20, gain_above: 80, allelic_min_diff: 40 }

export function MethodsScreen({
  overview,
  window,
  overrides,
  onChangeWindow,
  onChangeOverrides,
}: {
  overview: ImprintingOverview
  window: Window
  overrides: ThresholdOverrides
  onChangeWindow: (w: Window) => void
  onChangeOverrides: (o: ThresholdOverrides) => void
}) {
  const th = overview.thresholds
  const [draft, setDraft] = useState({ loss_below: String(th.loss_below), gain_above: String(th.gain_above), allelic_min_diff: String(th.allelic_min_diff) })
  const num = (s: string) => (s.trim() === '' ? NaN : Number(s))
  const parsed = { loss_below: num(draft.loss_below), gain_above: num(draft.gain_above), allelic_min_diff: num(draft.allelic_min_diff) }
  const valid =
    Object.values(parsed).every((v) => Number.isFinite(v)) && parsed.loss_below >= 0 && parsed.loss_below < parsed.gain_above && parsed.gain_above <= 100 && parsed.allelic_min_diff > 0 && parsed.allelic_min_diff <= 100
  const isCustom = th.loss_below !== DEFAULTS.loss_below || th.gain_above !== DEFAULTS.gain_above || th.allelic_min_diff !== DEFAULTS.allelic_min_diff
  const pv = overview.provenance

  const field = (key: keyof typeof draft, label: string, help: string) => (
    <label style={{ display: 'block' }}>
      <div style={{ fontSize: 12.5, fontWeight: 600 }}>{label}</div>
      <input
        value={draft[key]}
        onChange={(e) => setDraft({ ...draft, [key]: e.target.value })}
        inputMode="decimal"
        style={{ width: 90, marginTop: 6, padding: '6px 8px', fontSize: 13.5, fontFamily: mono, borderRadius: 6, border: `1px solid ${colors.borderInput}` }}
      />
      <div style={{ fontSize: 11.5, color: colors.textTertiary, marginTop: 4, maxWidth: 220, lineHeight: 1.45 }}>{help}</div>
    </label>
  )

  return (
    <div style={{ padding: 28, display: 'flex', flexDirection: 'column', gap: 20 }}>
      <Card>
        <CardHeading
          title="Thresholds"
          sub="These are research defaults chosen for this panel, not clinically validated cut-offs. Change them to see how a call depends on them."
        />
        <div style={{ display: 'flex', gap: 28, flexWrap: 'wrap', alignItems: 'flex-start' }}>
          {field('loss_below', 'Loss below (%)', 'Both alleles under this are a loss of methylation.')}
          {field('gain_above', 'Gain above (%)', 'Both alleles over this are a gain of methylation.')}
          {field('allelic_min_diff', 'Allelic difference (points)', 'Haplotypes at least this far apart are an allele-specific imprint.')}
          <div style={{ display: 'flex', flexDirection: 'column', gap: 8, marginTop: 22 }}>
            <button
              disabled={!valid}
              onClick={() => onChangeOverrides(parsed)}
              style={{ padding: '7px 16px', borderRadius: 6, border: 'none', background: valid ? colors.accentButtonBg : colors.statusUnknown, color: 'white', fontWeight: 600, fontSize: 12.5, cursor: valid ? 'pointer' : 'default' }}
            >
              Apply
            </button>
            <button
              disabled={!isCustom && !Object.keys(overrides).length}
              onClick={() => {
                setDraft({ loss_below: '20', gain_above: '80', allelic_min_diff: '40' })
                onChangeOverrides({})
              }}
              style={{ padding: '6px 16px', borderRadius: 6, border: `1px solid ${colors.borderInput}`, background: 'white', fontSize: 12.5, cursor: 'pointer' }}
            >
              Reset
            </button>
          </div>
        </div>
        <div style={{ fontSize: 12, color: colors.textSecondary, marginTop: 16, lineHeight: 1.6 }}>
          Fixed: at least {th.min_cpgs} CpGs and {th.min_mean_coverage}× mean depth ({th.min_hap_mean_coverage}× per haplotype) to make a call · unphased regions fall back to bands (intact imprint {th.imprint_low}–{th.imprint_high}%).
          {isCustom && <strong style={{ color: STATUS_COLOR.borderline }}> Custom thresholds are in use.</strong>}
        </div>
        {!valid && <div style={{ fontSize: 12, color: colors.statusPathogenic, marginTop: 8 }}>Need 0 ≤ loss &lt; gain ≤ 100 and a positive allelic difference.</div>}
      </Card>

      <Card>
        <CardHeading title="Analysis window" sub="Where the methylation is pooled for each DMR." />
        <div style={{ display: 'flex', gap: 8 }}>
          {(['core', 'dmr'] as const).map((w) => (
            <button
              key={w}
              onClick={() => onChangeWindow(w)}
              style={{
                padding: '7px 14px',
                borderRadius: 6,
                border: `1px solid ${window === w ? colors.accentButtonBg : colors.borderInput}`,
                background: window === w ? colors.accentButtonBg : 'white',
                color: window === w ? 'white' : colors.textPrimary,
                fontSize: 12.5,
                fontWeight: 600,
                cursor: 'pointer',
              }}
            >
              {w === 'core' ? 'Core (recommended)' : 'Whole published DMR'}
            </button>
          ))}
        </div>
        <div style={{ fontSize: 12, color: colors.textSecondary, marginTop: 10, lineHeight: 1.55 }}>
          Core windows are the short, well-defined stretch of the DMR used for the seven clinical loci (for example the SNRPN promoter CpG island). Elsewhere, and in “Whole published DMR”, the whole DMR is used.
        </div>
      </Card>

      <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 20 }}>
        <Card>
          <CardHeading title="How a region is called" />
          <div style={{ display: 'flex', flexDirection: 'column', gap: 10, fontSize: 12.5, lineHeight: 1.5 }}>
            {(['imprinted', 'loss', 'gain', 'borderline', 'low_coverage'] as const).map((s) => (
              <div key={s} style={{ display: 'grid', gridTemplateColumns: '140px 1fr', gap: 10, alignItems: 'start' }}>
                <Chip color={STATUS_COLOR[s]} solid={s !== 'low_coverage'}>
                  {STATUS_LABEL[s]}
                </Chip>
                <span style={{ color: colors.textBody }}>{STATUS_MEANING[s]}</span>
              </div>
            ))}
          </div>
        </Card>
        <Card>
          <CardHeading title="Region categories" sub="From the 2024 imprinting-disorder consensus." />
          <div style={{ display: 'flex', flexDirection: 'column', gap: 10, fontSize: 12.5, lineHeight: 1.5 }}>
            {(['CA', 'NC', 'X'] as const).map((c) => (
              <div key={c}>
                <strong>{CATEGORY_LABEL[c]}</strong>
                <div style={{ color: colors.textBody }}>{CATEGORY_HELP[c]}</div>
              </div>
            ))}
          </div>
        </Card>
      </div>

      <Card>
        <CardHeading title="Data and panel" />
        <dl style={{ display: 'grid', gridTemplateColumns: '190px 1fr', gap: '8px 16px', margin: 0, fontSize: 12.5 }}>
          <Term k="Sample" v={`${overview.label}${overview.note ? ` — ${overview.note}` : ''}`} />
          <Term k="Methylation caller" v={pv.pb_cpg_tools_version ? `pb-CpG-tools ${pv.pb_cpg_tools_version}` : 'version not recorded'} />
          <Term k="Pileup mode" v={[pv.pileup_mode, pv.modsites_mode].filter(Boolean).join(' · ') || 'not recorded'} />
          <Term k="Minimum site coverage" v={pv.min_coverage != null ? `${pv.min_coverage}×` : 'not recorded'} />
          <Term k="Tracks" v={pv.tracks.join(', ')} />
          <Term k="Reference" v={pv.genome_build} />
          <Term k="DMR panel" v={`${overview.panel.title} (${overview.panel.journal}, ${overview.panel.pmcid}, ${overview.panel.license})`} />
          <Term k="Pooling" v="Σ modified calls ÷ Σ (modified + unmodified calls) across the CpGs in the window, so deeper sites count for more." />
        </dl>
      </Card>

      <Card>
        <CardHeading title="What this does not tell you" />
        <div style={{ fontSize: 12.5, color: colors.textBody, lineHeight: 1.6 }}>{overview.caveats}</div>
      </Card>
    </div>
  )
}

const STATUS_MEANING = {
  imprinted: 'The two haplotypes differ by at least the allelic threshold, or (unphased) the combined value is in the intact band.',
  loss: 'Both haplotypes, or the combined value, are below the loss threshold.',
  gain: 'Both haplotypes, or the combined value, are above the gain threshold.',
  borderline: 'Between the bands, or haplotypes that disagree with the combined value. Look at the CpG-level plot.',
  low_coverage: 'Too few CpGs or too little depth for a call.',
} as const

function Term({ k, v }: { k: string; v: string }) {
  return (
    <>
      <dt style={{ color: colors.textSecondary, fontWeight: 600 }}>{k}</dt>
      <dd style={{ margin: 0 }}>{v}</dd>
    </>
  )
}
