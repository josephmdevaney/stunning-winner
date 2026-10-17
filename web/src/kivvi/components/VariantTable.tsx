import { useState } from 'react'
import { colors, fmtPct, mono } from '../theme'
import type { Variants } from '../api'
import { SectionTitle } from './ui'
import { cardStyle } from './styles'

const PAGE = 50

// Kivvi's small-variant calls. Columns are assembled alleles (Kivvi repurposes
// the VCF sample columns for that); a filled dot means the allele carries the
// variant. Variants on repeat units that weren't assembled into a complete
// allele carry no dot and are hidden unless asked for.
export function VariantTable({ variants }: { variants: Variants }) {
  const [onlyAssembled, setOnlyAssembled] = useState(true)
  const [shown, setShown] = useState(PAGE)

  const rows = onlyAssembled ? variants.records.filter((r) => r.on_assembled_allele) : variants.records
  const grid = `1.4fr 1.4fr ${variants.allele_labels.map(() => '1fr').join(' ')} 0.6fr 0.9fr`

  return (
    <div style={{ ...cardStyle, padding: 0, overflow: 'hidden' }}>
      <div style={{ padding: '16px 20px', borderBottom: `1px solid ${colors.borderSoft}`, display: 'flex', justifyContent: 'space-between', alignItems: 'center', gap: 16 }}>
        <SectionTitle
          title="Small variants"
          sub={`${variants.n_total} called · ${variants.n_on_assembled} on assembled alleles`}
        />
        <label style={{ fontSize: 12.5, color: colors.textSecondary, display: 'flex', alignItems: 'center', gap: 6, cursor: 'pointer' }}>
          <input
            type="checkbox"
            checked={onlyAssembled}
            onChange={(e) => {
              setOnlyAssembled(e.target.checked)
              setShown(PAGE)
            }}
          />
          Assembled alleles only
        </label>
      </div>

      <div style={{ display: 'grid', gridTemplateColumns: grid, gap: 10, padding: '10px 20px', fontSize: 10.5, fontWeight: 600, color: colors.textSecondary, letterSpacing: '0.04em', textTransform: 'uppercase', borderBottom: `1px solid ${colors.borderSofter}` }}>
        <div>Position</div>
        <div>Change</div>
        {variants.allele_labels.map((l) => (
          <div key={l}>{l}</div>
        ))}
        <div>Units</div>
        <div>Min support</div>
      </div>

      {rows.slice(0, shown).map((r, i) => (
        <div key={`${r.pos}-${r.ref}-${r.alt}-${i}`} style={{ display: 'grid', gridTemplateColumns: grid, gap: 10, padding: '8px 20px', fontSize: 12.5, alignItems: 'center', borderBottom: `1px solid ${colors.borderSoftest}`, fontFamily: mono }}>
          <div>
            {r.chrom}:{r.pos}
          </div>
          <div style={{ overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }} title={`${r.ref}>${r.alt}`}>
            {r.ref}→{r.alt}
          </div>
          {r.genotypes.map((g, j) => (
            <div key={j} style={{ color: g === '0' || g === '.' ? colors.textTertiary : colors.textPrimary }}>
              {g === '0' || g === '.' ? '○' : '●'}
            </div>
          ))}
          <div>{r.n_units}</div>
          <div>{fmtPct(r.min_vaf)}</div>
        </div>
      ))}

      {rows.length > shown && (
        <div style={{ padding: '12px 20px' }}>
          <button
            onClick={() => setShown((n) => n + PAGE)}
            style={{ border: `1px solid ${colors.borderInput}`, background: 'white', borderRadius: 6, padding: '6px 12px', fontSize: 12.5, cursor: 'pointer' }}
          >
            Show {Math.min(PAGE, rows.length - shown)} more ({rows.length - shown} remaining)
          </button>
        </div>
      )}
      {rows.length === 0 && <div style={{ padding: '16px 20px', color: colors.textSecondary, fontSize: 13 }}>No variants to show.</div>}
    </div>
  )
}
