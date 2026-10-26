import { colors, mono } from '../theme'
import { ProductSwitcher, type Product } from '../components/ProductSwitcher'
import type { SampleInfo } from './api'

export type EpigeneticsScreen = 'overview' | 'region' | 'episignatures' | 'signature' | 'methods'

const NAV_ITEMS: { screen: EpigeneticsScreen; label: string }[] = [
  { screen: 'overview', label: 'Imprinting' },
  { screen: 'region', label: 'Region Detail' },
  { screen: 'episignatures', label: 'Episignatures' },
  { screen: 'signature', label: 'Signature Detail' },
  { screen: 'methods', label: 'Methods' },
]

export function EpigeneticsSidebar({
  product,
  onSwitchProduct,
  screen,
  onNavigate,
  regionSelected,
  signatureSelected,
  samples,
  sample,
  onSelectSample,
  footer,
}: {
  product: Product
  onSwitchProduct: (p: Product) => void
  screen: EpigeneticsScreen
  onNavigate: (s: EpigeneticsScreen) => void
  regionSelected: boolean
  signatureSelected: boolean
  samples: SampleInfo[]
  sample: string
  onSelectSample: (s: string) => void
  footer: string
}) {
  return (
    <div
      style={{
        width: 232,
        flexShrink: 0,
        background: colors.sidebarBg,
        display: 'flex',
        flexDirection: 'column',
        padding: '20px 14px',
        gap: 4,
      }}
    >
      <ProductSwitcher product={product} onChange={onSwitchProduct} />

      <div style={{ display: 'flex', alignItems: 'center', gap: 10, padding: '8px 8px 20px 8px' }}>
        {/* Methylation mark: a strand with two filled (methylated) and two open (unmethylated) CpG lollipops. */}
        <svg width={28} height={28} viewBox="0 0 32 32">
          <line x1={3} y1={22} x2={29} y2={22} stroke={colors.accent} strokeWidth={2.5} />
          <line x1={8} y1={22} x2={8} y2={12} stroke={colors.accent} strokeWidth={2} />
          <line x1={14} y1={22} x2={14} y2={8} stroke={colors.accent} strokeWidth={2} />
          <line x1={20} y1={22} x2={20} y2={13} stroke={colors.accent} strokeWidth={2} />
          <line x1={26} y1={22} x2={26} y2={9} stroke={colors.accent} strokeWidth={2} />
          <circle cx={8} cy={9} r={3} fill={colors.accent} />
          <circle cx={14} cy={5} r={3} fill="none" stroke={colors.accent} strokeWidth={1.8} />
          <circle cx={20} cy={10} r={3} fill={colors.accent} />
          <circle cx={26} cy={6} r={3} fill="none" stroke={colors.accent} strokeWidth={1.8} />
        </svg>
        <div style={{ color: 'white', fontWeight: 600, fontSize: 15, letterSpacing: '0.2px' }}>Epigenetics</div>
      </div>

      <div
        style={{
          color: colors.sidebarTextMuted,
          fontSize: 11,
          fontWeight: 600,
          letterSpacing: '0.08em',
          textTransform: 'uppercase',
          padding: '10px 10px 6px',
        }}
      >
        Genome
      </div>

      {NAV_ITEMS.map((item) => (
        <NavRow
          key={item.screen}
          active={screen === item.screen}
          label={item.label}
          disabled={(item.screen === 'region' && !regionSelected) || (item.screen === 'signature' && !signatureSelected)}
          onClick={() => onNavigate(item.screen)}
        />
      ))}

      <div style={{ marginTop: 'auto', padding: '14px 12px', borderRadius: 10, background: colors.sidebarPanelBg }}>
        <div style={{ color: colors.sidebarTextMuted, fontSize: 10.5, fontWeight: 600, letterSpacing: '0.06em', textTransform: 'uppercase' }}>
          Active sample
        </div>
        {samples.length > 1 ? (
          <select
            value={sample}
            onChange={(e) => onSelectSample(e.target.value)}
            style={{
              width: '100%',
              marginTop: 4,
              background: colors.sidebarBg,
              color: 'white',
              fontSize: 13.5,
              fontWeight: 600,
              border: `1px solid ${colors.sidebarHover}`,
              borderRadius: 6,
              padding: '4px 6px',
              cursor: 'pointer',
            }}
          >
            {samples.map((s) => (
              <option key={s.key} value={s.key}>
                {s.label}
              </option>
            ))}
          </select>
        ) : (
          <div style={{ color: 'white', fontSize: 13.5, fontWeight: 600, marginTop: 4 }}>{samples[0]?.label ?? '—'}</div>
        )}
        <div style={{ color: colors.sidebarTextSub, fontSize: 11.5, fontFamily: mono, marginTop: 6 }}>{footer}</div>
      </div>
    </div>
  )
}

function NavRow({ active, label, disabled, onClick }: { active: boolean; label: string; disabled?: boolean; onClick: () => void }) {
  return (
    <div
      onClick={disabled ? undefined : onClick}
      style={{
        display: 'flex',
        alignItems: 'center',
        gap: 10,
        padding: '10px 12px',
        borderRadius: 8,
        cursor: disabled ? 'default' : 'pointer',
        opacity: disabled ? 0.45 : 1,
        color: colors.sidebarTextItem,
        fontSize: 13.5,
        fontWeight: 500,
      }}
      onMouseEnter={(e) => {
        if (!disabled) e.currentTarget.style.background = colors.sidebarHover
      }}
      onMouseLeave={(e) => (e.currentTarget.style.background = 'transparent')}
    >
      <div style={{ width: 6, height: 6, borderRadius: 2, flexShrink: 0, background: active ? colors.accent : colors.sidebarDotInactive }} />
      {label}
    </div>
  )
}
