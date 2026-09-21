import { colors, mono } from '../theme'
import { ProductSwitcher, type Product } from '../components/ProductSwitcher'

export type ParaphaseScreen = 'overview' | 'region'

const NAV_ITEMS: { screen: ParaphaseScreen; label: string }[] = [
  { screen: 'overview', label: 'Overview' },
  { screen: 'region', label: 'Region Detail' },
]

export function ParaphaseSidebar({
  product,
  onSwitchProduct,
  screen,
  onNavigate,
  regionSelected,
  samples,
  sample,
  onSelectSample,
}: {
  product: Product
  onSwitchProduct: (p: Product) => void
  screen: ParaphaseScreen
  onNavigate: (s: ParaphaseScreen) => void
  regionSelected: boolean
  samples: string[]
  sample: string
  onSelectSample: (s: string) => void
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
        {/* Branching/paralog mark: one trunk splitting into gene-copy haplotypes. */}
        <svg width={28} height={28} viewBox="0 0 32 32">
          <path d="M16 28 V16" stroke={colors.accent} strokeWidth={2.5} fill="none" />
          <path d="M16 16 L8 6" stroke={colors.accent} strokeWidth={2.5} fill="none" />
          <path d="M16 16 L24 6" stroke={colors.accent} strokeWidth={2.5} fill="none" />
          <circle cx={16} cy={28} r={2.5} fill={colors.accent} />
          <circle cx={8} cy={5} r={2.5} fill={colors.accent} opacity={0.55} />
          <circle cx={24} cy={5} r={2.5} fill={colors.accent} opacity={0.55} />
        </svg>
        <div style={{ color: 'white', fontWeight: 600, fontSize: 15, letterSpacing: '0.2px' }}>Paraphase</div>
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
        Sample
      </div>

      {NAV_ITEMS.map((item) => (
        <NavRow
          key={item.screen}
          active={screen === item.screen}
          label={item.label}
          disabled={item.screen === 'region' && !regionSelected}
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
              <option key={s} value={s}>
                {s}
              </option>
            ))}
          </select>
        ) : (
          <div style={{ color: 'white', fontSize: 13.5, fontWeight: 600, marginTop: 4 }}>{sample || '—'}</div>
        )}
        <div style={{ color: colors.sidebarTextSub, fontSize: 11.5, fontFamily: mono, marginTop: 6 }}>
          interpret_paraphase.py · 14 medical regions
        </div>
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
      <div
        style={{
          width: 6,
          height: 6,
          borderRadius: 2,
          flexShrink: 0,
          background: active ? colors.accent : colors.sidebarDotInactive,
        }}
      />
      {label}
    </div>
  )
}
