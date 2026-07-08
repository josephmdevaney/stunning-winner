import { colors, mono } from './theme'
import { ProductSwitcher, type Product } from '../components/ProductSwitcher'

export type MitoScreen = 'overview' | 'detail' | 'qc'

const NAV_ITEMS: { screen: MitoScreen; label: string }[] = [
  { screen: 'overview', label: 'Overview' },
  { screen: 'detail', label: 'Variant Detail' },
  { screen: 'qc', label: 'Sample QC' },
]

export function MitoSidebar({
  product,
  onSwitchProduct,
  screen,
  onNavigate,
  sample,
  instrumentLine,
}: {
  product: Product
  onSwitchProduct: (p: Product) => void
  screen: MitoScreen
  onNavigate: (s: MitoScreen) => void
  sample: string
  instrumentLine: string
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
        {/* Mitochondrion-style mark: outer ring + filled inner circle + top dot. */}
        <svg width={28} height={28} viewBox="0 0 32 32">
          <circle cx={16} cy={16} r={13} fill="none" stroke={colors.accent} strokeWidth={2.5} />
          <circle cx={16} cy={16} r={6.5} fill={colors.accent} opacity={0.55} />
          <circle cx={16} cy={3} r={2} fill={colors.accent} />
        </svg>
        <div style={{ color: 'white', fontWeight: 600, fontSize: 15, letterSpacing: '0.2px' }}>MitoLR</div>
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
        <NavRow key={item.screen} active={screen === item.screen} label={item.label} onClick={() => onNavigate(item.screen)} />
      ))}

      <div style={{ marginTop: 'auto', padding: '14px 12px', borderRadius: 10, background: colors.sidebarPanelBg }}>
        <div style={{ color: colors.sidebarTextMuted, fontSize: 10.5, fontWeight: 600, letterSpacing: '0.06em', textTransform: 'uppercase' }}>
          Active sample
        </div>
        <div style={{ color: 'white', fontSize: 13.5, fontWeight: 600, marginTop: 4 }}>{sample}</div>
        <div style={{ color: colors.sidebarTextSub, fontSize: 11.5, fontFamily: mono, marginTop: 2 }}>{instrumentLine}</div>
      </div>
    </div>
  )
}

function NavRow({ active, label, onClick }: { active: boolean; label: string; onClick: () => void }) {
  return (
    <div
      onClick={onClick}
      style={{
        display: 'flex',
        alignItems: 'center',
        gap: 10,
        padding: '10px 12px',
        borderRadius: 8,
        cursor: 'pointer',
        color: colors.sidebarTextItem,
        fontSize: 13.5,
        fontWeight: 500,
      }}
      onMouseEnter={(e) => (e.currentTarget.style.background = colors.sidebarHover)}
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
