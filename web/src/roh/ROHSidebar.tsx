import { colors, mono } from '../theme'
import { ProductSwitcher, type Product } from '../components/ProductSwitcher'

export type ROHScreenKey = 'overview' | 'chromosome' | 'regions'

const NAV_ITEMS: { screen: ROHScreenKey; label: string }[] = [
  { screen: 'overview', label: 'Overview' },
  { screen: 'chromosome', label: 'Chromosome Detail' },
  { screen: 'regions', label: 'Regions' },
]

export function ROHSidebar({
  product,
  onSwitchProduct,
  screen,
  onNavigate,
  samples,
  sample,
  onSelectSample,
  genomeBuild,
  caller,
}: {
  product: Product
  onSwitchProduct: (p: Product) => void
  screen: ROHScreenKey
  onNavigate: (s: ROHScreenKey) => void
  samples: string[]
  sample: string
  onSelectSample: (s: string) => void
  genomeBuild?: string
  caller?: string
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
        {/* Chromosome mark: a bar with a centromere pinch and two ROH blocks. */}
        <svg width={28} height={28} viewBox="0 0 32 32">
          <rect x={4} y={13} width={24} height={6} rx={3} fill="none" stroke={colors.accent} strokeWidth={2} />
          <rect x={8} y={14.5} width={7} height={3} rx={1.5} fill={colors.accent} opacity={0.8} />
          <rect x={19} y={14.5} width={6} height={3} rx={1.5} fill={colors.accent} opacity={0.5} />
        </svg>
        <div style={{ color: 'white', fontWeight: 600, fontSize: 15, letterSpacing: '0.2px' }}>ROH</div>
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
        {caller && <div style={{ color: colors.sidebarTextSub, fontSize: 11.5, fontFamily: mono, marginTop: 6 }}>{caller}</div>}
        {genomeBuild && <div style={{ color: colors.sidebarTextMuted, fontSize: 11, fontFamily: mono, marginTop: 2 }}>{genomeBuild}</div>}
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
