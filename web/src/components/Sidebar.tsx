import { colors } from '../theme'
import type { Screen } from '../RepeatScopeApp'
import { ProductSwitcher, type Product } from './ProductSwitcher'

const NAV_ITEMS: { screen: Screen; label: string }[] = [
  { screen: 'overview', label: 'Overview' },
  { screen: 'locus', label: 'Locus Detail' },
  { screen: 'roh', label: 'Runs of Homozygosity' },
  { screen: 'qc', label: 'Sample QC' },
  { screen: 'runqc', label: 'Run QC' },
  { screen: 'catalog', label: 'Repeat Catalog' },
]

export function Sidebar({
  product,
  onSwitchProduct,
  screen,
  onNavigate,
  availableScreens,
  sample,
  samples,
  onSelectSample,
  caller,
}: {
  product: Product
  onSwitchProduct: (p: Product) => void
  screen: Screen
  onNavigate: (s: Screen) => void
  availableScreens: Screen[]
  sample: string
  samples: string[]
  onSelectSample: (s: string) => void
  caller: string
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
        <svg width={28} height={28} viewBox="0 0 32 32">
          <polygon points="16,2 28,9 28,23 16,30 4,23 4,9" fill="none" stroke={colors.accent} strokeWidth={2.5} />
          <polygon points="16,10 22,13.5 22,20.5 16,24 10,20.5 10,13.5" fill={colors.accent} opacity={0.55} />
        </svg>
        <div style={{ color: 'white', fontWeight: 600, fontSize: 15, letterSpacing: '0.2px' }}>RepeatScope</div>
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

      {NAV_ITEMS.filter((item) => availableScreens.includes(item.screen)).map((item) => (
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
          <div style={{ color: 'white', fontSize: 13.5, fontWeight: 600, marginTop: 4 }}>{sample}</div>
        )}
        <div style={{ color: colors.sidebarTextSub, fontSize: 11.5, fontFamily: "'IBM Plex Mono'", marginTop: 6 }}>
          PacBio HiFi · {caller}
        </div>
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
