import { colors, mono } from './theme'
import { ProductSwitcher, type Product } from '../components/ProductSwitcher'
import type { SampleInfo } from './api'

export type StructVarScreen = 'overview' | 'svs'

const NAV_ITEMS: { screen: StructVarScreen; label: string }[] = [
  { screen: 'overview', label: 'Overview' },
  { screen: 'svs', label: 'Structural variants' },
]

export function StructVarSidebar({
  product, onSwitchProduct, screen, onNavigate, samples, sample, onSelectSample, reference,
}: {
  product: Product
  onSwitchProduct: (p: Product) => void
  screen: StructVarScreen
  onNavigate: (s: StructVarScreen) => void
  samples: SampleInfo[]
  sample: string
  onSelectSample: (s: string) => void
  reference?: string
}) {
  return (
    <div style={{ width: 232, flexShrink: 0, background: colors.sidebarBg, display: 'flex', flexDirection: 'column', padding: '20px 14px', gap: 4 }}>
      <ProductSwitcher product={product} onChange={onSwitchProduct} />

      <div style={{ display: 'flex', alignItems: 'center', gap: 10, padding: '8px 8px 20px 8px' }}>
        {/* A strand with a deleted stretch and an inserted block. */}
        <svg width={28} height={28} viewBox="0 0 32 32">
          <line x1={2} y1={18} x2={11} y2={18} stroke={colors.accent} strokeWidth={2.5} />
          <line x1={14} y1={18} x2={19} y2={18} stroke={colors.accent} strokeWidth={2.5} strokeDasharray="1 3" opacity={0.6} />
          <line x1={22} y1={18} x2={30} y2={18} stroke={colors.accent} strokeWidth={2.5} />
          <rect x={22} y={8} width={7} height={6} rx={1.5} fill={colors.accent} opacity={0.9} />
        </svg>
        <div style={{ color: 'white', fontWeight: 600, fontSize: 15, letterSpacing: '0.2px' }}>Structural variants</div>
      </div>

      <div style={{ color: colors.sidebarTextMuted, fontSize: 11, fontWeight: 600, letterSpacing: '0.08em', textTransform: 'uppercase', padding: '10px 10px 6px' }}>
        Call set
      </div>

      {NAV_ITEMS.map((item) => (
        <NavRow key={item.screen} active={screen === item.screen} label={item.label} onClick={() => onNavigate(item.screen)} />
      ))}

      <div style={{ marginTop: 'auto', padding: '14px 12px', borderRadius: 10, background: colors.sidebarPanelBg }}>
        <div style={{ color: colors.sidebarTextMuted, fontSize: 10.5, fontWeight: 600, letterSpacing: '0.06em', textTransform: 'uppercase' }}>Active call set</div>
        {samples.length > 1 ? (
          <select
            value={sample}
            onChange={(e) => onSelectSample(e.target.value)}
            style={{ width: '100%', marginTop: 4, background: colors.sidebarBg, color: 'white', fontSize: 13.5, fontWeight: 600, border: `1px solid ${colors.sidebarHover}`, borderRadius: 6, padding: '4px 6px', cursor: 'pointer' }}
          >
            {samples.map((s) => (
              <option key={s.key} value={s.key}>
                {s.label}
              </option>
            ))}
          </select>
        ) : (
          <div style={{ color: 'white', fontSize: 13.5, fontWeight: 600, marginTop: 4 }}>{samples[0]?.label ?? (sample || '—')}</div>
        )}
        {reference && <div style={{ color: colors.sidebarTextMuted, fontSize: 11, fontFamily: mono, marginTop: 6 }}>{reference}</div>}
      </div>
    </div>
  )
}

function NavRow({ active, label, onClick }: { active: boolean; label: string; onClick: () => void }) {
  return (
    <div
      onClick={onClick}
      style={{ display: 'flex', alignItems: 'center', gap: 10, padding: '10px 12px', borderRadius: 8, cursor: 'pointer', color: colors.sidebarTextItem, fontSize: 13.5, fontWeight: 500 }}
      onMouseEnter={(e) => (e.currentTarget.style.background = colors.sidebarHover)}
      onMouseLeave={(e) => (e.currentTarget.style.background = 'transparent')}
    >
      <div style={{ width: 6, height: 6, borderRadius: 2, flexShrink: 0, background: active ? colors.accent : colors.sidebarDotInactive }} />
      {label}
    </div>
  )
}
