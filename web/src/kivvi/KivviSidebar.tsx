import { colors, mono } from './theme'
import { ProductSwitcher, type Product } from '../components/ProductSwitcher'
import type { SampleInfo } from './api'

export type KivviScreen = 'overview' | 'd4z4' | 'fshd2' | 'kiv2' | 'qc'

const NAV_ITEMS: { screen: KivviScreen; label: string }[] = [
  { screen: 'overview', label: 'Overview' },
  { screen: 'd4z4', label: 'D4Z4 (FSHD)' },
  { screen: 'fshd2', label: 'FSHD2 evidence' },
  { screen: 'kiv2', label: 'LPA KIV2' },
  { screen: 'qc', label: 'Sample QC' },
]

export function KivviSidebar({
  product,
  onSwitchProduct,
  screen,
  onNavigate,
  samples,
  sample,
  onSelectSample,
  kivviVersion,
  reference,
}: {
  product: Product
  onSwitchProduct: (p: Product) => void
  screen: KivviScreen
  onNavigate: (s: KivviScreen) => void
  samples: SampleInfo[]
  sample: string
  onSelectSample: (s: string) => void
  kivviVersion?: string | null
  reference?: string
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
        {/* Tandem-repeat mark: a row of repeat-unit blocks, one flagged. */}
        <svg width={28} height={28} viewBox="0 0 32 32">
          {[0, 1, 2, 3].map((i) => (
            <rect key={i} x={2 + i * 7.5} y={11} width={6} height={10} rx={1.5} fill={colors.accent} opacity={i === 2 ? 1 : 0.5} />
          ))}
        </svg>
        <div style={{ color: 'white', fontWeight: 600, fontSize: 15, letterSpacing: '0.2px' }}>Kivvi</div>
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
              <option key={s.key} value={s.key}>
                {s.label}
              </option>
            ))}
          </select>
        ) : (
          <div style={{ color: 'white', fontSize: 13.5, fontWeight: 600, marginTop: 4 }}>{sample || '—'}</div>
        )}
        {kivviVersion && <div style={{ color: colors.sidebarTextSub, fontSize: 11.5, fontFamily: mono, marginTop: 6 }}>Kivvi v{kivviVersion}</div>}
        {reference && <div style={{ color: colors.sidebarTextMuted, fontSize: 11, fontFamily: mono, marginTop: 2 }}>{reference}</div>}
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
