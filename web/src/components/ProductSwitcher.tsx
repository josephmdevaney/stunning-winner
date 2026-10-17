import { colors } from '../theme'

export type Product = 'genome_wide' | 'puretarget_v2' | 'carrier' | 'mitolr' | 'paraphase' | 'kivvi'

const PRODUCTS: { key: Product; label: string }[] = [
  { key: 'genome_wide', label: 'RepeatScope' },
  { key: 'puretarget_v2', label: 'PureTarget V2' },
  { key: 'carrier', label: 'Carrier' },
  { key: 'mitolr', label: 'MitoLR' },
  { key: 'paraphase', label: 'Paraphase' },
  { key: 'kivvi', label: 'Kivvi' },
]

// Top-level section switcher pinned to every sidebar. All sections share
// one deployable and one visual system; this is how you move between them.
export function ProductSwitcher({ product, onChange }: { product: Product; onChange: (p: Product) => void }) {
  return (
    <div style={{ padding: '0 8px 12px' }}>
      <div
        style={{
          color: colors.sidebarTextMuted,
          fontSize: 10.5,
          fontWeight: 600,
          letterSpacing: '0.08em',
          textTransform: 'uppercase',
          marginBottom: 6,
        }}
      >
        Suite
      </div>
      <div style={{ display: 'flex', flexDirection: 'column', gap: 2, padding: 4, background: colors.sidebarPanelBg, borderRadius: 8 }}>
        {PRODUCTS.map((p) => (
          <button
            key={p.key}
            onClick={() => onChange(p.key)}
            style={{
              padding: '7px 10px',
              borderRadius: 6,
              border: 'none',
              cursor: 'pointer',
              textAlign: 'left',
              fontSize: 12.5,
              fontWeight: 600,
              fontFamily: "'IBM Plex Sans', sans-serif",
              background: product === p.key ? colors.accentButtonBg : 'transparent',
              color: product === p.key ? 'white' : colors.sidebarTextItem,
            }}
          >
            {p.label}
          </button>
        ))}
      </div>
    </div>
  )
}
