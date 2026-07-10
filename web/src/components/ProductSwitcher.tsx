import { colors } from '../theme'

export type Product = 'repeatscope' | 'mitolr'

const PRODUCTS: { key: Product; label: string }[] = [
  { key: 'repeatscope', label: 'RepeatScope' },
  { key: 'mitolr', label: 'MitoLR' },
]

// Segmented control pinned to the top of either product's sidebar. The two tools
// share one deployable and one visual system; this is how you move between them.
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
      <div style={{ display: 'flex', gap: 4, padding: 4, background: colors.sidebarPanelBg, borderRadius: 8 }}>
        {PRODUCTS.map((p) => (
          <button
            key={p.key}
            onClick={() => onChange(p.key)}
            style={{
              flex: 1,
              padding: '6px 4px',
              borderRadius: 6,
              border: 'none',
              cursor: 'pointer',
              fontSize: 12,
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
