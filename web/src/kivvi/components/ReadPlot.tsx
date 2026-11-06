import { useState } from 'react'
import { colors } from '../theme'
import { kivviApi, type Target } from '../api'
import { SectionTitle } from './ui'
import { cardStyle } from './styles'

// Kivvi's own read-support plot (reads drawn onto each assembled allele). The SVG
// is several MB for D4Z4, so it's only fetched when asked for.
export function ReadPlot({ target, sample, available }: { target: Target; sample: string; available: boolean }) {
  const [open, setOpen] = useState(false)
  if (!available) return null
  return (
    <div style={cardStyle}>
      <SectionTitle title="Read support" sub="Kivvi's plot of reads laid onto each assembled allele." />
      {open ? (
        <div style={{ overflow: 'auto', maxHeight: 720, border: `1px solid ${colors.borderSoft}`, borderRadius: 8 }}>
          <img src={kivviApi.svgUrl(target, sample)} alt={`Kivvi ${target} read-support plot for ${sample}`} style={{ display: 'block', maxWidth: '100%' }} />
        </div>
      ) : (
        <button
          onClick={() => setOpen(true)}
          style={{ border: 'none', background: colors.accentButtonBg, color: 'white', borderRadius: 8, padding: '9px 16px', fontSize: 13, fontWeight: 600, cursor: 'pointer' }}
        >
          Load plot
        </button>
      )}
    </div>
  )
}
