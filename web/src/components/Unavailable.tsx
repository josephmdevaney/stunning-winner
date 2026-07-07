import { colors } from '../theme'
import type { UnavailableMetric } from '../api'

// Original design's QC screens simulated genome-wide/alignment metrics that a
// repeat-calling TRGT VCF simply doesn't contain. Rather than fabricate numbers,
// those panels report what's unavailable and why -- rendered here as a distinct,
// clearly-labeled block instead of pretending the design's chart is live data.
export function UnavailablePanel({ metrics }: { metrics: UnavailableMetric[] }) {
  if (metrics.length === 0) return null
  return (
    <div
      style={{
        background: 'oklch(0.98 0.005 90)',
        border: `1px dashed ${colors.borderInput}`,
        borderRadius: 12,
        padding: '16px 20px',
      }}
    >
      <div style={{ fontSize: 12, fontWeight: 600, color: colors.textSecondary, textTransform: 'uppercase', letterSpacing: '0.04em', marginBottom: 10 }}>
        Not available from this VCF
      </div>
      <div style={{ display: 'flex', flexDirection: 'column', gap: 8 }}>
        {metrics.map((m) => (
          <div key={m.label} style={{ fontSize: 13 }}>
            <span style={{ fontWeight: 600 }}>{m.label}</span>
            <span style={{ color: colors.textTertiary }}> — {m.reason}</span>
          </div>
        ))}
      </div>
    </div>
  )
}
