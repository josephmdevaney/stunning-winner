import { colors, fmtPct } from '../theme'

const UNIT_W = 12

// Distinct repeat-unit IDs get distinct hues so unit-level structure (e.g. a run
// of identical units) is visible at a glance.
function unitColor(id: string): string {
  const n = parseInt(id, 10)
  const hue = Number.isFinite(n) ? (Math.abs(n) * 47) % 360 : 0
  return `oklch(0.78 0.1 ${hue})`
}

function methylColor(v: number): string {
  return `color-mix(in oklch, ${colors.methylated} ${Math.round(v * 100)}%, ${colors.unmethylated})`
}

function Row({ label, children }: { label: string; children: React.ReactNode }) {
  return (
    <div style={{ display: 'flex', alignItems: 'center', gap: 10 }}>
      <div style={{ width: 70, flexShrink: 0, fontSize: 10.5, fontWeight: 600, color: colors.textSecondary, textTransform: 'uppercase', letterSpacing: '0.04em' }}>
        {label}
      </div>
      <div style={{ display: 'flex', flexWrap: 'wrap', gap: 2 }}>{children}</div>
    </div>
  )
}

// One allele drawn as its ordered repeat units: a type strip (unit ID) and, when
// available, a methylation strip (per-unit methylated-read fraction).
export function AlleleStrip({ units, methylation }: { units: string[]; methylation?: number[] }) {
  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: 4 }}>
      <Row label="Unit type">
        {units.map((id, i) => (
          <div key={i} title={`Unit ${i + 1} · repeat unit ID ${id}`} style={{ width: UNIT_W, height: 14, borderRadius: 2, background: unitColor(id) }} />
        ))}
      </Row>
      {methylation && methylation.length > 0 && (
        <Row label="Methylation">
          {methylation.map((v, i) => (
            <div key={i} title={`Unit ${i + 1} · ${fmtPct(v)} methylated`} style={{ width: UNIT_W, height: 14, borderRadius: 2, background: methylColor(v) }} />
          ))}
        </Row>
      )}
    </div>
  )
}
