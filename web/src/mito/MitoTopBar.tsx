import { colors } from './theme'

export function MitoTopBar({
  title,
  sample,
  runDate,
  reference,
}: {
  title: string
  sample: string
  runDate: string
  reference: string
}) {
  return (
    <div
      style={{
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'space-between',
        padding: '16px 28px',
        borderBottom: `1px solid ${colors.border}`,
        background: 'white',
      }}
    >
      <div>
        <div style={{ fontSize: 19, fontWeight: 600 }}>{title}</div>
        <div style={{ fontSize: 12.5, color: colors.textSecondary, marginTop: 2 }}>
          Sample {sample} · Run {runDate} · {reference}
        </div>
      </div>
      <div style={{ display: 'flex', alignItems: 'center', gap: 10 }}>
        <div
          style={{
            display: 'flex',
            alignItems: 'center',
            gap: 8,
            padding: '8px 14px',
            border: `1px solid ${colors.borderInputAlt}`,
            borderRadius: 8,
            fontSize: 13,
            color: colors.textSecondary,
          }}
        >
          <svg width={14} height={14} viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={2}>
            <circle cx={11} cy={11} r={7} />
            <line x1={21} y1={21} x2={16.65} y2={16.65} />
          </svg>
          Search variant or locus
        </div>
        <div
          style={{
            padding: '9px 16px',
            borderRadius: 8,
            background: colors.accentButtonBg,
            color: 'white',
            fontSize: 13,
            fontWeight: 600,
            cursor: 'pointer',
          }}
        >
          Export Report
        </div>
      </div>
    </div>
  )
}
