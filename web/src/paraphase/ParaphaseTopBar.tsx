import { colors } from '../theme'

export function ParaphaseTopBar({ title, sample }: { title: string; sample: string }) {
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
          Sample {sample} · Paraphase Annotation · GRCh38
        </div>
      </div>
    </div>
  )
}
