import { colors } from '../theme'

// The original design showed a fabricated per-read CpG dot track. TRGT's VCF only
// reports one real number per allele (AM: mean 5mCpG methylation level) -- a per-read
// track would require the spanning-reads BAM, which this pipeline doesn't have. This
// renders the real per-allele value instead of inventing per-read detail.
function AlleleRow({ label, value, color }: { label: string; value: number | null; color: string }) {
  const pct = value != null ? Math.round(value * 100) : null
  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: 6 }}>
      <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: 12.5 }}>
        <span style={{ color: colors.textSecondary }}>{label}</span>
        <span style={{ fontWeight: 600, fontFamily: 'IBM Plex Mono' }}>{pct != null ? `${pct}%` : 'N/A'}</span>
      </div>
      <div style={{ height: 10, borderRadius: 5, background: colors.trackBg, overflow: 'hidden' }}>
        {pct != null && <div style={{ height: '100%', width: `${pct}%`, background: color }} />}
      </div>
    </div>
  )
}

export function MethylationPanel({
  a1Methylation,
  a2Methylation,
  genotyped,
}: {
  a1Methylation: number | null
  a2Methylation: number | null
  genotyped: boolean
}) {
  if (!genotyped) {
    return <div style={{ fontSize: 13, color: colors.textTertiary, padding: '20px 0' }}>No genotype call at this locus.</div>
  }
  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: 16, padding: '8px 0 4px' }}>
      <AlleleRow label="Allele 1 mean 5mCpG methylation" value={a1Methylation} color={colors.cluster1} />
      <AlleleRow label="Allele 2 mean 5mCpG methylation" value={a2Methylation} color={colors.cluster2} />
      {a1Methylation == null && a2Methylation == null && (
        <div style={{ fontSize: 12, color: colors.textTertiary }}>
          No CpG sites or methylation data available at this locus in the source VCF.
        </div>
      )}
    </div>
  )
}
