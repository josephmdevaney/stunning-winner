import { useEffect, useState } from 'react'
import { colors } from '../../theme'
import type { RohOverviewResponse, RohRegion } from '../api'
import { PerChromosomeTrackCard } from '../components'

export function ROHChromosomeScreen({ overview, allRegions }: { overview: RohOverviewResponse; allRegions: RohRegion[] }) {
  const [selChrom, setSelChrom] = useState('')

  // Default to the chromosome carrying a UPD flag, else the first with regions;
  // keep the user's choice if it is still valid for the newly loaded sample.
  useEffect(() => {
    const withRegions = overview.by_chromosome.filter((c) => c.centromere_mid_bp != null && c.n_regions > 0)
    if (withRegions.length === 0) return
    setSelChrom((prev) => {
      if (prev && withRegions.some((c) => c.chrom === prev)) return prev
      const updChrom = overview.upd_flags[0]?.chrom
      if (updChrom && withRegions.some((c) => c.chrom === updChrom)) return updChrom
      return withRegions[0].chrom
    })
  }, [overview])

  return (
    <div style={{ padding: 28, display: 'flex', flexDirection: 'column', gap: 20 }}>
      {selChrom ? (
        <PerChromosomeTrackCard overview={overview} allRegions={allRegions} selChrom={selChrom} onSelectChrom={setSelChrom} />
      ) : (
        <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, padding: '60px 24px', textAlign: 'center', color: colors.textSecondary, fontSize: 14 }}>
          No chromosome-level ROH regions to display for this sample.
        </div>
      )}
    </div>
  )
}
