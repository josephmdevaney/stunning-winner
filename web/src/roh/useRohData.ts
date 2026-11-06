import { useEffect, useState } from 'react'
import { rohApi, type RohOverviewResponse, type RohRegion } from './api'

// Loads a sample's ROH overview and its full region list once, so every screen
// of the suite shares them (the ideogram and chromosome track both need every
// region, not just a page of the regions table).
export function useRohData(sample: string): {
  overview: RohOverviewResponse | null
  allRegions: RohRegion[]
  error: string | null
} {
  const [overview, setOverview] = useState<RohOverviewResponse | null>(null)
  const [allRegions, setAllRegions] = useState<RohRegion[]>([])
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    if (!sample) return
    let current = true
    setOverview(null)
    setAllRegions([])
    setError(null)
    rohApi
      .overview(sample)
      .then(async (o) => {
        if (!current) return
        setOverview(o)
        const r = await rohApi.regions({ sample, chrom: 'all', minLenMb: 0, page: 0, pageSize: Math.max(1, o.stats.n_regions_total) })
        if (current) setAllRegions(r.rows)
      })
      .catch((e) => current && setError(String(e)))
    return () => {
      current = false
    }
  }, [sample])

  return { overview, allRegions, error }
}
