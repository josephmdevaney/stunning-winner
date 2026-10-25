import { useEffect, useState } from 'react'
import { epigeneticsApi, type ImprintingOverview, type ThresholdOverrides, type Window } from './api'

// The overview for one sample. A refetch (new thresholds, other window) keeps the previous result on screen with
// `loading` set, so the table doesn't flash empty on every change.
export function useImprinting(sample: string, window: Window, overrides: ThresholdOverrides) {
  const [overview, setOverview] = useState<ImprintingOverview | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [loading, setLoading] = useState(false)
  const key = JSON.stringify(overrides)

  useEffect(() => {
    if (!sample) return
    let current = true
    setLoading(true)
    epigeneticsApi
      .imprinting(sample, window, overrides)
      .then((o) => {
        if (!current) return
        setOverview(o)
        setError(null)
      })
      .catch((e) => current && setError(String(e.message ?? e)))
      .finally(() => current && setLoading(false))
    return () => {
      current = false
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps -- `overrides` is compared by value through `key`
  }, [sample, window, key])

  return { overview, error, loading }
}
