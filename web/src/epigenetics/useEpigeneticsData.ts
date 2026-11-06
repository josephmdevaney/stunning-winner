import { useEffect, useState } from 'react'
import { epigeneticsApi, type ArrayResponse, type EpisignatureResponse, type ImprintingOverview, type ThresholdOverrides, type Window } from './api'

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

// EpiSignaLR scores for one sample. Fetched only once the Episignatures screens are opened (`enabled`), because the
// first request reads a few thousand CpGs and the imprinting screens do not need it.
export function useEpisignatures(sample: string, enabled: boolean) {
  const [data, setData] = useState<EpisignatureResponse | null>(null)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    if (!sample || !enabled) return
    let current = true
    setData(null)
    setError(null)
    epigeneticsApi
      .episignatures(sample)
      .then((d) => current && setData(d))
      .catch((e) => current && setError(String(e.message ?? e)))
    return () => {
      current = false
    }
  }, [sample, enabled])

  return { data, error }
}

// EpiSignaLR's array-derived scores, fetched with the other episignature screens. Slower than they are (every
// genome retrains 34 SVMs), so it has its own state and the sequencing-derived scores do not wait for it.
export function useArraySignatures(sample: string, enabled: boolean) {
  const [data, setData] = useState<ArrayResponse | null>(null)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    if (!sample || !enabled) return
    let current = true
    setData(null)
    setError(null)
    epigeneticsApi
      .arraySignatures(sample)
      .then((d) => current && setData(d))
      .catch((e) => current && setError(String(e.message ?? e)))
    return () => {
      current = false
    }
  }, [sample, enabled])

  return { data, error }
}
