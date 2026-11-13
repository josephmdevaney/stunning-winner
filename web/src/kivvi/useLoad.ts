import { useEffect, useRef, useState } from 'react'

// Loads a screen's data whenever `key` (the sample) changes. Responses that
// arrive after the key has changed again are dropped, so a fast sample switch
// can't show the previous sample's data.
export function useLoad<T>(load: () => Promise<T>, key: string): { data: T | null; error: string | null } {
  const loadRef = useRef(load)
  loadRef.current = load
  const [data, setData] = useState<T | null>(null)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    let current = true
    setData(null)
    setError(null)
    loadRef.current()
      .then((d) => current && setData(d))
      .catch((e) => current && setError(String(e)))
    return () => {
      current = false
    }
  }, [key])

  return { data, error }
}
