import { useEffect, useState } from 'react'
import { rohApi, type RohRegionsResponse } from '../api'
import { RegionsTable } from '../components'

const PAGE_SIZE = 50

export function ROHRegionsScreen({ sample }: { sample: string }) {
  const [regions, setRegions] = useState<RohRegionsResponse | null>(null)
  const [chrom, setChrom] = useState('all')
  const [minLenMb, setMinLenMb] = useState('')
  const [page, setPage] = useState(0)

  useEffect(() => {
    let current = true
    setRegions(null)
    rohApi
      .regions({ sample, chrom, minLenMb: minLenMb ? Number(minLenMb) : undefined, page, pageSize: PAGE_SIZE })
      .then((r) => current && setRegions(r))
    return () => {
      current = false
    }
  }, [sample, chrom, minLenMb, page])

  return (
    <div style={{ padding: 28 }}>
      <RegionsTable
        regions={regions}
        chrom={chrom}
        minLenMb={minLenMb}
        onChangeChrom={(c) => {
          setChrom(c)
          setPage(0)
        }}
        onChangeMinLen={(v) => {
          setMinLenMb(v)
          setPage(0)
        }}
        page={page}
        pageSize={PAGE_SIZE}
        onChangePage={setPage}
      />
    </div>
  )
}
