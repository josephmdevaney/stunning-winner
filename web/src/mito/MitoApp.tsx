import { useEffect, useState } from 'react'
import { MitoSidebar, type MitoScreen } from './MitoSidebar'
import { MitoTopBar } from './MitoTopBar'
import type { Product } from '../components/ProductSwitcher'
import { colors } from './theme'
import { mitoApi, type MitoMeta, type SampleInfo } from './api'
import { MitoOverviewScreen } from './screens/MitoOverviewScreen'
import { MitoVariantDetailScreen } from './screens/MitoVariantDetailScreen'
import { MitoQCScreen } from './screens/MitoQCScreen'

const TITLES: Record<MitoScreen, string> = {
  overview: 'Genome-wide Overview',
  detail: 'Variant Detail',
  qc: 'Sample QC Dashboard',
}

export default function MitoApp({
  product,
  onSwitchProduct,
}: {
  product: Product
  onSwitchProduct: (p: Product) => void
}) {
  const [screen, setScreen] = useState<MitoScreen>('overview')
  const [selectedVariantId, setSelectedVariantId] = useState('')

  const [samples, setSamples] = useState<SampleInfo[]>([])
  const [sample, setSample] = useState('')
  const [meta, setMeta] = useState<MitoMeta | null>(null)

  useEffect(() => {
    mitoApi.samples().then((s) => {
      setSamples(s.samples)
      setSample(s.default)
    })
  }, [])

  useEffect(() => {
    if (!sample) return
    mitoApi.meta(sample).then(setMeta)
  }, [sample])

  const onSelectSample = (s: string) => {
    // Variant ids are sample-specific, so drop the selection and return to the
    // overview when switching samples.
    setSample(s)
    setSelectedVariantId('')
    setScreen('overview')
  }

  const goToVariant = (id: string) => {
    setSelectedVariantId(id)
    setScreen('detail')
  }

  return (
    <div style={{ display: 'flex', width: '100%', minHeight: '100vh', background: colors.pageBg, color: colors.textPrimary }}>
      <MitoSidebar
        product={product}
        onSwitchProduct={onSwitchProduct}
        screen={screen}
        onNavigate={setScreen}
        samples={samples}
        sample={sample}
        onSelectSample={onSelectSample}
        caller={meta?.caller ?? ''}
        instrumentLine={meta?.instrument_line ?? ''}
        tissue={meta?.tissue}
      />

      <div style={{ flex: 1, display: 'flex', flexDirection: 'column', minWidth: 0 }}>
        <MitoTopBar
          title={TITLES[screen]}
          sample={meta?.sample ?? sample}
          runDate={meta?.run_date ?? '—'}
          reference={meta?.reference ?? ''}
        />

        <div style={{ flex: 1, overflowY: 'auto' }}>
          {sample && screen === 'overview' && <MitoOverviewScreen sample={sample} onSelectVariant={goToVariant} />}
          {sample && screen === 'detail' && <MitoVariantDetailScreen sample={sample} variantId={selectedVariantId} />}
          {sample && screen === 'qc' && <MitoQCScreen sample={sample} />}
        </div>
      </div>
    </div>
  )
}
