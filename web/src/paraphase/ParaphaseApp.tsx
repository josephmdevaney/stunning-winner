import { useEffect, useState } from 'react'
import { ParaphaseSidebar, type ParaphaseScreen } from './ParaphaseSidebar'
import { ParaphaseTopBar } from './ParaphaseTopBar'
import type { Product } from '../components/ProductSwitcher'
import { colors } from '../theme'
import { paraphaseApi } from './api'
import { ParaphaseOverviewScreen } from './screens/ParaphaseOverviewScreen'
import { ParaphaseRegionScreen } from './screens/ParaphaseRegionScreen'

const TITLES: Record<ParaphaseScreen, string> = {
  overview: 'Paraphase Annotation Overview',
  region: 'Region Detail',
}

export default function ParaphaseApp({
  product,
  onSwitchProduct,
}: {
  product: Product
  onSwitchProduct: (p: Product) => void
}) {
  const [screen, setScreen] = useState<ParaphaseScreen>('overview')
  const [selectedRegion, setSelectedRegion] = useState('')

  const [samples, setSamples] = useState<string[]>([])
  const [sample, setSample] = useState('')

  useEffect(() => {
    paraphaseApi.samples().then((s) => {
      setSamples(s.samples)
      setSample(s.samples[0] ?? '')
    })
  }, [])

  const onSelectSample = (s: string) => {
    setSample(s)
    setSelectedRegion('')
    setScreen('overview')
  }

  const goToRegion = (region: string) => {
    setSelectedRegion(region)
    setScreen('region')
  }

  return (
    <div style={{ display: 'flex', width: '100%', minHeight: '100vh', background: colors.pageBg, color: colors.textPrimary }}>
      <ParaphaseSidebar
        product={product}
        onSwitchProduct={onSwitchProduct}
        screen={screen}
        onNavigate={setScreen}
        regionSelected={!!selectedRegion}
        samples={samples}
        sample={sample}
        onSelectSample={onSelectSample}
      />

      <div style={{ flex: 1, display: 'flex', flexDirection: 'column', minWidth: 0 }}>
        <ParaphaseTopBar title={TITLES[screen]} sample={sample} />

        <div style={{ flex: 1, overflowY: 'auto' }}>
          {!sample && (
            <div style={{ padding: 28, color: colors.textSecondary, fontSize: 13.5 }}>
              No samples with a Paraphase JSON on file yet.
            </div>
          )}
          {sample && screen === 'overview' && <ParaphaseOverviewScreen sample={sample} onSelectRegion={goToRegion} />}
          {sample && screen === 'region' && selectedRegion && <ParaphaseRegionScreen sample={sample} region={selectedRegion} />}
        </div>
      </div>
    </div>
  )
}
