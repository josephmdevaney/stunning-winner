import { useEffect, useState } from 'react'
import { ParaphaseSidebar, type ParaphaseScreen } from './ParaphaseSidebar'
import { ParaphaseTopBar } from './ParaphaseTopBar'
import type { Product } from '../components/ProductSwitcher'
import { colors } from '../theme'
import { paraphaseApi } from './api'
import { ParaphaseCallsScreen } from './screens/ParaphaseCallsScreen'
import { ParaphaseOverviewScreen } from './screens/ParaphaseOverviewScreen'
import { ParaphaseRegionScreen } from './screens/ParaphaseRegionScreen'
import { ParaphaseRegionsScreen } from './screens/ParaphaseRegionsScreen'
import { ParaphaseVersionsScreen } from './screens/ParaphaseVersionsScreen'

const TITLES: Record<ParaphaseScreen, string> = {
  regions: 'Paraphase Regions',
  calls: 'Region Calls',
  overview: 'Paraphase Annotation Overview',
  region: 'Annotation Detail',
  versions: 'PacBio Tool Versions',
}

export default function ParaphaseApp({
  product,
  onSwitchProduct,
}: {
  product: Product
  onSwitchProduct: (p: Product) => void
}) {
  const [screen, setScreen] = useState<ParaphaseScreen>('regions')
  const [selectedCall, setSelectedCall] = useState('')
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
    setSelectedCall('')
    setSelectedRegion('')
    setScreen(screen === 'versions' ? 'versions' : 'regions')
  }

  const goToCall = (region: string) => {
    setSelectedCall(region)
    setScreen('calls')
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
        callSelected={!!selectedCall}
        regionSelected={!!selectedRegion}
        samples={samples}
        sample={sample}
        onSelectSample={onSelectSample}
      />

      <div style={{ flex: 1, display: 'flex', flexDirection: 'column', minWidth: 0 }}>
        <ParaphaseTopBar
          title={TITLES[screen]}
          subtitle={
            screen === 'versions'
              ? 'Which release of each PacBio tool made your data, against the newest on GitHub'
              : `Sample ${sample} · ${screen === 'regions' || screen === 'calls' ? "Paraphase's own JSON" : 'interpret_paraphase.py annotation'} · GRCh38`
          }
        />

        <div style={{ flex: 1, overflowY: 'auto' }}>
          {!sample && screen !== 'versions' && (
            <div style={{ padding: 28, color: colors.textSecondary, fontSize: 13.5 }}>
              No samples with a Paraphase JSON on file yet.
            </div>
          )}
          {screen === 'versions' && <ParaphaseVersionsScreen />}
          {sample && screen === 'regions' && <ParaphaseRegionsScreen sample={sample} onSelectRegion={goToCall} />}
          {sample && screen === 'calls' && selectedCall && <ParaphaseCallsScreen sample={sample} region={selectedCall} />}
          {sample && screen === 'overview' && <ParaphaseOverviewScreen sample={sample} onSelectRegion={goToRegion} />}
          {sample && screen === 'region' && selectedRegion && <ParaphaseRegionScreen sample={sample} region={selectedRegion} />}
        </div>
      </div>
    </div>
  )
}
