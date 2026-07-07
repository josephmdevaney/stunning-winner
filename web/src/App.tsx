import { useEffect, useState } from 'react'
import { Sidebar } from './components/Sidebar'
import { TopBar } from './components/TopBar'
import { colors } from './theme'
import { api, type ModesResponse } from './api'
import { OverviewScreen } from './screens/OverviewScreen'
import { LocusDetailScreen } from './screens/LocusDetailScreen'
import { SampleQCScreen } from './screens/SampleQCScreen'
import { RunQCScreen } from './screens/RunQCScreen'
import { RepeatCatalogScreen } from './screens/RepeatCatalogScreen'

export type Screen = 'overview' | 'locus' | 'qc' | 'runqc' | 'catalog'

const TITLES: Record<Screen, string> = {
  overview: 'Overview',
  locus: 'Locus Detail',
  qc: 'Sample QC Dashboard',
  runqc: 'Run QC Dashboard',
  catalog: 'Repeat Catalog',
}

export default function App() {
  const [screen, setScreen] = useState<Screen>('overview')
  const [selectedLocusId, setSelectedLocusId] = useState<string>('')

  const [modesInfo, setModesInfo] = useState<ModesResponse | null>(null)
  const [mode, setMode] = useState('')
  const [samples, setSamples] = useState<string[]>([])
  const [sample, setSample] = useState('')
  const [caller, setCaller] = useState('')

  useEffect(() => {
    api.modes().then((m) => {
      setModesInfo(m)
      setMode(m.default)
    })
  }, [])

  // Whenever the mode changes, load that mode's sample list and pick its default.
  useEffect(() => {
    if (!mode) return
    api.samples(mode).then((s) => {
      setSamples(s.samples)
      setSample(s.default)
    })
  }, [mode])

  useEffect(() => {
    if (!mode || !sample) return
    api.meta(mode, sample).then((m) => setCaller(m.caller))
  }, [mode, sample])

  const goToLocus = (id: string) => {
    setSelectedLocusId(id)
    setScreen('locus')
  }

  // Clearing `sample`/`selectedLocusId` here (rather than in the `[mode]` effect
  // above) matters: React batches state updates from the same handler into one
  // render, so the still-mounted screen below never sees a (new mode, stale
  // sample) pair -- e.g. requesting Carrier's overview with Genome-wide's sample,
  // which 404s until the effect above catches up a render later.
  const selectMode = (newMode: string) => {
    setSample('')
    setSelectedLocusId('')
    setMode(newMode)
  }

  const activeModeLabel = modesInfo?.modes.find((m) => m.key === mode)?.label ?? ''

  return (
    <div className="app-shell" style={{ display: 'flex', width: '100%', minHeight: '100vh', background: colors.pageBg, color: colors.textPrimary }}>
      <Sidebar
        screen={screen}
        onNavigate={setScreen}
        modes={modesInfo?.modes ?? []}
        mode={mode}
        onSelectMode={selectMode}
        sample={sample}
        samples={samples}
        onSelectSample={setSample}
        caller={caller}
      />

      <div className="app-main" style={{ flex: 1, display: 'flex', flexDirection: 'column', minWidth: 0 }}>
        <TopBar title={TITLES[screen]} sample={sample} modeLabel={activeModeLabel} />

        <div className="app-screenbody" style={{ flex: 1, overflowY: 'auto' }}>
          {!sample && mode && (
            <div style={{ padding: 28, color: colors.textSecondary }}>
              No samples available for {activeModeLabel} yet.
            </div>
          )}
          {sample && screen === 'overview' && <OverviewScreen mode={mode} sample={sample} onSelectLocus={goToLocus} />}
          {sample && screen === 'locus' && <LocusDetailScreen mode={mode} sample={sample} locusId={selectedLocusId} />}
          {sample && screen === 'qc' && <SampleQCScreen mode={mode} sample={sample} />}
          {sample && screen === 'runqc' && <RunQCScreen mode={mode} sample={sample} />}
          {sample && screen === 'catalog' && <RepeatCatalogScreen mode={mode} sample={sample} />}
        </div>
      </div>
    </div>
  )
}
