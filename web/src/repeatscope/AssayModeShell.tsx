import { useEffect, useState } from 'react'
import { Sidebar } from '../components/Sidebar'
import { TopBar } from '../components/TopBar'
import type { Product } from '../components/ProductSwitcher'
import { colors } from '../theme'
import { api, type LocusSummary, type ModesResponse } from '../api'
import { OverviewScreen } from '../screens/OverviewScreen'
import { LocusDetailScreen } from '../screens/LocusDetailScreen'
import { GeneFindingScreen } from '../screens/GeneFindingScreen'
import { SampleQCScreen } from '../screens/SampleQCScreen'
import { RunQCScreen } from '../screens/RunQCScreen'
import { RepeatCatalogScreen } from '../screens/RepeatCatalogScreen'

export type Screen = 'overview' | 'locus' | 'gene-finding' | 'qc' | 'runqc' | 'catalog'

const TITLES: Record<Screen, string> = {
  overview: 'Overview',
  locus: 'Locus Detail',
  'gene-finding': 'Gene Finding',
  qc: 'Sample QC Dashboard',
  runqc: 'Run QC Dashboard',
  catalog: 'Repeat Catalog',
}

export type RepeatScopeMode = 'genome_wide' | 'puretarget_v2' | 'carrier'

// The shell shared by RepeatScope's three separate assay-mode apps (see
// GenomeWideApp / PureTargetV2App / CarrierApp) -- each mounts this with its
// own fixed `mode` and `screens`, so switching modes in the Suite switcher
// unmounts one app and mounts another rather than one app internally
// reassigning its mode. That means every piece of per-mode state below
// (`screen`, `sample`, etc.) is naturally fresh per mode via useState's
// initializer -- no explicit reset-on-mode-change effect needed.
export function AssayModeShell({
  mode,
  screens,
  product,
  onSwitchProduct,
}: {
  mode: RepeatScopeMode
  screens: Screen[]
  product: Product
  onSwitchProduct: (p: Product) => void
}) {
  const [screen, setScreen] = useState<Screen>('overview')
  const [selectedLocusId, setSelectedLocusId] = useState<string>('HD_HTT')
  const [selectedGene, setSelectedGene] = useState<string>('')

  const [modesInfo, setModesInfo] = useState<ModesResponse | null>(null)
  const [samples, setSamples] = useState<string[]>([])
  const [sample, setSample] = useState('')
  const [caller, setCaller] = useState('')
  const [sex, setSex] = useState<string | null>(null)
  const [loci, setLoci] = useState<LocusSummary[]>([])

  useEffect(() => {
    api.modes().then(setModesInfo)
  }, [])

  useEffect(() => {
    api.samples(mode).then((s) => {
      setSamples(s.samples)
      setSample(s.default)
    })
  }, [mode])

  useEffect(() => {
    if (!sample) return
    api.meta(mode, sample).then((m) => {
      setCaller(m.caller)
      setSex(m.sex)
    })
  }, [mode, sample])

  // Independent of the OverviewScreen's own overview fetch: the TopBar search
  // box lives above the screen switch below, so it needs the sample's loci
  // list regardless of which screen is currently showing.
  useEffect(() => {
    if (!sample) return
    api.overview(mode, sample).then((o) => setLoci(o.loci))
  }, [mode, sample])

  const goToLocus = (id: string) => {
    setSelectedLocusId(id)
    setScreen('locus')
  }

  const goToGeneFinding = (gene: string) => {
    setSelectedGene(gene)
    setScreen('gene-finding')
  }

  const activeModeLabel = modesInfo?.modes.find((m) => m.key === mode)?.label ?? ''

  return (
    <div className="app-shell" style={{ display: 'flex', width: '100%', minHeight: '100vh', background: colors.pageBg, color: colors.textPrimary }}>
      <Sidebar
        product={product}
        onSwitchProduct={onSwitchProduct}
        screen={screen}
        onNavigate={setScreen}
        availableScreens={screens}
        sample={sample}
        samples={samples}
        onSelectSample={setSample}
        caller={caller}
      />

      <div className="app-main" style={{ flex: 1, display: 'flex', flexDirection: 'column', minWidth: 0 }}>
        <TopBar
          title={TITLES[screen]}
          sample={sample}
          modeLabel={activeModeLabel}
          sex={sex}
          loci={loci}
          onSelectLocus={goToLocus}
          onSelectGeneFinding={goToGeneFinding}
        />

        <div className="app-screenbody" style={{ flex: 1, overflowY: 'auto' }}>
          {!sample && (
            <div style={{ padding: 28, color: colors.textSecondary }}>
              No samples available for {activeModeLabel} yet.
            </div>
          )}
          {sample && screen === 'overview' && (
            <OverviewScreen mode={mode} sample={sample} onSelectLocus={goToLocus} onSelectGeneFinding={goToGeneFinding} />
          )}
          {sample && screen === 'locus' && <LocusDetailScreen mode={mode} sample={sample} locusId={selectedLocusId} />}
          {sample && screen === 'gene-finding' && <GeneFindingScreen mode={mode} sample={sample} gene={selectedGene} />}
          {sample && screen === 'qc' && <SampleQCScreen mode={mode} sample={sample} />}
          {sample && screen === 'runqc' && <RunQCScreen mode={mode} sample={sample} />}
          {sample && screen === 'catalog' && <RepeatCatalogScreen mode={mode} sample={sample} />}
        </div>
      </div>
    </div>
  )
}
