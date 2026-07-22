import { useEffect, useState } from 'react'
import { Sidebar } from './components/Sidebar'
import { TopBar } from './components/TopBar'
import type { Product } from './components/ProductSwitcher'
import { colors } from './theme'
import { api, type ModesResponse } from './api'
import { OverviewScreen } from './screens/OverviewScreen'
import { LocusDetailScreen } from './screens/LocusDetailScreen'
import { GeneFindingScreen } from './screens/GeneFindingScreen'
import { SampleQCScreen } from './screens/SampleQCScreen'
import { RunQCScreen } from './screens/RunQCScreen'
import { RepeatCatalogScreen } from './screens/RepeatCatalogScreen'
import { ROHScreen } from './screens/ROHScreen'

export type Screen = 'overview' | 'locus' | 'gene-finding' | 'qc' | 'runqc' | 'catalog' | 'roh'

const TITLES: Record<Screen, string> = {
  overview: 'Overview',
  locus: 'Locus Detail',
  'gene-finding': 'Gene Finding',
  qc: 'Sample QC Dashboard',
  runqc: 'Run QC Dashboard',
  catalog: 'Repeat Catalog',
  roh: 'Runs of Homozygosity',
}

// RepeatScope's three assay modes (Genome-wide, PureTarget V2, Carrier) are
// each a top-level "Suite" section (see ProductSwitcher) rather than a nested
// picker, so this app only ever renders for one of those three `product`
// values -- App.tsx's ternary narrows out 'mitolr' before reaching here.
type RepeatScopeProduct = Exclude<Product, 'mitolr'>

// Runs of Homozygosity (bcftools roh, keyed by sample independent of assay
// mode -- see ROHScreen) and the Repeat Catalog (the full STRchive reference,
// not a targeted panel) are only meaningful for Genome-wide; PureTarget V2 and
// Carrier are targeted panels where neither applies.
const SCREENS_BY_PRODUCT: Record<RepeatScopeProduct, Screen[]> = {
  genome_wide: ['overview', 'locus', 'roh', 'qc', 'runqc', 'catalog'],
  puretarget_v2: ['overview', 'locus', 'qc', 'runqc'],
  carrier: ['overview', 'locus', 'qc', 'runqc'],
}

export default function RepeatScopeApp({
  product,
  onSwitchProduct,
}: {
  product: RepeatScopeProduct
  onSwitchProduct: (p: Product) => void
}) {
  const [screen, setScreen] = useState<Screen>('overview')
  const [selectedLocusId, setSelectedLocusId] = useState<string>('HD_HTT')
  const [selectedGene, setSelectedGene] = useState<string>('')

  const [modesInfo, setModesInfo] = useState<ModesResponse | null>(null)
  // `mode` mirrors `product` (they share the same key -- genome_wide /
  // puretarget_v2 / carrier) since the backend's mode keys ARE the Suite
  // section identifiers; kept as separate state only because the samples/meta
  // effects below key off it the same way they always have.
  const [mode, setMode] = useState(product)
  const [samples, setSamples] = useState<string[]>([])
  const [sample, setSample] = useState('')
  const [caller, setCaller] = useState('')
  const [sex, setSex] = useState<string | null>(null)

  useEffect(() => {
    api.modes().then(setModesInfo)
  }, [])

  // Clearing `sample`/`selectedLocusId`/`selectedGene` here (rather than in the
  // `[mode]` effect below) matters: React batches state updates from the same
  // handler into one render, so the still-mounted screen below never sees a (new
  // mode, stale sample) pair -- e.g. requesting Carrier's overview with
  // Genome-wide's sample, which 404s until the effect below catches up a render
  // later.
  useEffect(() => {
    setSample('')
    setSelectedLocusId('')
    setSelectedGene('')
    setMode(product)
    setScreen((prev) => (SCREENS_BY_PRODUCT[product].includes(prev) ? prev : 'overview'))
  }, [product])

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
    api.meta(mode, sample).then((m) => {
      setCaller(m.caller)
      setSex(m.sex)
    })
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
        availableScreens={SCREENS_BY_PRODUCT[product]}
        sample={sample}
        samples={samples}
        onSelectSample={setSample}
        caller={caller}
      />

      <div className="app-main" style={{ flex: 1, display: 'flex', flexDirection: 'column', minWidth: 0 }}>
        <TopBar title={TITLES[screen]} sample={sample} modeLabel={activeModeLabel} sex={sex} />

        <div className="app-screenbody" style={{ flex: 1, overflowY: 'auto' }}>
          {!sample && mode && (
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
          {sample && screen === 'roh' && <ROHScreen sample={sample} />}
        </div>
      </div>
    </div>
  )
}
