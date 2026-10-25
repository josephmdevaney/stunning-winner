import { useEffect, useState } from 'react'
import type { Product } from '../components/ProductSwitcher'
import { colors } from '../theme'
import { epigeneticsApi, type SampleInfo, type ThresholdOverrides, type Window } from './api'
import { EpigeneticsSidebar, type EpigeneticsScreen } from './EpigeneticsSidebar'
import { EpigeneticsTopBar } from './EpigeneticsTopBar'
import { ImprintingScreen } from './screens/ImprintingScreen'
import { MethodsScreen } from './screens/MethodsScreen'
import { RegionScreen } from './screens/RegionScreen'
import { useImprinting } from './useEpigeneticsData'
import { Message } from './ui'

const TITLES: Record<EpigeneticsScreen, string> = {
  overview: 'Imprinting',
  region: 'Region Detail',
  methods: 'Methods',
}

export default function EpigeneticsApp({ product, onSwitchProduct }: { product: Product; onSwitchProduct: (p: Product) => void }) {
  const [screen, setScreen] = useState<EpigeneticsScreen>('overview')
  const [samples, setSamples] = useState<SampleInfo[]>([])
  const [sample, setSample] = useState('')
  const [loaded, setLoaded] = useState(false)
  const [samplesError, setSamplesError] = useState<string | null>(null)
  const [window, setWindow] = useState<Window>('core')
  const [overrides, setOverrides] = useState<ThresholdOverrides>({})
  const [regionId, setRegionId] = useState<string | null>(null)

  useEffect(() => {
    epigeneticsApi
      .samples()
      .then((s) => {
        setSamples(s.samples)
        setSample(s.default)
        setLoaded(true)
      })
      .catch((e) => setSamplesError(String(e.message ?? e)))
  }, [])

  const { overview, error, loading } = useImprinting(sample, window, overrides)

  const openRegion = (id: string) => {
    setRegionId(id)
    setScreen('region')
  }

  const sampleInfo = samples.find((s) => s.key === sample)
  const subtitle = overview
    ? `${overview.provenance.genome_build} · ${overview.provenance.pb_cpg_tools_version ? `pb-CpG-tools ${overview.provenance.pb_cpg_tools_version}` : 'methylation pileup'}${loading ? ' · updating…' : ''}`
    : 'GRCh38'

  return (
    <div style={{ display: 'flex', width: '100%', minHeight: '100vh', background: colors.pageBg, color: colors.textPrimary }}>
      <EpigeneticsSidebar
        product={product}
        onSwitchProduct={onSwitchProduct}
        screen={screen}
        onNavigate={setScreen}
        regionSelected={regionId != null}
        samples={samples}
        sample={sample}
        onSelectSample={(s) => {
          setSample(s)
          setRegionId(null)
          if (screen === 'region') setScreen('overview')
        }}
        footer={sampleInfo ? sampleInfo.tracks.join(' · ') : 'GRCh38'}
      />

      <div style={{ flex: 1, display: 'flex', flexDirection: 'column', minWidth: 0 }}>
        <EpigeneticsTopBar title={TITLES[screen]} subtitle={subtitle} />

        <div style={{ flex: 1, overflowY: 'auto' }}>
          {samplesError && <Message tone="error">Failed to load samples: {samplesError}</Message>}
          {loaded && !sample && <Message>No genomes are registered yet. Add a sample's pb-CpG-tools output to backend/config/epigenetics.local.json (see epigenetics.local.json.example).</Message>}
          {sample && error && !overview && <Message tone="error">No methylation analysis could be loaded for "{sample}": {error}</Message>}
          {sample && !error && !overview && <Message>Reading methylation…</Message>}
          {overview && screen === 'overview' && <ImprintingScreen overview={overview} onOpenRegion={openRegion} />}
          {overview && screen === 'region' && (
            <RegionScreen sample={sample} window={window} overview={overview} regionId={regionId ?? overview.regions[0].id} onSelectRegion={setRegionId} />
          )}
          {overview && screen === 'methods' && (
            <MethodsScreen overview={overview} window={window} overrides={overrides} onChangeWindow={setWindow} onChangeOverrides={setOverrides} />
          )}
        </div>
      </div>
    </div>
  )
}
