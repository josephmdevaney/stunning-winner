import { useEffect, useState } from 'react'
import type { Product } from '../components/ProductSwitcher'
import { colors } from '../theme'
import { epigeneticsApi, type SampleInfo, type ThresholdOverrides, type Window } from './api'
import { EpigeneticsSidebar, type EpigeneticsScreen } from './EpigeneticsSidebar'
import { EpigeneticsTopBar } from './EpigeneticsTopBar'
import { ImprintingScreen } from './screens/ImprintingScreen'
import { MethodsScreen } from './screens/MethodsScreen'
import { EpisignatureScreen } from './screens/EpisignatureScreen'
import { EpisignaturesScreen } from './screens/EpisignaturesScreen'
import { RegionScreen } from './screens/RegionScreen'
import { useEpisignatures, useImprinting } from './useEpigeneticsData'
import { Message } from './ui'

const TITLES: Record<EpigeneticsScreen, string> = {
  overview: 'Imprinting',
  region: 'Region Detail',
  episignatures: 'Episignatures',
  signature: 'Signature Detail',
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
  const [modelId, setModelId] = useState<string | null>(null)

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
  const wantsEpisignatures = screen === 'episignatures' || screen === 'signature'
  const episignatures = useEpisignatures(sample, wantsEpisignatures)

  const openRegion = (id: string) => {
    setRegionId(id)
    setScreen('region')
  }

  const openModel = (id: string) => {
    setModelId(id)
    setScreen('signature')
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
        signatureSelected={modelId != null}
        samples={samples}
        sample={sample}
        onSelectSample={(s) => {
          setSample(s)
          setRegionId(null)
          setModelId(null)
          if (screen === 'region') setScreen('overview')
          if (screen === 'signature') setScreen('episignatures')
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
          {sample && screen === 'episignatures' && <EpisignaturesScreen data={episignatures.data} error={episignatures.error} onOpenModel={openModel} />}
          {sample && screen === 'signature' && !episignatures.data && !episignatures.error && <Message>Scoring episignatures…</Message>}
          {screen === 'signature' && episignatures.error && !episignatures.data && <Message tone="error">Episignatures could not be scored: {episignatures.error}</Message>}
          {screen === 'signature' && episignatures.data?.installed && (
            <EpisignatureScreen data={episignatures.data} modelId={modelId ?? episignatures.data.models[0].id} onSelectModel={setModelId} />
          )}
          {overview && screen === 'methods' && (
            <MethodsScreen overview={overview} window={window} overrides={overrides} onChangeWindow={setWindow} onChangeOverrides={setOverrides} />
          )}
        </div>
      </div>
    </div>
  )
}
