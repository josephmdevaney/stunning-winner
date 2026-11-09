import { useEffect, useState } from 'react'
import { KivviSidebar, type KivviScreen } from './KivviSidebar'
import { KivviTopBar } from './KivviTopBar'
import type { Product } from '../components/ProductSwitcher'
import { colors } from './theme'
import { kivviApi, type KivviMeta, type SampleInfo } from './api'
import { KivviOverviewScreen } from './screens/KivviOverviewScreen'
import { D4Z4Screen } from './screens/D4Z4Screen'
import { FSHD2Screen } from './screens/FSHD2Screen'
import { Kiv2Screen } from './screens/Kiv2Screen'
import { KivviQCScreen } from './screens/KivviQCScreen'

const TITLES: Record<KivviScreen, string> = {
  overview: 'Overview',
  d4z4: 'D4Z4 Repeat (FSHD)',
  fshd2: 'FSHD2 Evidence',
  kiv2: 'LPA KIV2 Repeat',
  qc: 'Sample QC',
}

export default function KivviApp({
  product,
  onSwitchProduct,
}: {
  product: Product
  onSwitchProduct: (p: Product) => void
}) {
  const [screen, setScreen] = useState<KivviScreen>('overview')
  const [samples, setSamples] = useState<SampleInfo[]>([])
  const [sample, setSample] = useState('')
  const [meta, setMeta] = useState<KivviMeta | null>(null)
  const [loaded, setLoaded] = useState(false)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    kivviApi
      .samples()
      .then((s) => {
        setSamples(s.samples)
        setSample(s.default)
        setLoaded(true)
      })
      .catch((e) => setError(String(e)))
  }, [])

  useEffect(() => {
    if (!sample) return
    kivviApi.meta(sample).then(setMeta).catch((e) => setError(String(e)))
  }, [sample])

  const subtitle = meta ? `${meta.reference}${meta.input_note ? ` · ${meta.input_note}` : ''}` : 'GRCh38'

  return (
    <div style={{ display: 'flex', width: '100%', minHeight: '100vh', background: colors.pageBg, color: colors.textPrimary }}>
      <KivviSidebar
        product={product}
        onSwitchProduct={onSwitchProduct}
        screen={screen}
        onNavigate={setScreen}
        samples={samples}
        sample={sample}
        onSelectSample={setSample}
        kivviVersion={meta?.kivvi_version}
        reference={meta?.reference}
      />

      <div style={{ flex: 1, display: 'flex', flexDirection: 'column', minWidth: 0 }}>
        <KivviTopBar title={TITLES[screen]} sample={sample} subtitle={subtitle} />

        <div style={{ flex: 1, overflowY: 'auto' }}>
          {error && <div style={{ padding: 28, color: colors.statusPathogenic }}>Failed to load: {error}</div>}
          {loaded && !sample && <NoSamples />}
          {sample && screen === 'overview' && <KivviOverviewScreen sample={sample} onNavigate={setScreen} />}
          {sample && screen === 'd4z4' && <D4Z4Screen sample={sample} />}
          {sample && screen === 'fshd2' && <FSHD2Screen sample={sample} />}
          {sample && screen === 'kiv2' && <Kiv2Screen sample={sample} />}
          {sample && screen === 'qc' && <KivviQCScreen sample={sample} />}
        </div>
      </div>
    </div>
  )
}

function NoSamples() {
  return (
    <div style={{ padding: 28, maxWidth: 720, lineHeight: 1.6 }}>
      <div style={{ fontSize: 16, fontWeight: 600 }}>No Kivvi results found</div>
      <div style={{ color: colors.textSecondary, marginTop: 8, fontSize: 13.5 }}>
        Run Kivvi on a WGS BAM (<code>kivvi -b sample.bam -o &lt;root&gt;/&lt;sample&gt;/d4z4 -p &lt;sample&gt; d4z4</code>, likewise for{' '}
        <code>kiv2</code>) and point <code>KIVVI_RESULTS_ROOT</code> at the <code>&lt;root&gt;</code> directory. Each sample appears here as soon as{' '}
        <code>&lt;root&gt;/&lt;sample&gt;/&lt;target&gt;/&lt;sample&gt;.kivvi.&lt;target&gt;.json</code> exists.
      </div>
    </div>
  )
}
