import { useEffect, useState } from 'react'
import { ROHSidebar, type ROHScreenKey } from './ROHSidebar'
import { ROHTopBar } from './ROHTopBar'
import type { Product } from '../components/ProductSwitcher'
import { colors } from '../theme'
import { rohApi } from './api'
import { useRohData } from './useRohData'
import { ROHOverviewScreen } from './screens/ROHOverviewScreen'
import { ROHChromosomeScreen } from './screens/ROHChromosomeScreen'
import { ROHRegionsScreen } from './screens/ROHRegionsScreen'

const TITLES: Record<ROHScreenKey, string> = {
  overview: 'Runs of Homozygosity',
  chromosome: 'Chromosome Detail',
  regions: 'ROH Regions',
}

const SEX_LABEL: Record<string, string> = { M: 'Male (XY)', F: 'Female (XX)' }

export default function ROHApp({
  product,
  onSwitchProduct,
}: {
  product: Product
  onSwitchProduct: (p: Product) => void
}) {
  const [screen, setScreen] = useState<ROHScreenKey>('overview')
  const [samples, setSamples] = useState<string[]>([])
  const [sample, setSample] = useState('')
  const [loaded, setLoaded] = useState(false)
  const [samplesError, setSamplesError] = useState<string | null>(null)

  useEffect(() => {
    rohApi
      .samples()
      .then((s) => {
        setSamples(s.samples)
        setSample(s.default)
        setLoaded(true)
      })
      .catch((e) => setSamplesError(String(e)))
  }, [])

  const { overview, allRegions, error } = useRohData(sample)

  const sexLabel = overview?.sex_inference.sex ? SEX_LABEL[overview.sex_inference.sex] : null
  const subtitle = overview ? `${overview.genome_build} · ${overview.caller}${sexLabel ? ` · ${sexLabel}` : ''}` : 'GRCh38'

  return (
    <div style={{ display: 'flex', width: '100%', minHeight: '100vh', background: colors.pageBg, color: colors.textPrimary }}>
      <ROHSidebar
        product={product}
        onSwitchProduct={onSwitchProduct}
        screen={screen}
        onNavigate={setScreen}
        samples={samples}
        sample={sample}
        onSelectSample={setSample}
        genomeBuild={overview?.genome_build}
        caller={overview?.caller}
      />

      <div style={{ flex: 1, display: 'flex', flexDirection: 'column', minWidth: 0 }}>
        <ROHTopBar title={TITLES[screen]} sample={sample} subtitle={subtitle} />

        <div style={{ flex: 1, overflowY: 'auto' }}>
          {samplesError && <Message tone="error">Failed to load samples: {samplesError}</Message>}
          {loaded && !sample && <Message>No ROH results are registered yet. Add a sample's bcftools roh output to the ROH sample registry.</Message>}
          {sample && error && <Message tone="error">No ROH analysis could be loaded for sample "{sample}": {error}</Message>}
          {sample && !error && !overview && <Message>Loading ROH analysis…</Message>}
          {sample && overview && screen === 'overview' && <ROHOverviewScreen overview={overview} allRegions={allRegions} />}
          {sample && overview && screen === 'chromosome' && <ROHChromosomeScreen overview={overview} allRegions={allRegions} />}
          {sample && overview && screen === 'regions' && <ROHRegionsScreen key={sample} sample={sample} />}
        </div>
      </div>
    </div>
  )
}

function Message({ children, tone }: { children: React.ReactNode; tone?: 'error' }) {
  return (
    <div style={{ padding: 28 }}>
      <div
        style={{
          background: 'white',
          border: `1px solid ${colors.border}`,
          borderRadius: 12,
          padding: '60px 24px',
          textAlign: 'center',
          color: tone === 'error' ? colors.statusPathogenic : colors.textSecondary,
          fontSize: 14,
        }}
      >
        {children}
      </div>
    </div>
  )
}
