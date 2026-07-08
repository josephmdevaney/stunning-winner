import { useEffect, useState } from 'react'
import { MitoSidebar, type MitoScreen } from './MitoSidebar'
import { MitoTopBar } from './MitoTopBar'
import type { Product } from '../components/ProductSwitcher'
import { colors } from './theme'
import { mitoApi, type MitoMeta } from './api'
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
  const [selectedVariantId, setSelectedVariantId] = useState('m3243AG')
  // Configurable coverage-chart Y-axis max (QC dashboard). A per-report setting.
  const [coverageYMax, setCoverageYMax] = useState(6000)
  const [meta, setMeta] = useState<MitoMeta | null>(null)

  useEffect(() => {
    mitoApi.meta().then(setMeta)
  }, [])

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
        sample={meta?.sample ?? 'MT-0417-P2'}
        instrumentLine={meta?.instrument_line ?? 'PacBio HiFi · mtDNA amplicon'}
      />

      <div style={{ flex: 1, display: 'flex', flexDirection: 'column', minWidth: 0 }}>
        <MitoTopBar
          title={TITLES[screen]}
          sample={meta?.sample ?? 'MT-0417-P2'}
          runDate={meta?.run_date ?? '—'}
          reference={meta?.reference ?? 'rCRS (NC_012920.1)'}
        />

        <div style={{ flex: 1, overflowY: 'auto' }}>
          {screen === 'overview' && <MitoOverviewScreen onSelectVariant={goToVariant} />}
          {screen === 'detail' && <MitoVariantDetailScreen variantId={selectedVariantId} />}
          {screen === 'qc' && <MitoQCScreen coverageYMax={coverageYMax} onChangeCoverageYMax={setCoverageYMax} />}
        </div>
      </div>
    </div>
  )
}
