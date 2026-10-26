import { useState } from 'react'
import GenomeWideApp from './repeatscope/GenomeWideApp'
import PureTargetV2App from './repeatscope/PureTargetV2App'
import CarrierApp from './repeatscope/CarrierApp'
import MitoApp from './mito/MitoApp'
import ParaphaseApp from './paraphase/ParaphaseApp'
import KivviApp from './kivvi/KivviApp'
import ROHApp from './roh/ROHApp'
import CohortApp from './cohort/CohortApp'
import type { Product } from './components/ProductSwitcher'

// The suite ships seven separate tools behind one deployable, combined into a
// single Suite switcher and sharing one visual system: three RepeatScope
// assay-mode apps (GenomeWideApp, PureTargetV2App, CarrierApp -- each its own
// component, not one app internally reassigning its mode), MitoLR (MitoApp),
// Paraphase Annotation (ParaphaseApp), Kivvi (KivviApp), and Runs of
// Homozygosity (ROHApp). Switching `product` swaps which tool is mounted.
export default function App() {
  const [product, setProduct] = useState<Product>('genome_wide')

  switch (product) {
    case 'mitolr':
      return <MitoApp product={product} onSwitchProduct={setProduct} />
    case 'paraphase':
      return <ParaphaseApp product={product} onSwitchProduct={setProduct} />
    case 'kivvi':
      return <KivviApp product={product} onSwitchProduct={setProduct} />
    case 'roh':
      return <ROHApp product={product} onSwitchProduct={setProduct} />
    case 'cohort':
      return <CohortApp product={product} onSwitchProduct={setProduct} />
    case 'puretarget_v2':
      return <PureTargetV2App product={product} onSwitchProduct={setProduct} />
    case 'carrier':
      return <CarrierApp product={product} onSwitchProduct={setProduct} />
    default:
      return <GenomeWideApp product={product} onSwitchProduct={setProduct} />
  }
}
