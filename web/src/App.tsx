import { useState } from 'react'
import RepeatScopeApp from './RepeatScopeApp'
import MitoApp from './mito/MitoApp'
import type { Product } from './components/ProductSwitcher'

// The suite ships four sibling sections behind one deployable, sharing the
// same visual system: RepeatScope's three assay modes (genome-wide,
// PureTarget V2, Carrier) and MitoLR (mitochondrial variant review). A sidebar
// switcher moves between all four; RepeatScope's three modes share one app
// (RepeatScopeApp), which reacts to `product` to pick its active assay mode.
export default function App() {
  const [product, setProduct] = useState<Product>('genome_wide')

  return product === 'mitolr' ? (
    <MitoApp product={product} onSwitchProduct={setProduct} />
  ) : (
    <RepeatScopeApp product={product} onSwitchProduct={setProduct} />
  )
}
