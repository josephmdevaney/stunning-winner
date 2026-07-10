import { useState } from 'react'
import RepeatScopeApp from './RepeatScopeApp'
import MitoApp from './mito/MitoApp'
import type { Product } from './components/ProductSwitcher'

// The suite ships two sibling products behind one deployable, sharing the same
// visual system: RepeatScope (nuclear tandem-repeat expansions) and MitoLR
// (mitochondrial variant review). A sidebar switcher moves between them.
export default function App() {
  const [product, setProduct] = useState<Product>('repeatscope')

  return product === 'repeatscope' ? (
    <RepeatScopeApp product={product} onSwitchProduct={setProduct} />
  ) : (
    <MitoApp product={product} onSwitchProduct={setProduct} />
  )
}
