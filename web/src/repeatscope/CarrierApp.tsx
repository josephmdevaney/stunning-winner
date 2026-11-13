import type { Product } from '../components/ProductSwitcher'
import { AssayModeShell, type Screen } from './AssayModeShell'

const SCREENS: Screen[] = ['overview', 'locus', 'qc', 'runqc']

export default function CarrierApp({
  product,
  onSwitchProduct,
}: {
  product: Product
  onSwitchProduct: (p: Product) => void
}) {
  return <AssayModeShell mode="carrier" screens={SCREENS} product={product} onSwitchProduct={onSwitchProduct} />
}
