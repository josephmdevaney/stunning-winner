import type { Product } from '../components/ProductSwitcher'
import { AssayModeShell, type Screen } from './AssayModeShell'

const SCREENS: Screen[] = ['overview', 'locus', 'qc', 'runqc']

export default function PureTargetV2App({
  product,
  onSwitchProduct,
}: {
  product: Product
  onSwitchProduct: (p: Product) => void
}) {
  return <AssayModeShell mode="puretarget_v2" screens={SCREENS} product={product} onSwitchProduct={onSwitchProduct} />
}
