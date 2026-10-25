import type { Product } from '../components/ProductSwitcher'
import { AssayModeShell, type Screen } from './AssayModeShell'

const SCREENS: Screen[] = ['overview', 'locus', 'qc', 'runqc', 'catalog']

// Genome-wide is the only assay mode with the Repeat Catalog (the full STRchive
// reference, not a targeted panel) -- PureTarget V2 and Carrier are targeted
// panels where it doesn't apply (see PureTargetV2App/CarrierApp). Runs of
// Homozygosity is its own suite (see roh/ROHApp).
export default function GenomeWideApp({
  product,
  onSwitchProduct,
}: {
  product: Product
  onSwitchProduct: (p: Product) => void
}) {
  return <AssayModeShell mode="genome_wide" screens={SCREENS} product={product} onSwitchProduct={onSwitchProduct} />
}
