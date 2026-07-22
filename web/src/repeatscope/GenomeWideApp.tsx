import type { Product } from '../components/ProductSwitcher'
import { AssayModeShell, type Screen } from './AssayModeShell'

const SCREENS: Screen[] = ['overview', 'locus', 'roh', 'qc', 'runqc', 'catalog']

// Genome-wide is the only assay mode with Runs of Homozygosity (bcftools roh
// is keyed by sample, independent of assay mode) and the Repeat Catalog (the
// full STRchive reference, not a targeted panel) -- PureTarget V2 and Carrier
// are targeted panels where neither applies (see PureTargetV2App/CarrierApp).
export default function GenomeWideApp({
  product,
  onSwitchProduct,
}: {
  product: Product
  onSwitchProduct: (p: Product) => void
}) {
  return <AssayModeShell mode="genome_wide" screens={SCREENS} product={product} onSwitchProduct={onSwitchProduct} />
}
