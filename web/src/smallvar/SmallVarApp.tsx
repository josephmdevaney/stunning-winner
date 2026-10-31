import { useEffect, useState } from 'react'
import { SmallVarSidebar, type SmallVarScreen } from './SmallVarSidebar'
import { SmallVarTopBar } from './SmallVarTopBar'
import type { Product } from '../components/ProductSwitcher'
import { colors, fmtInt, mono } from './theme'
import { smallvarApi, DEFAULT_FILTERS, type BuildStatus, type FilterState, type Meta, type SampleInfo } from './api'
import { Button } from './components/ui'
import { cardStyle } from './components/styles'
import { OverviewScreen } from './screens/OverviewScreen'
import { VariantsScreen } from './screens/VariantsScreen'

const TITLES: Record<SmallVarScreen, string> = { overview: 'Overview', variants: 'Variants' }
const POLL_MS = 1500

export default function SmallVarApp({ product, onSwitchProduct }: { product: Product; onSwitchProduct: (p: Product) => void }) {
  const [screen, setScreen] = useState<SmallVarScreen>('overview')
  const [samples, setSamples] = useState<SampleInfo[]>([])
  const [sample, setSample] = useState('')
  const [loaded, setLoaded] = useState(false)
  const [meta, setMeta] = useState<Meta | null>(null)
  const [status, setStatus] = useState<BuildStatus | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [retry, setRetry] = useState(0)
  // A funnel step on the overview opens the variant list with that filter applied.
  const [preset, setPreset] = useState<{ filters: FilterState; nonce: number } | null>(null)

  useEffect(() => {
    smallvarApi
      .samples()
      .then((s) => {
        setSamples(s.samples)
        setSample(s.default)
        setLoaded(true)
      })
      .catch((e) => setError(String(e)))
  }, [])

  // Asking for status starts the index build; keep asking until it is ready.
  useEffect(() => {
    if (!sample) return
    let cancelled = false
    let timer: ReturnType<typeof setTimeout> | undefined
    setStatus(null)
    setMeta(null)
    const tick = () => {
      smallvarApi
        .status(sample)
        .then((s) => {
          if (cancelled) return
          setStatus(s)
          if (s.state === 'ready') smallvarApi.meta(sample).then((m) => !cancelled && setMeta(m))
          else if (s.state !== 'error') timer = setTimeout(tick, POLL_MS)
        })
        .catch((e) => !cancelled && setError(String(e)))
    }
    tick()
    return () => {
      cancelled = true
      if (timer) clearTimeout(timer)
    }
  }, [sample, retry])

  const openVariants = (filters: Partial<FilterState>) => {
    setPreset({ filters: { ...DEFAULT_FILTERS, ...filters }, nonce: Date.now() })
    setScreen('variants')
  }

  const ready = status?.state === 'ready' && meta
  const sampleLabel = meta?.label ?? samples.find((s) => s.key === sample)?.label ?? sample
  const subtitle = ready
    ? `${sampleLabel} · ${meta.reference} · ${fmtInt(meta.n_variants)} calls${meta.caller ? ` · ${meta.caller}` : ''}`
    : `${sampleLabel || '—'} · GRCh38`

  return (
    <div style={{ display: 'flex', width: '100%', minHeight: '100vh', background: colors.pageBg, color: colors.textPrimary }}>
      <SmallVarSidebar
        product={product} onSwitchProduct={onSwitchProduct} screen={screen} onNavigate={setScreen}
        samples={samples} sample={sample} onSelectSample={setSample} reference={meta?.reference}
      />
      <div style={{ flex: 1, display: 'flex', flexDirection: 'column', minWidth: 0 }}>
        <SmallVarTopBar title={TITLES[screen]} subtitle={subtitle} />
        <div style={{ flex: 1, overflowY: 'auto' }}>
          {error && <div style={{ padding: 28, color: colors.statusPathogenic }}>Failed to load: {error}</div>}
          {loaded && !sample && <NoSamples />}
          {sample && status && status.state !== 'ready' && <Building status={status} onRetry={() => setRetry((n) => n + 1)} />}
          {sample && ready && screen === 'overview' && <OverviewScreen sample={sample} onOpenVariants={openVariants} />}
          {sample && ready && screen === 'variants' && (
            <VariantsScreen key={`${sample}-${preset?.nonce ?? 0}`} sample={sample} initial={preset?.filters ?? DEFAULT_FILTERS} />
          )}
        </div>
      </div>
    </div>
  )
}

function Building({ status, onRetry }: { status: BuildStatus; onRetry: () => void }) {
  if (status.state === 'error') {
    return (
      <div style={{ padding: 28, maxWidth: 760 }}>
        <div style={{ ...cardStyle, borderColor: colors.statusPathogenic }}>
          <div style={{ fontSize: 15, fontWeight: 600, color: colors.statusPathogenic }}>The variant index could not be built</div>
          <div style={{ fontFamily: mono, fontSize: 12.5, marginTop: 8, wordBreak: 'break-word' }}>{status.error}</div>
          <div style={{ marginTop: 14 }}>
            <Button onClick={onRetry} primary>
              Try again
            </Button>
          </div>
        </div>
      </div>
    )
  }
  const p = status.progress
  return (
    <div style={{ padding: 28, maxWidth: 760 }}>
      <div style={cardStyle}>
        <div style={{ fontSize: 15, fontWeight: 600 }}>Building the variant index</div>
        <div style={{ color: colors.textSecondary, fontSize: 13, marginTop: 6, lineHeight: 1.6 }}>
          The VCF is read once, annotated with ClinVar, gnomAD and gene bodies, and indexed so that every filter after this is instant. A whole
          genome takes a couple of minutes; the index is kept and reused until the VCF or an annotation resource changes.
        </div>
        <div style={{ fontFamily: mono, fontSize: 12.5, marginTop: 14, color: colors.textBody }}>
          {p ? `${p.phase} · ${fmtInt(p.records)} records · ${p.seconds}s` : 'starting…'}
        </div>
      </div>
    </div>
  )
}

function NoSamples() {
  return (
    <div style={{ padding: 28, maxWidth: 720, lineHeight: 1.6 }}>
      <div style={{ fontSize: 16, fontWeight: 600 }}>No small-variant VCFs found</div>
      <div style={{ color: colors.textSecondary, marginTop: 8, fontSize: 13.5 }}>
        Put a single-sample small-variant VCF (DeepVariant output, bgzipped or plain) at <code>&lt;root&gt;/&lt;sample&gt;/&lt;sample&gt;.vcf.gz</code> and point{' '}
        <code>SMALLVAR_RESULTS_ROOT</code> at <code>&lt;root&gt;</code> (default <code>~/smallvar_data/samples</code>). An optional <code>sample.json</code> beside
        it can give a label, <code>sex</code> and <code>sex_source</code>, and <code>caller</code>.
      </div>
    </div>
  )
}
