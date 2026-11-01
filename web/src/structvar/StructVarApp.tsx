import { useEffect, useState } from 'react'
import { StructVarSidebar, type StructVarScreen } from './StructVarSidebar'
import { StructVarTopBar } from './StructVarTopBar'
import type { Product } from '../components/ProductSwitcher'
import { colors, fmtInt, mono } from './theme'
import { structvarApi, DEFAULT_FILTERS, type BuildStatus, type FilterState, type Meta, type SampleInfo } from './api'
import { Button } from './components/ui'
import { cardStyle } from './components/styles'
import { OverviewScreen } from './screens/OverviewScreen'
import { SvsScreen } from './screens/SvsScreen'

const TITLES: Record<StructVarScreen, string> = { overview: 'Overview', svs: 'Structural variants' }
const POLL_MS = 1500

export default function StructVarApp({ product, onSwitchProduct }: { product: Product; onSwitchProduct: (p: Product) => void }) {
  const [screen, setScreen] = useState<StructVarScreen>('overview')
  const [samples, setSamples] = useState<SampleInfo[]>([])
  const [sample, setSample] = useState('')
  const [loaded, setLoaded] = useState(false)
  const [meta, setMeta] = useState<Meta | null>(null)
  const [status, setStatus] = useState<BuildStatus | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [retry, setRetry] = useState(0)
  // A funnel step on the overview opens the SV list with that filter applied.
  const [preset, setPreset] = useState<{ filters: FilterState; nonce: number } | null>(null)

  useEffect(() => {
    structvarApi
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
      structvarApi
        .status(sample)
        .then((s) => {
          if (cancelled) return
          setStatus(s)
          if (s.state === 'ready') structvarApi.meta(sample).then((m) => !cancelled && setMeta(m))
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

  const openSvs = (filters: Partial<FilterState>) => {
    setPreset({ filters: { ...DEFAULT_FILTERS, ...filters }, nonce: Date.now() })
    setScreen('svs')
  }

  const ready = status?.state === 'ready' && meta
  const label = meta?.label ?? samples.find((s) => s.key === sample)?.label ?? sample
  const subtitle = ready
    ? `${label} · ${meta.reference} · ${fmtInt(meta.n_svs)} calls${meta.source ? ` · ${meta.source}` : ''}`
    : `${label || '—'} · GRCh38`

  return (
    <div style={{ display: 'flex', width: '100%', minHeight: '100vh', background: colors.pageBg, color: colors.textPrimary }}>
      <StructVarSidebar product={product} onSwitchProduct={onSwitchProduct} screen={screen} onNavigate={setScreen} samples={samples} sample={sample} onSelectSample={setSample} reference={meta?.reference} />
      <div style={{ flex: 1, display: 'flex', flexDirection: 'column', minWidth: 0 }}>
        <StructVarTopBar title={TITLES[screen]} subtitle={subtitle} />
        <div style={{ flex: 1, overflowY: 'auto' }}>
          {error && <div style={{ padding: 28, color: colors.statusPathogenic }}>Failed to load: {error}</div>}
          {loaded && !sample && <NoSamples />}
          {sample && status && status.state !== 'ready' && <Building status={status} onRetry={() => setRetry((n) => n + 1)} />}
          {sample && ready && screen === 'overview' && <OverviewScreen sample={sample} onOpenSvs={openSvs} />}
          {sample && ready && screen === 'svs' && <SvsScreen key={`${sample}-${preset?.nonce ?? 0}`} sample={sample} initial={preset?.filters ?? DEFAULT_FILTERS} />}
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
          <div style={{ fontSize: 15, fontWeight: 600, color: colors.statusPathogenic }}>The SV index could not be built</div>
          <div style={{ fontFamily: mono, fontSize: 12.5, marginTop: 8, wordBreak: 'break-word', lineHeight: 1.6 }}>{status.error}</div>
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
        <div style={{ fontSize: 15, fontWeight: 600 }}>Building the SV index</div>
        <div style={{ color: colors.textSecondary, fontSize: 13, marginTop: 6, lineHeight: 1.6 }}>
          The VCF is read once and each call is annotated with gene impact, population matches and genome tracks. It takes a few seconds; the index is kept
          and reused until the VCF or an annotation resource changes.
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
      <div style={{ fontSize: 16, fontWeight: 600 }}>No structural-variant VCFs found</div>
      <div style={{ color: colors.textSecondary, marginTop: 8, fontSize: 13.5 }}>
        Put an SV VCF (sawfish, pbsv or sniffles; single sample or a joint call across a family) at <code>&lt;root&gt;/&lt;callset&gt;/&lt;name&gt;.vcf.gz</code> and point{' '}
        <code>STRUCTVAR_RESULTS_ROOT</code> at <code>&lt;root&gt;</code> (default <code>~/sv_data/samples</code>). An optional <code>sample.json</code> beside it can give a
        label, the <code>caller</code>, and the <code>proband</code> and <code>parents</code> (sample names) when the family structure should not be inferred.
      </div>
    </div>
  )
}
