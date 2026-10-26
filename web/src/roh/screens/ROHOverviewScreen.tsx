import { useState } from 'react'
import { colors, mono } from '../../theme'
import { StatCard, StatGrid } from '../../components/StatCard'
import type { RohOverviewResponse, RohRegion } from '../api'
import { ConsanguinityCard, IdeogramCard, SexInferenceBanner } from '../components'

export function ROHOverviewScreen({ overview, allRegions }: { overview: RohOverviewResponse; allRegions: RohRegion[] }) {
  const [ideogramFilter, setIdeogramFilter] = useState('all')
  const s = overview.stats

  return (
    <div style={{ padding: 28, display: 'flex', flexDirection: 'column', gap: 20 }}>
      <SexInferenceBanner sexInference={overview.sex_inference} />

      <StatGrid>
        <StatCard
          label="Autosomal F (inbreeding coefficient)"
          value={`${s.f_value_pct}%`}
          sub={`ROH ≥${s.f_value_min_length_mb} Mb, ${overview.genome_build} autosomes`}
        />
        <StatCard label="Autosomal ROH total" value={`${(s.total_autosomal_roh_bp / 1_000_000).toFixed(1)} Mb`} sub={`${s.n_regions_autosomal.toLocaleString()} regions`} />
        <StatCard
          label="UPD screening flags"
          value={String(s.n_upd_flags)}
          sub="Heuristic, not diagnostic"
          valueStyle={s.n_upd_flags > 0 ? { color: colors.statusIntermediate } : undefined}
          dotColor={s.n_upd_flags > 0 ? colors.statusIntermediate : undefined}
        />
        <StatCard label="AR candidate genes" value={String(s.n_ar_candidate_genes)} sub={`Within ROH ≥${s.ar_candidate_min_length_mb} Mb`} />
      </StatGrid>

      {overview.upd_flags.length > 0 && (
        <div style={{ background: colors.interpBg.intermediate, border: `1px solid ${colors.statusIntermediate}`, borderRadius: 12, padding: '16px 20px' }}>
          <div style={{ display: 'flex', alignItems: 'center', gap: 8, fontSize: 14, fontWeight: 600, color: colors.interpColor.intermediate }}>
            <div style={{ width: 8, height: 8, borderRadius: '50%', background: colors.statusIntermediate }} />
            Possible UPD — {overview.upd_flags.length} region{overview.upd_flags.length > 1 ? 's' : ''} flagged
          </div>
          {overview.upd_flags.map((f) => (
            <div key={`${f.chrom}:${f.start}`} style={{ fontSize: 12.5, color: colors.textBody, marginTop: 6, fontFamily: mono }}>
              {f.chrom}:{f.start.toLocaleString()}–{f.end.toLocaleString()} ({f.length_mb} Mb, {f.position}) — {f.upd_flag}
            </div>
          ))}
          <div style={{ fontSize: 11.5, color: colors.textSecondary, marginTop: 8, lineHeight: 1.5 }}>{overview.upd_methodology_note}</div>
        </div>
      )}

      <div style={{ display: 'grid', gridTemplateColumns: '1.55fr 1fr', gap: 20, alignItems: 'start' }}>
        <IdeogramCard overview={overview} allRegions={allRegions} filter={ideogramFilter} onFilterChange={setIdeogramFilter} />

        <div style={{ display: 'flex', flexDirection: 'column', gap: 16 }}>
          <ConsanguinityCard overview={overview} allRegions={allRegions} />

          <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, padding: '18px 20px' }}>
            <div style={{ fontSize: 13, fontWeight: 600, marginBottom: 4 }}>AR candidate genes</div>
            <div style={{ fontSize: 11.5, color: colors.textSecondary, marginBottom: 12 }}>
              Autosomal-recessive OMIM genes overlapping an ROH ≥{s.ar_candidate_min_length_mb} Mb (autozygosity-mapping candidates).
            </div>
            {overview.ar_candidate_genes.length === 0 ? (
              <div style={{ fontSize: 13, color: colors.textSecondary, padding: '12px 0' }}>None found for this sample.</div>
            ) : (
              <div style={{ display: 'flex', flexDirection: 'column', gap: 8 }}>
                {overview.ar_candidate_genes.map((g) => (
                  <div key={g.symbol} style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', fontSize: 12.5, borderBottom: `1px solid ${colors.borderSoftest}`, paddingBottom: 8 }}>
                    <div>
                      <span style={{ fontWeight: 700, fontFamily: mono }}>{g.symbol}</span>
                      <span style={{ color: colors.textSecondary, marginLeft: 8 }}>{g.disease}</span>
                    </div>
                    <div style={{ color: colors.textTertiary, fontFamily: mono, whiteSpace: 'nowrap' }}>
                      {g.chrom} · {g.roh_length_mb} Mb
                    </div>
                  </div>
                ))}
              </div>
            )}
          </div>
        </div>
      </div>

      <div style={{ background: 'white', border: `1px solid ${colors.border}`, borderRadius: 12, padding: '18px 20px' }}>
        <div style={{ fontSize: 11, color: colors.textSecondary, lineHeight: 1.5 }}>{overview.non_autosomal_note}</div>
      </div>
    </div>
  )
}
