# Carrier panel fixture

Two kinds of samples, clearly separated by naming:

## `m84036_250714_174336_s1.bc2113` -- real data

Real PacBio HiFi carrier-panel output for one sample (`bc2113`, from the
`Carrier_Coriell_16plex_RevioSRPQ` Coriell 16-plex Revio run, PureTarget Carrier
Panel pipeline `ptcp` v3.1.0), provided by the user. Not simulated, not hand-built.

- **`bc2113.trgt.vcf`** -- real `trgt` v5.1.0 output, genotyped with proper male
  hemizygous chrX calling (`--sex M`, real single-allele `GT` values like `1`, not
  a diploid-padded `1/1`). Includes a genuine full-mutation finding: `FXS_FMR1` at
  731 CGG repeats, hemizygous, well past the 201-repeat pathogenic threshold --
  exactly the kind of real, unfabricated finding this project is built to surface
  honestly. Parsing this file's hemizygous `GT` values surfaced a real bug (a
  single-allele `GT` was silently treated as an ungenotyped diploid call); see the
  `hemizygous` handling in `repeatscope/trgt.py` and its regression tests.
- **`paraphase/m84036_250714_174336_s1.bc2113.paraphase.json`** -- a trimmed export
  of PacBio Paraphase's real per-gene JSON output for this sample, covering the 6
  non-repeat carrier targets Paraphase itself calls (`gba`, `cyp21`, `hba`, `hbb`,
  `rpgr`, `smn1`; TNXB and F8's inversion are out of scope for Paraphase's own copy
  number caller -- see below). Only the summary fields RepeatScope actually
  consumes are kept (`total_cn`, `gene_cn`, `genotype`, `fusions_called`,
  `sv_called`, `smn1_cn`/`smn2_cn`, ...); the full real output also includes
  per-read alignment provenance (thousands of CCS read IDs and per-site genotype
  strings per haplotype) that isn't needed here and was dropped to keep the fixture
  reviewable.
- **`paraphase/m84036_250714_174336_s1.bc2113.f8inversion.vcf`** -- real output
  from PacBio's dedicated F8 intron-1/22 inversion caller (a separate, purpose-built
  tool from Paraphase's general copy-number caller). Reports `GT=0`/`SUPP=0` at
  both breakpoints for this sample: no inversion detected.

See `repeatscope/paraphase.py` for how these are turned into Overview findings --
every finding is read directly off a field one of these two real tools already
asserts (never a pathogenic-variant call re-derived from raw sequence here).

`sample_sex.json`'s "M" for this sample is real, from the run's own
`sample_sheet.csv` (`bam_id,sex` columns) -- and matches this VCF's own hemizygous
chrX calling, unlike HG002's fixture (see that fixture's README).

## `m84036_250809_035451_s4.bc2192` -- real data

A second real PacBio HiFi carrier-panel sample (`bc2192`, from the same Coriell
16-plex Revio run as `bc2113`), provided by the user with substantially richer
per-gene data than `bc2113` had. Also not simulated, not hand-built.

- **`bc2192.trgt.vcf`** -- real `trgt` v5.1.0 output, genotyped `--sex F`
  (confirmed against the F8 inversion caller's own command line). All 5 repeat
  loci genotype Normal.
- **`paraphase/m84036_250809_035451_s4.bc2192.paraphase.json`** -- trimmed the
  same way as bc2113's, above.
- **`paraphase/m84036_250809_035451_s4.bc2192.f8inversion.vcf`** -- real F8
  inversion caller output: `GT=0/0`/`SUPP=0` at both breakpoints (no inversion).
  Because this sample is female, chrX is reported as a diploid `0/0`, not the
  hemizygous `0` bc2113's male sample has -- parsing that surfaced a real bug
  (`_classify_f8` only recognized hemizygous reference GTs, so a diploid `0/0`
  was misread as an inversion); see `_gt_has_alt_allele` in `repeatscope/paraphase.py`
  and its regression test.
- **`paraphase/m84036_250809_035451_s4.bc2192.havanno.json`** -- real output from
  PacBio's "havanno" per-haplotype variant annotation (references PacBio's
  Gauchian GBA1 caller and a curated CYP21A2 variant list), kept in full since it's
  already free of per-read bloat. For each Paraphase gene, lists the named
  pathogenic (and pseudogene-marker) variants found on each haplotype-phase call.
  This sample has two real, clinically notable findings surfaced entirely from
  named variants the tool itself asserts:
  - **CYP21A2**: 3 of 8 haplotype-phase calls (`cyp21_hap2`/`hap5`/`hap7`) all
    carry the identical set of 9 named variants classically associated with
    congenital adrenal hyperplasia (21-hydroxylase deficiency), including
    `p.Gln319Ter` and the `Val238Cluster` variants.
  - **GBA**: 4 of 7 haplotype-phase calls carry overlapping sets of 13-14 named
    variants associated with Gaucher disease, including several `Gauchian_Exon9-11`
    markers characteristic of a GBA1/GBAP1 recombinant allele.

  `repeatscope/paraphase.py`'s `_classify_havanno_gene` reports these as
  "Full expansion" (pathogenic) findings, listing the deduplicated union of named
  variants and how many of the gene's haplotype-phase calls carry them --
  collapsing redundant phasing views of what's likely 1-2 physical haplotypes,
  without guessing which. Genes with no named pathogenic variants (F8, HBA, HBB,
  RPGR, SMN1) fall back to the same `total_cn`/copy-number classifier bc2113 uses.
- **`paraphase/m84036_250809_035451_s4.bc2192.sv.vcf`** -- real `sawfish` v2.0.3
  structural-variant output, trimmed to its 9 chr16 records (this sample's sawfish
  output only covers the HBA1/HBA2 region; see `carrier_screen.non_repeat_targets.json`'s
  own HBA1/2 entry, which already documents "large deletion via sawfish" as that
  target's mechanism). Paraphase's own HBA copy-number caller couldn't resolve a
  genotype for this sample (`genotype: null`, `total_cn: 5` against an expected 4),
  but sawfish independently found a real heterozygous 321bp deletion at
  `chr16:163481-163802` overlapping the region. `_classify_hba` now falls back to
  sawfish when Paraphase's own genotype is unresolved, and reports this as
  "Uncertain significance" -- the deletion is real, but nothing here attributes it
  to HBA1 vs. HBA2 specifically, so it is not guessed as either benign or
  pathogenic.

`sample_sex.json`'s "F" for this sample is real, from the F8 inversion caller's
own `--sex F` command-line invocation recorded in its VCF header.

## `m84036_250809_035451_s4.bc2108` -- real data

A third real sample from the same Coriell Nanobind 96-plex Revio run
(`Carrier_Nanobind_96plex_RevioSRPQ`), pulled from the run's full output set (all
96 samples' TRGT/Paraphase/havanno/F8-inversion output, made available in bulk).
Selected specifically because its `havanno.json` shows two clean, textbook
heterozygous-carrier findings -- unlike bc2192's atypically elevated CYP21A2 copy
number, this sample's copy numbers are entirely normal, so the pathogenic finding
is unambiguous rather than confounded by a structural complication:

- **CYP21A2**: `total_cn: 4` (normal), 2 of the locus's 4 haplotype-phase calls
  (`cyp21_hap1`/`hap2`) carry the identical 9-variant CAH-associated set also seen
  in bc2192 -- `_havanno_zygosity` resolves this to "Heterozygous" (2 of 4 total
  copies, not all of them).
- **GBA**: `total_cn: 4` (normal), 2 of 4 haplotype-phase calls
  (`gba_hap1`/`hap2`) carry the identical 14-variant Gaucher-associated set --
  also resolves to "Heterozygous".
- **HBA1/2**: unlike bc2113 (`aa/aa`, normal) and bc2192 (genotype unresolved),
  Paraphase itself resolved a genotype here: `-a/aa` -- one copy carrying the
  classic single-alpha-globin-gene deletion (silent-carrier alpha thalassemia),
  the other normal. Its own `sv_called` field independently names the deletion:
  `"hba_3p7delhap1": "3p7del,chr16:172240-176044"` (the well-known "-alpha-3.7"
  deletion allele). `_classify_hba`'s existing non-`aa/aa`-genotype branch reports
  this as "Uncertain significance", surfacing the real genotype and deletion
  verbatim rather than asserting a specific clinical read on a single-gene
  deletion this project has no dedicated classifier for.

Only this sample's `trgt.vcf`, `f8inversion.vcf`, `havanno.json`, and a trimmed
`paraphase.json` were pulled in (no sawfish VCF -- not needed since Paraphase's
own HBA genotype/`sv_called` already resolved a finding here, unlike bc2192).

`sample_sex.json`'s "F" is real, from this sample's own F8 inversion caller
`--sex F` command-line invocation, same as bc2192.

## `CR-01` .. `CR-10` -- synthetic

Hand-built, spec-accurate synthetic `trgt` v5.1.0 VCFs (same approach as
`fixtures/synthetic_trgt.vcf`) covering just this panel's 5 repeat loci
(`FRDA_FXN`, `PRTS_ARX`, `EIEE1_ARX`, `FXS_FMR1`, `FRAXE_AFF2`), generated by
`backend/scripts/generate_carrier_panel_vcfs.py`. No real carrier-panel BAM was
available for a 16-sample panel beyond `bc2113`, so these exist to give the Carrier
mode UI a realistic spread of scenarios (clean negative screens, intermediate
carriers, full-pathogenic carriers, a no-call) that a single real sample can't
cover alone. They have no corresponding Paraphase data, so their non-repeat targets
correctly show "Not applicable" in the Overview -- the same honest fallback used
before any real Paraphase data existed.

`sample_sex.json` labels all ten "F": their 5 repeat loci (including the chrX ones)
were generated as diploid two-allele calls, which is only consistent with a female
(XX) sample -- see `generate_carrier_panel_vcfs.py`.
