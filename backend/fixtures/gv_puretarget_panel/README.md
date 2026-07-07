# GV PureTarget panel fixture

`GV-01.trgt.vcf` is real `trgt` v3.0.0 output from PacBio's commercial **PureTarget
Repeat Expansion Panel 2.0** -- a targeted 38-locus clinical repeat-expansion panel,
run with `--preset targeted` against `PureTarget_repeat_expansion_panel_2.0.repeat_definition.GRCh38.bed`
(per the VCF's own `##trgtCommand` header line). Provided by the user as one of 48
samples (GV-01 through GV-48) from the same panel run.

Since PureTarget's panel covers 38 of STRchive's 80 loci, the other 42 correctly show
up as "No call" in RepeatScope rather than fabricated values -- this is real, partial
panel coverage, not a bug.
