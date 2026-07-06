# Genome-wide QC tab

_Started 2026-07-06 12:43 UTC_

---

## User

Can we make another tab that is a general QC of the genome run include information in the sample QC but make it genome wide

## Assistant

_[tool: read_file]_

_[tool: read_file]_

_[tool: dc_write]_

Now wiring it into the nav and screen switch.

_[tool: str_replace_edit]_

_[tool: ready_for_verification]_

Added a "Run QC" tab: genome-wide coverage distribution, GC-bias scatter, per-chromosome coverage across all 24 chromosomes, alignment/read metrics, and a per-sample summary table for the run.

## User

Save this design as a PDF: Repeat Expansion Analyzer.dc.html

## Assistant

_[tool: read_file]_

_[tool: write_file]_

_[tool: show_html]_

_[tool: get_webview_logs]_

_[tool: open_for_print]_

Print tab opened — it will auto-trigger the print dialog once fonts/images finish loading. Note: only the currently active screen (whichever tab was last selected) prints, since the others aren't rendered in the DOM.

## User

Can we make the run QC specific for one sample?

## User

make the mean coverage have a Y axis of depth

## User

make the mean coverage have a Y axis of depth

## Assistant

_[tool: dc_js_str_replace]_

_[tool: ready_for_verification]_

Added a depth-labeled Y axis (0–40x gridlines) to the per-chromosome coverage chart.

## User

Beautiful Cluade!

## Assistant

Glad it landed! Let me know if you want more tweaks.

## User

can you make a PDF of just the Run QC

## Assistant

_[tool: str_replace_edit]_

_[tool: open_for_print]_

Print tab opened, now defaulted to the Run QC screen — the print dialog should trigger automatically.

## User

move alligment & read metrics above genome-wide distribution

## Assistant

_[tool: read_file]_

_[tool: dc_html_str_replace]_

_[tool: ready_for_verification]_

Done — Alignment & read metrics now sits above the coverage distribution/GC bias row.

## User

can a PDF be made of run qc tab

## Assistant

_[tool: read_file]_

_[tool: write_file]_

_[tool: open_for_print]_

Yes — opened the print tab defaulted to the Run QC screen; use the print dialog's "Save as PDF" destination.

## User

can we use TRGT v5.1.0 https://github.com/PacificBiosciences/trgt/releases/tag/v5.1.0

## Assistant

_[tool: str_replace_edit]_

_[tool: ready_for_verification]_

Updated the sidebar tag to TRGT v5.1.0.

## User

make the caller for all TRGT v5.1.0

## Assistant

_[tool: grep]_

_[tool: dc_html_str_replace]_

_[tool: str_replace_edit]_

_[tool: ready_for_verification]_

All TRGT caller references now read v5.1.0.

## User

can we add another tab that is for the rest of the repeats from this bed file

## Assistant

_[tool: read_file]_

_[tool: read_file]_

_[tool: dc_write]_

Now wiring the new tab into the main nav/screen switch.

_[tool: str_replace_edit]_

_[tool: ready_for_verification]_

Added a "Repeat Catalog" tab covering the remaining ~171K genome-wide repeat regions from the bed file (beyond the 80 STRchive disease loci) — with chromosome/motif/length filters and pagination since the full set is too large to render at once.

## User

can we add calls for the sample to the repeat catalog

## Assistant

_[tool: dc_html_str_replace]_

_[tool: dc_js_str_replace]_

_[tool: dc_js_str_replace]_

Also updating the sub-copy on the stat card and note to make clear these are HG-2214-R3's genotype calls.

_[tool: dc_js_str_replace]_

_[tool: ready_for_verification]_

Added per-locus genotype calls (A1, A2, depth) for HG-2214-R3 to the Repeat Catalog table, seeded deterministically per region so they stay stable across pagination/filtering.

