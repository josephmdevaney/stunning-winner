# Paraphase suite

The Paraphase suite shows what PacBio's [Paraphase](https://github.com/PacificBiosciences/paraphase) found in a sample, in three parts:

| Part | Source | Speed |
|---|---|---|
| **Regions / Region Calls** | The sample's Paraphase JSON, read directly | Instant, no network |
| **Annotation / Annotation Detail** | PacBio's `interpret_paraphase.py`, run against Ensembl VEP and ClinVar | Minutes to over an hour per sample |
| **Tool Versions** | Version headers in the loaded data, against each PacBio tool's newest GitHub release | Instant (releases cached) |

Backend: `backend/paraphase_annotate/` (routes under `/api/paraphase`) and `backend/pacbio_tools/` (routes under `/api/tools`). Frontend: `web/src/paraphase/`. Tests: `backend/tests/test_paraphase_calls.py`, `test_paraphase_regions.py`, `test_paraphase_runner.py`, `test_pacbio_tools.py`.

## Regions and Region Calls (the JSON view)

`paraphase_annotate/regions.py` turns a Paraphase JSON into one row per region: genes, total and gene copy number, haplotype count (and how many Paraphase inferred at two copies), median depth and its ratio to genome depth, a one-line headline for regions that have one (SMN1/SMN2 copy numbers, HBA genotype, RCCX phasing), and flags.

Everything is a field Paraphase wrote. Nothing is interpreted: no pathogenicity and no carrier call. A flag only restates the JSON:

- **Fusion**: `fusions_called` lists a fusion haplotype.
- **Structural variant**: `sv_called` is non-empty (HBA, F8).
- **Phasing failed**: `phasing_success` is false.
- **Copy number not reported**: the SMN region has no `smn1_cn`.
- **No data**: no depth and no phased haplotypes.

Region Calls adds the haplotype table (what each copy is assigned to, boundary, truncation, haplogroup, unique reads, variants), fusion calls, alleles and the region's scalar fields. Variants are capped at 25 per haplotype with a count; Paraphase's per-region VCF is the source for accurate small-variant calls.

"Medical" means the 14 regions `interpret_paraphase.py` annotates by default (`runner.MEDICAL_REGIONS`).

## Two JSON layouts

Paraphase v4.0.0 (the Rust rewrite) moved every region-specific field (`smn1_cn`, `fusions_called`, `gene_cn`, `annotated_alleles`, ...) under a `region_specific_info` sub-dict. Earlier output has them flat. `paraphase_annotate/calls.py` reads both, and `to_v4()` converts the flat layout to the nested one. The rule is the one in Paraphase's `docs/json.md`: the documented general fields stay where they are and everything else is region-specific, so a field a later Paraphase adds is handled correctly without a code change.

## Annotation

The annotation run is PacBio's own script, vendored (never committed: Paraphase's license grants binary redistribution only) by:

```
python backend/scripts/setup_paraphase_annotate.py            # the pinned tag (v4.1.0)
python backend/scripts/setup_paraphase_annotate.py --latest   # newest GitHub release
python backend/scripts/setup_paraphase_annotate.py --ref v3.5.0 --force
```

`paraphase_annotate/vendor/VERSION.json` records the tag, fetch time and checksums. A different tag is re-fetched, not skipped.

**v4.1.0 on any Paraphase output.** The v4.x script reads the nested layout, so each run writes a normalized copy of the JSON to `<results>/<sample>/input/`, beside a symlink to the real `<sample>_paraphase_vcfs/` directory (the script finds VCFs that way). Nothing under the sample's own directory is written.

**Python.** The v4.1.0 script needs Python 3.12 or later (nested double-quoted f-strings), with `requests` and `pyyaml`. The runner finds one (`$PARAPHASE_ANNOTATE_PYTHON`, then python3.14/3.13/3.12 on PATH and the usual Homebrew, conda and framework locations); the API server itself can stay on an older Python. If none is found the run is "blocked" with the reason.

**Provenance.** Each run writes `run_meta.json`: the script release, the Python used, whether the JSON had to be adapted, which Paraphase produced the data (read from its VCF header; the JSON records no version), and which default regions the JSON doesn't contain (the script skips those with a warning). The Annotation screen shows this, and offers a re-run when a newer script than the one that made the results is vendored.

**Re-runs keep the old results.** A run can take over an hour of rate-limited lookups, so starting a new one moves the previous results to `<results>/<sample>/archive/<UTC timestamp>/` instead of deleting them.

**Hung requests.** PacBio's script looks each variant up on Ensembl and ClinVar through `requests` with no timeout, so a connection that goes quiet half way is waited on forever (a real run sat on one for 17 hours, using two seconds of CPU). The runner therefore starts the script through `paraphase_annotate/launch.py`, which gives every request that has no timeout of its own a 15-second connect and a 90-second read timeout (`PARAPHASE_REQUEST_CONNECT_TIMEOUT` and `PARAPHASE_REQUEST_READ_TIMEOUT` change them). The script's own retry logic then sees a timeout like any other failed request, so a bad connection is retried and, if it stays bad, the run fails and says so in its log. The values used are recorded in `run_meta.json`. `socket.setdefaulttimeout` does not do this: `requests` passes urllib3 an explicit `timeout=None`, which overrides it.

## Tool Versions

`pacbio_tools/` compares, for each PacBio tool, the version found in the data against the newest GitHub release.

- **Observed versions** come from the files, never assumed: `##trgtVersion` in TRGT VCFs, `##source="sawfish x.y.z"` in sawfish VCFs, `##paraphase_version` in Paraphase VCFs, the recorded caller string for MitoRSaw, and the vendored script's tag. A tool whose version no loaded file records (pbmm2: it is in the BAM header, which this suite doesn't read) shows "not recorded". PTCP's version is documented for the Carrier cohort rather than read from a file, and is labelled that way.
- **Newest releases** come from three places, freshest first: a live GitHub query, an on-disk cache of the last one (6 hours), and `pacbio_tools/releases_snapshot.json` (committed, so an offline machine still has something to compare against; the page says when it is using it). "Check GitHub now" forces a live query. Set `GITHUB_TOKEN` to raise GitHub's unauthenticated rate limit. One tool failing to fetch keeps its previous release and is reported; it doesn't blank the page.
- **Gap** is major / minor / patch behind. A tool is "Behind" when any sample's data came from an older release.

The registry (`pacbio_tools/registry.py`) lists the tools in this suite's data flow (Paraphase, TRGT, sawfish, MitoRSaw, PTCP, pbmm2) and, separately, other PacBio tools (HiPhase, pb-CpG-tools, pbsv, HiFiCNV, pb-StarPhase, pbaa, pbfusion, pbtk) so the whole toolchain's newest releases are in one place. Adding a tool is one entry there, plus a reader in `observed.py` and a line in `router.collect()` if its version can be read from a file.

Refresh the committed snapshot with `python backend/scripts/refresh_pacbio_releases.py`.

## Validation, and what wasn't validated

- **Offline equivalence (done).** The v4.1.0 script, run on a v3.5.0 sample's JSON after `to_v4()`, with the lookups answered locally so nothing left the machine, completed all 14 default regions and reproduced everything derived from the JSON in the existing v3.5.0 results: all 80 summary cells, all 56 haplotype rows (gene, transcript, boundaries, boundary check) and the known-nonfunctional call.
- **Not validated.** The live Ensembl/ClinVar annotation with the v4.1.0 script has not been run end to end, so the pathogenic-variant columns it produces are untested with this version. No real Paraphase v4.x JSON was on hand either (the data on hand was made by Paraphase 3.5.0); the v4 layout was implemented from Paraphase's documentation and the v4.1.0 script's source, and exercised with converted v3 output and synthetic v4 JSON.
- **Unit tests** cover the schema adapter (including a field the code has never seen), the runner (stand-in script: adapted input, symlink, provenance, archive-on-rerun, blocked states, v3 script compatibility), the region table (both layouts give the same table) and the versions package (parsing, readers, status rules, release fallbacks).

## Limits

- Paraphase v4.1.0's `--file-prefix` option can rename the output files. The sample registry (`samples.py`) points at the JSON path directly, but the annotation script and `vcfs_dir()` expect the default `<sample>.paraphase.json` and `<sample>_paraphase_vcfs/` names.
- New v4.1.0 targets (CATSPER2, CLCNKB, TPSAB1) appear in the Regions table when a JSON contains them, but they are not in the default annotation set, which stays at the script's own 14 regions.
- The Carrier cohort's per-region Paraphase VCFs aren't loaded here, so Paraphase's version for that cohort shows as not recorded.
- Version numbers are compared by their first `x.y.z`; a build suffix (`0.2.7-8ee9b5b`) is ignored, so two builds of one release look identical.
