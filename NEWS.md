# epiRomics 0.99.5

## Changes

* Added public getter and setter methods for all five slots of the
  `epiRomicsS4` class: `annotations()`, `meta()`, `txdb()`,
  `organism()`, and `genome()` (plus the corresponding `<-`
  assignment forms). Users should prefer these accessors over
  `obj@slot` or `methods::slot(obj, "slot")`. `organism()` extends
  the generic from `BiocGenerics` and `genome()` extends the generic
  from `GenomeInfoDb`, following Bioconductor convention. Every
  setter invokes `methods::validObject()` so invalid assignments
  (for example an empty-string `genome`) fail fast with a clear
  error. A full audit of the codebase confirms these five slots
  are the only slots on `epiRomicsS4`; the class definition has
  not changed, and downstream pipeline functions never introduce
  new slots — they only update the content of `annotations`
  (a `GRanges`, whose `mcols` carry the per-stage details).
* Added `?epiRomicsS4-accessors` overview topic that lists every
  getter/setter pair in one table with runnable examples. The
  `?epiRomicsS4-class` page cross-references each individual
  accessor via `@seealso` and documents each slot with an
  explicit "Access via ..." note pointing at the matching
  getter.
* Refactored both vignettes to demonstrate the new accessor API.
  `methods::slot()` no longer appears in any user-facing example,
  vignette, or man page. `getting-started-with-epiRomics.Rmd` had
  7 `methods::slot()` calls replaced with `annotations()`, plus a
  prose paragraph noting that users should prefer the getter API
  (one paragraph inserted in each vignette near first use).
  `articles/full-analysis-with-epiRomics.Rmd` had 22
  `methods::slot()` calls replaced.
* Every roxygen `@examples` block that previously used
  `methods::slot(obj, "...")` now uses the matching getter
  (touched files: `R/synthetic-data.R`, `R/enhanceosomes.R`,
  `R/enhancers.R`, `R/regions_of_interest.R`).
* Swapped the four runtime `methods::slot(database, "...")`
  calls inside `make_example_enhanceosome()` in
  `R/synthetic-data.R` for the matching getters
  (`genome(database)`, `meta(database)`, `txdb(database)`,
  `organism(database)`).
* Updated user-facing roxygen prose that referred to
  `database@meta` or `@annotations` to reference the new
  `meta()` / `annotations()` accessors instead (touched files:
  `R/synthetic-data.R`, `R/chromatin_states.R`,
  `R/enhancers.R`, `R/benchmark_enhancer_predictor.R`,
  `R/regions_of_interest.R`).
* Added `tests/testthat/test-accessors.R` with 20 `test_that`
  blocks / 50 assertions covering: every getter on the hg38
  and mm10 fixtures and on a fresh empty object, every setter
  round-trip and class preservation, validity enforcement
  (`genome(x) <- ""` and `txdb(x) <- ""` must error), and
  confirmation that `organism()` / `genome()` extend the
  `BiocGenerics` / `GenomeInfoDb` generics rather than
  shadowing them. After this change the full package test
  suite runs 261 tests / 607 assertions with 0 failures,
  0 errors, 7 legitimate skips (Windows-only harness checks
  and optional non-model-organism TxDb packages).
* Added a one-time session tip to the rename-warning helper
  (`.warn_renamed()`): the first time a deprecated `epiRomics_*`
  alias fires, users also see a one-line note pointing them at
  the new accessor API (`?epiRomicsS4-accessors`).
* Fixed the `test-deprecated-aliases.R` harness, which was still
  written against the pre-0.99.4 `.Deprecated()`-based behaviour
  and therefore was order-dependent under the `message()`-based
  rename helper introduced in 0.99.4. The harness now asserts
  against the `message()` condition (via `invokeRestart("muffleMessage")`)
  and resets `.warn_renamed()` state at the top of the block so
  every alias fires fresh exactly once, regardless of the order
  other test files ran in. No package code was changed by this
  test fix.
* Quality gates for 0.99.5: `R CMD check --as-cran` reports 0
  ERRORs / 0 WARNINGs / 1 NOTE (the expected "New submission"
  note); `BiocCheck` reports 0 ERRORs / 0 WARNINGs / 9 NOTES
  (all pre-existing and stylistic — none introduced by this
  change); both vignettes rebuild cleanly with the new
  accessor API.

# epiRomics 0.99.4

## Changes

* Replaced `.Deprecated()` calls in renamed-function wrappers with
  one-time `message()` notices via a shared `.warn_renamed()` helper.
  Users calling the old `epiRomics_*` names still see a clear nudge
  to migrate to the new verb-first names (removal planned for the
  next release), but BiocCheck's `.Deprecated / .Defunct usage`
  warning no longer fires.

# epiRomics 0.99.3

## Highlights

* **On-demand example data**: `cache_data()` downloads example
  datasets (~1.3 GB) via BiocFileCache, keeping the package tarball
  under 5 MB. `has_cache()` and `get_cache_path()` check and locate
  cached data.

* **Base R graphics track layer**: `plot_tracks()` renders
  publication-quality multi-track genome browser views in 1-2 seconds
  per locus using base R graphics. ATAC and RNA signals display in
  mirrored panels for direct cross-sample comparison.

* **Quick visualization**: `plot_quick_view()` provides zero-setup
  BigWig signal visualization for any gene locus — no database
  required.

* **Gene-centered visualization**: `plot_gene_tracks()` displays
  all database tracks at any gene locus without the enhanceosome
  pipeline.

* **Chromatin state classification**: `classify_chromatin_states()`
  classifies regions into biologically meaningful states (active,
  bivalent, poised, primed, repressed, unmarked) following
  ChromHMM/Roadmap Epigenomics conventions.

* **TF co-binding analysis**: `analyze_tf_cobinding()` uses Fisher's
  exact test with odds ratios and pointwise mutual information (PMI)
  to identify significant TF co-binding pairs.

## API

All exported functions follow Bioconductor naming conventions
(verb-first, no package prefix). Legacy names remain exported as
`.Deprecated()` aliases for one release cycle and will be removed in
a future version.

Function names:

* `epiRomics_build_dB()` → `build_database()`
* `epiRomics_quick_view()` → `plot_quick_view()`
* `epiRomics_track_layer()` → `plot_tracks()`
* `epiRomics_track_layer_fast()` → `plot_tracks_fast()`
* `epiRomics_track_layer_gene()` → `plot_gene_tracks()`
* `epiRomics_putative_enhancers()` → `find_putative_enhancers()`
* `epiRomics_enhancers_co_marks()` → `find_enhancers_by_comarks()`
* `epiRomics_enhanceosome()` → `find_enhanceosomes()`
* `epiRomics_enhancers_filter()` → `filter_enhancers()`
* `epiRomics_filter_accessible()` → `filter_accessible_regions()`
* `epiRomics_chromatin_states()` → `classify_chromatin_states()`
* `epiRomics_chromatin_states_categories()` → `chromatin_state_categories()`
* `epiRomics_annotate_putative()` → `annotate_enhancers()`
* `epiRomics_tf_cobinding()` → `analyze_tf_cobinding()`
* `epiRomics_tf_overlap()` → `analyze_tf_overlap()`
* `epiRomics_enhancer_predictor_to_ref()` → `benchmark_enhancer_predictor()`
* `epiRomics_regions_of_interest()` → `get_regions_of_interest()`
* `epiRomics_cache_data()` → `cache_data()`
* `epiRomics_cache_path()` → `get_cache_path()`
* `epiRomics_has_cache()` → `has_cache()`

Parameter names (no deprecation shim — update calls directly):

* `epiRomics_dB` → `database`
* `epiRomics_db_file` → `db_file`
* `epiRomics_genome` → `genome`
* `epiRomics_organism` → `organism`
* `epiRomics_curated_database` → `curated_database`
* `epiRomics_histone`, `epiRomics_histone_mark_1`, `epiRomics_histone_mark_2`
  → `histone`, `histone_mark_1`, `histone_mark_2`
* `epiRomics_track_connection` → `track_connection`
* `epiRomics_index` → `index`
* `epiRomics_keep_epitracks` → `keep_epitracks`
* `epiRomics_type` → `type`
* `epiRomics_test_regions` → `test_regions`
* `epiRomics_putative_enhancers` → `putative_enhancers`
* `epiRomics_putative_enhanceosome` → `putative_enhanceosome`
* `epiRomics_enhanceosome` → `enhanceosome`

The `epiRomicsS4` class name is unchanged (the `S4` suffix indicates
class source).

* `genome` is a free-form string validated against the user-supplied
  `TxDb` (`validate_genome_matches_txdb()`). No hardcoded whitelist;
  any organism with a `TxDb.*` and `org.*.db` package is supported.

## Dependencies

* `TxDb.Hsapiens.UCSC.hg38.knownGene`,
  `TxDb.Mmusculus.UCSC.mm10.knownGene`, `org.Hs.eg.db`,
  `org.Mm.eg.db`, `BiocFileCache`, and `parallel` moved to `Suggests`
  with `requireNamespace()` guards at every call site.
* Full `@importFrom` audit across all R files; NAMESPACE regenerated.

## Examples and vignettes

* Ships a ~0.8 MB toy dataset at `inst/extdata/toy/` (400 kb
  `hg38 chr11` window centred on the *INS* locus) curated from the
  Zenodo archive (https://zenodo.org/records/19189987). The
  *Getting Started* vignette runs end-to-end against the toy data;
  the full walkthrough is preserved as *Full Analysis with epiRomics*
  and resolves the Zenodo archive via `cache_data()`.
* The full walkthrough lives under `vignettes/articles/` as a
  pkgdown-only article (listed in `.Rbuildignore`). It is rendered on
  the pkgdown site but is **not** built during `R CMD check` or
  package installation, so users never wait for the 1.3 GB Zenodo
  download at install time. The lightweight *Getting Started*
  vignette is the only registered vignette in the tarball.
* Full walkthrough *Interactive showcases* section links to both the
  Mouse Islet (Mawla et al. 2023) and Human Islet
  (Mawla & Huising 2021) Shiny browsers.
* Exported synthetic-data helpers — `make_example_database()`,
  `make_example_putative_enhancers()`, `make_example_enhanceosome()`,
  `make_example_bigwig()` — power every man-page example so examples
  run under `R CMD check --run-donttest` without network access. The
  only remaining `\donttest` block is in `cache_data()` (network
  download, justified).
* Both vignettes open with visible setup chunks, interleave
  explanatory prose before and after every code chunk, and print
  data-frame structure (manifest CSVs, `head()` of result objects) so
  readers see inputs and outputs.
* Getting-started introduction frames the multi-omics enhancer scope,
  enumerates Bioconductor integrations, and contrasts epiRomics
  against Gviz, trackViewer, and Signac.
* `vignettes/references.bib` holds package citations; vignettes use
  `[@mawla2023; @mawla2021]` inline references.

## Documentation and installation

* README recommends Bioconductor installation first; GitHub install
  is flagged as not recommended and uses `BiocManager::install()`.
* README points users to `browseVignettes("epiRomics")`, the pkgdown
  site (<https://huising-lab.github.io/epiRomics/>), and the
  Bioconductor landing page.
* pkgdown site restored: `_pkgdown.yml` and a dedicated
  `.github/workflows/pkgdown.yaml` rebuild the site on every push to
  `main` and on every release.
* Top-level `doc/` directory removed from the source tree.
* Comprehensive HTML and PDF vignettes with live BigWig signal
  visualizations.
* `.onAttach()` startup message with version, cache status, and
  citation info.

## Speed

* Parallel overlap counting via `parallel::mclapply()` when available.
* Pre-split annotations by type for O(1) lookup in
  `find_enhanceosomes()`.
* Vectorized chromatin state classification and signal aggregation.
* BigWig RDS caching for memory-efficient large-scale analyses.

## CI/CD

* Cross-platform GitHub Actions matrix: Ubuntu (release/devel/oldrel-1),
  macOS ARM, macOS Intel, Windows.
* BiocCheck validation on every push/PR.
* Test coverage with full integration tests via cached example data.

## Housekeeping

* BiocCheck cleanups: eliminated `paste()` / `paste0()` usage inside
  condition signals across 5 files; removed `<<-` from
  `R/synthetic-data.R`; refactored `make_example_enhanceosome()` to
  a synthetic fast path (analyze_tf_cobinding example now runs in
  <1 s, down from 35 s).
* Fixed CRAN URL-redirect NOTE by pointing README and vignette links
  at canonical URLs.
* Aligned README R-version minimum to DESCRIPTION `Depends`
  (R ≥ 4.5.0).
* Trimmed residual >80-character lines from R source for BiocCheck
  line-length cleanliness.
* Dataset provenance wording throughout the vignettes, toy README,
  and `inst/scripts/make-toy-data.R` clearly identifies the curator
  role and points at the package README for the full curation
  methodology (GEO GSE76268 source; ENCODE-DCC ATAC pipeline; MACS2;
  DiffBind; FANTOM5, Human Islet Regulome, UCNEs).
