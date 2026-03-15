# epiRomics 0.99.0

First major release.

## Highlights

* **On-demand example data**: New `epiRomics_cache_data()` downloads example
  datasets (~1.3 GB) via BiocFileCache, keeping the package tarball under 5 MB.
  `epiRomics_has_cache()` and `epiRomics_cache_path()` check and locate cached
  data.

* **Base R graphics track layer**: `epiRomics_track_layer()` renders
  publication-quality multi-track genome browser views in 1-2 seconds per locus
  using base R graphics. ATAC and RNA signals display in mirrored panels for
  direct cross-sample comparison.

* **Quick visualization**: New `epiRomics_quick_view()` provides zero-setup
  BigWig signal visualization for any gene locus — no database required.

* **Gene-centered visualization**: New `epiRomics_track_layer_gene()` displays
  all database tracks at any gene locus without the enhanceosome pipeline.

* **Chromatin state classification**: `epiRomics_chromatin_states()` classifies
  regions into biologically meaningful states (active, bivalent, poised, primed,
  repressed, unmarked) following ChromHMM/Roadmap Epigenomics conventions.

* **TF co-binding analysis**: `epiRomics_tf_cobinding()` uses Fisher's exact
  test with odds ratios and pointwise mutual information (PMI) to identify
  significant TF co-binding pairs.

## Speed Optimizations

* Parallel overlap counting via `parallel::mclapply()` when available
* Pre-split annotations by type for O(1) lookup in `epiRomics_enhanceosome()`
* Vectorized chromatin state classification and signal aggregation
* BigWig RDS caching for memory-efficient large-scale analyses

## CI/CD

* Cross-platform GitHub Actions matrix: Ubuntu (release/devel/oldrel-1),
  macOS ARM, macOS Intel, Windows
* BiocCheck validation on every push/PR
* Test coverage with full integration tests via cached example data

## Documentation

* Comprehensive HTML and PDF vignette with live BigWig signal visualizations
* pkgdown site deployment via GitHub Pages
* `.onAttach()` startup message with version, cache status, and citation info
