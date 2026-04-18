# Renamed function aliases.
#
# Every exported function was renamed from the `epiRomics_*` prefix
# convention to verb-first, descriptive names. These wrappers keep the
# prior API working for one release cycle. Each wrapper emits a
# message() notice (once per session, per old name) directing users to
# the new name, then forwards all arguments through `...`.
#
# `.Deprecated()` is intentionally not used here: this is a first-time
# Bioconductor submission whose user base predates Bioc (published on
# GitHub since 2021 with citations), and `.Deprecated()` would trigger
# BiocCheck without improving the user experience beyond the message()
# notice below. All wrappers defined here will be removed in the next
# release.

.warn_renamed <- local({
  state <- new.env(parent = emptyenv())
  state$seen <- character(0)
  function(old, new) {
    if (!(old %in% state$seen)) {
      message("'", old, "' has been renamed to '", new, "'. ",
              "Please use '", new, "' in new code. ",
              "'", old, "' will be removed in the next release of epiRomics.")
      state$seen <- c(state$seen, old)
    }
  }
})

#' Renamed epiRomics functions
#'
#' These functions have been renamed in epiRomics 0.99.1 to follow
#' Bioconductor naming conventions. They remain as thin wrappers that
#' emit a one-time rename notice and forward to the new name. They
#' will be removed in the next release.
#'
#' @param ... Arguments forwarded to the replacement function.
#' @name epiRomics-deprecated
#' @keywords internal
NULL

#' @rdname epiRomics-deprecated
#' @export
epiRomics_build_dB <- function(...) {
  .warn_renamed("epiRomics_build_dB", "build_database")
  build_database(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_cache_data <- function(...) {
  .warn_renamed("epiRomics_cache_data", "cache_data")
  cache_data(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_cache_path <- function(...) {
  .warn_renamed("epiRomics_cache_path", "get_cache_path")
  get_cache_path(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_has_cache <- function(...) {
  .warn_renamed("epiRomics_has_cache", "has_cache")
  has_cache(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_chromatin_states <- function(...) {
  .warn_renamed("epiRomics_chromatin_states", "classify_chromatin_states")
  classify_chromatin_states(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_chromatin_states_categories <- function(...) {
  .warn_renamed("epiRomics_chromatin_states_categories",
                "chromatin_state_categories")
  chromatin_state_categories(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_enhanceosome <- function(...) {
  .warn_renamed("epiRomics_enhanceosome", "find_enhanceosomes")
  find_enhanceosomes(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_enhancer_predictor_to_ref <- function(...) {
  .warn_renamed("epiRomics_enhancer_predictor_to_ref",
                "benchmark_enhancer_predictor")
  benchmark_enhancer_predictor(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_enhancers_co_marks <- function(...) {
  .warn_renamed("epiRomics_enhancers_co_marks", "find_enhancers_by_comarks")
  find_enhancers_by_comarks(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_enhancers_filter <- function(...) {
  .warn_renamed("epiRomics_enhancers_filter", "filter_enhancers")
  filter_enhancers(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_filter_accessible <- function(...) {
  .warn_renamed("epiRomics_filter_accessible", "filter_accessible_regions")
  filter_accessible_regions(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_annotate_putative <- function(...) {
  .warn_renamed("epiRomics_annotate_putative", "annotate_enhancers")
  annotate_enhancers(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_putative_enhancers <- function(...) {
  .warn_renamed("epiRomics_putative_enhancers", "find_putative_enhancers")
  find_putative_enhancers(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_quick_view <- function(...) {
  .warn_renamed("epiRomics_quick_view", "plot_quick_view")
  plot_quick_view(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_regions_of_interest <- function(...) {
  .warn_renamed("epiRomics_regions_of_interest", "get_regions_of_interest")
  get_regions_of_interest(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_tf_cobinding <- function(...) {
  .warn_renamed("epiRomics_tf_cobinding", "analyze_tf_cobinding")
  analyze_tf_cobinding(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_tf_overlap <- function(...) {
  .warn_renamed("epiRomics_tf_overlap", "analyze_tf_overlap")
  analyze_tf_overlap(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_track_layer <- function(...) {
  .warn_renamed("epiRomics_track_layer", "plot_tracks")
  plot_tracks(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_track_layer_fast <- function(...) {
  .warn_renamed("epiRomics_track_layer_fast", "plot_tracks_fast")
  plot_tracks_fast(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_track_layer_gene <- function(...) {
  .warn_renamed("epiRomics_track_layer_gene", "plot_gene_tracks")
  plot_gene_tracks(...)
}
