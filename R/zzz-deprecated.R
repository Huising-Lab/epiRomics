# Deprecated function aliases.
#
# Every exported function was renamed from the `epiRomics_*` prefix
# convention to verb-first, descriptive names. These wrappers keep the
# prior API working for one release cycle with a `.Deprecated()` notice
# pointing users at the new name. All arguments are forwarded through
# `...` so renamed parameters are also honored.
#
# All wrappers defined here will be removed in a future release.

#' Deprecated epiRomics functions
#'
#' These functions have been renamed in epiRomics 0.99.1 to follow
#' Bioconductor naming conventions. They remain as thin wrappers that
#' emit a deprecation warning and forward to the new name.
#'
#' @param ... Arguments forwarded to the replacement function.
#' @name epiRomics-deprecated
#' @keywords internal
NULL

#' @rdname epiRomics-deprecated
#' @export
epiRomics_build_dB <- function(...) {
  .Deprecated("build_database", old = "epiRomics_build_dB")
  build_database(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_cache_data <- function(...) {
  .Deprecated("cache_data", old = "epiRomics_cache_data")
  cache_data(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_cache_path <- function(...) {
  .Deprecated("get_cache_path", old = "epiRomics_cache_path")
  get_cache_path(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_has_cache <- function(...) {
  .Deprecated("has_cache", old = "epiRomics_has_cache")
  has_cache(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_chromatin_states <- function(...) {
  .Deprecated("classify_chromatin_states", old = "epiRomics_chromatin_states")
  classify_chromatin_states(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_chromatin_states_categories <- function(...) {
  .Deprecated("chromatin_state_categories",
              old = "epiRomics_chromatin_states_categories")
  chromatin_state_categories(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_enhanceosome <- function(...) {
  .Deprecated("find_enhanceosomes", old = "epiRomics_enhanceosome")
  find_enhanceosomes(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_enhancer_predictor_to_ref <- function(...) {
  .Deprecated("benchmark_enhancer_predictor",
              old = "epiRomics_enhancer_predictor_to_ref")
  benchmark_enhancer_predictor(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_enhancers_co_marks <- function(...) {
  .Deprecated("find_enhancers_by_comarks", old = "epiRomics_enhancers_co_marks")
  find_enhancers_by_comarks(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_enhancers_filter <- function(...) {
  .Deprecated("filter_enhancers", old = "epiRomics_enhancers_filter")
  filter_enhancers(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_filter_accessible <- function(...) {
  .Deprecated("filter_accessible_regions", old = "epiRomics_filter_accessible")
  filter_accessible_regions(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_annotate_putative <- function(...) {
  .Deprecated("annotate_enhancers", old = "epiRomics_annotate_putative")
  annotate_enhancers(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_putative_enhancers <- function(...) {
  .Deprecated("find_putative_enhancers", old = "epiRomics_putative_enhancers")
  find_putative_enhancers(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_quick_view <- function(...) {
  .Deprecated("plot_quick_view", old = "epiRomics_quick_view")
  plot_quick_view(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_regions_of_interest <- function(...) {
  .Deprecated("get_regions_of_interest", old = "epiRomics_regions_of_interest")
  get_regions_of_interest(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_tf_cobinding <- function(...) {
  .Deprecated("analyze_tf_cobinding", old = "epiRomics_tf_cobinding")
  analyze_tf_cobinding(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_tf_overlap <- function(...) {
  .Deprecated("analyze_tf_overlap", old = "epiRomics_tf_overlap")
  analyze_tf_overlap(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_track_layer <- function(...) {
  .Deprecated("plot_tracks", old = "epiRomics_track_layer")
  plot_tracks(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_track_layer_fast <- function(...) {
  .Deprecated("plot_tracks_fast", old = "epiRomics_track_layer_fast")
  plot_tracks_fast(...)
}

#' @rdname epiRomics-deprecated
#' @export
epiRomics_track_layer_gene <- function(...) {
  .Deprecated("plot_gene_tracks", old = "epiRomics_track_layer_gene")
  plot_gene_tracks(...)
}
