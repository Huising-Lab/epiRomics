# Deprecation alias contract tests.
#
# Every old `epiRomics_*` exported function must:
#   (1) still be defined and exported,
#   (2) emit a one-time rename message when called (epiRomics 0.99.4
#       replaced `.Deprecated()` with a `message()` notice from the
#       shared `.warn_renamed()` helper),
#   (3) point the user at the correct replacement function.
#
# These tests call each alias with zero arguments inside
# `tryCatch()` so argument-validation errors thrown by the delegated
# new function do not fail the test — we only assert the rename
# message and its replacement target.

alias_map <- list(
  epiRomics_build_dB                   = "build_database",
  epiRomics_cache_data                 = "cache_data",
  epiRomics_cache_path                 = "get_cache_path",
  epiRomics_has_cache                  = "has_cache",
  epiRomics_chromatin_states           = "classify_chromatin_states",
  epiRomics_chromatin_states_categories = "chromatin_state_categories",
  epiRomics_enhanceosome               = "find_enhanceosomes",
  epiRomics_enhancer_predictor_to_ref  = "benchmark_enhancer_predictor",
  epiRomics_enhancers_co_marks         = "find_enhancers_by_comarks",
  epiRomics_enhancers_filter           = "filter_enhancers",
  epiRomics_filter_accessible          = "filter_accessible_regions",
  epiRomics_annotate_putative          = "annotate_enhancers",
  epiRomics_putative_enhancers         = "find_putative_enhancers",
  epiRomics_quick_view                 = "plot_quick_view",
  epiRomics_regions_of_interest        = "get_regions_of_interest",
  epiRomics_tf_cobinding               = "analyze_tf_cobinding",
  epiRomics_tf_overlap                 = "analyze_tf_overlap",
  epiRomics_track_layer                = "plot_tracks",
  epiRomics_track_layer_fast           = "plot_tracks_fast",
  epiRomics_track_layer_gene           = "plot_gene_tracks"
)

test_that("every deprecated alias is exported and resolves to a function", {
  for (old_name in names(alias_map)) {
    expect_true(
      exists(old_name, envir = asNamespace("epiRomics"), inherits = FALSE),
      info = paste("missing alias:", old_name)
    )
    fn <- get(old_name, envir = asNamespace("epiRomics"))
    expect_true(is.function(fn), info = paste("not a function:", old_name))
  }
})

test_that("each deprecated alias emits a rename message naming its replacement", {
  # The rename helper deduplicates per old-name within a single session,
  # so each alias may have been called earlier in the suite. Reset the
  # internal `state$seen` tracker before asserting so every alias fires
  # fresh exactly once here.
  deprecated_env <- environment(
    get(".warn_renamed", envir = asNamespace("epiRomics"))
  )
  if (exists("state", envir = deprecated_env)) {
    assign("seen", character(0), envir = deprecated_env$state)
    # keep getter_hint_shown TRUE so we do not assert on the hint line
    assign("getter_hint_shown", TRUE, envir = deprecated_env$state)
  }

  for (old_name in names(alias_map)) {
    new_name <- alias_map[[old_name]]
    fn <- get(old_name, envir = asNamespace("epiRomics"))

    messaged <- FALSE
    msg_text <- NA_character_
    tryCatch(
      withCallingHandlers(
        fn(),
        message = function(m) {
          # Capture only the FIRST message (the rename notice);
          # later messages from the delegated function must not overwrite it.
          if (!messaged) {
            messaged <<- TRUE
            msg_text <<- conditionMessage(m)
          }
          invokeRestart("muffleMessage")
        }
      ),
      error = function(e) invisible(NULL) # delegated fn may fail on empty args
    )

    expect_true(messaged,
                info = paste("alias did not emit message:", old_name))
    expect_match(msg_text, old_name, fixed = TRUE,
                 info = paste("message did not cite old name:", old_name))
    expect_match(msg_text, new_name, fixed = TRUE,
                 info = paste("message did not cite replacement:", old_name,
                              "->", new_name))
  }
})

test_that("all 20 aliases are present in the package NAMESPACE exports", {
  ns_exports <- getNamespaceExports("epiRomics")
  missing <- setdiff(names(alias_map), ns_exports)
  expect_length(missing, 0)
})
