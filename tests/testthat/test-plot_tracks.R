library(testthat)
library(epiRomics)

# Shared fixture for the plot_tracks test suite. require_extdata() inside each
# test_that() block makes a missing tarball fail loudly; building the fixture
# up-front keeps each test focused on plot_tracks behavior.
.make_plot_tracks_fixture <- function() {
  database <- build_test_dB()
  track_connection <- make_track_connection()
  histone_marks <- database@meta$name[database@meta$type == "histone"]
  epiRomics_enhancers <- epiRomics::find_enhancers_by_comarks(
    database, histone_marks[1], histone_marks[2])
  putative_enhanceosome <- epiRomics::find_enhanceosomes(
    epiRomics_enhancers, database)
  base::list(
    database = database,
    track_connection = track_connection,
    putative_enhanceosome = putative_enhanceosome
  )
}

# Test: Function works with valid input
testthat::test_that("plot_tracks works with valid input", {
  require_extdata()
  fx <- .make_plot_tracks_fixture()
  result <- epiRomics::plot_tracks(
    putative_enhanceosome = fx$putative_enhanceosome,
    index = 1,
    database = fx$database,
    track_connection = fx$track_connection,
    keep_epitracks = TRUE
  )
  # Check that function runs without errors
  # Gviz::plotTracks can return different types, so we just check it's not an error
  testthat::expect_true(!inherits(result, "try-error"))
})

# Test: Function works with epitracks disabled
testthat::test_that("plot_tracks works with epitracks disabled", {
  require_extdata()
  fx <- .make_plot_tracks_fixture()
  result <- epiRomics::plot_tracks(
    putative_enhanceosome = fx$putative_enhanceosome,
    index = 1,
    database = fx$database,
    track_connection = fx$track_connection,
    keep_epitracks = FALSE
  )
  # Check that function runs without errors
  testthat::expect_true(!inherits(result, "try-error"))
})

# Test: Function performance is acceptable
testthat::test_that("plot_tracks performance is acceptable", {
  require_extdata()
  fx <- .make_plot_tracks_fixture()
  start_time <- Sys.time()
  result <- epiRomics::plot_tracks(
    putative_enhanceosome = fx$putative_enhanceosome,
    index = 1,
    database = fx$database,
    track_connection = fx$track_connection,
    keep_epitracks = TRUE
  )
  end_time <- Sys.time()
  execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Performance threshold: should complete within 30 seconds
  testthat::expect_true(execution_time < 30)
  testthat::expect_true(!inherits(result, "try-error"))
})

# Test 4: Function validates parameters correctly
testthat::test_that("plot_tracks validates parameters correctly", {
  # Parameter validation does not require toy extdata — the invalid-input
  # branches reject before any file I/O. Keep a minimal fixture so the
  # valid-enhanceosome/database references exist for the negative cases.
  require_extdata()
  fx <- .make_plot_tracks_fixture()
  putative_enhanceosome <- fx$putative_enhanceosome
  database <- fx$database
  track_connection <- fx$track_connection
  # Test with invalid index
  testthat::expect_error(
    epiRomics::plot_tracks(
      putative_enhanceosome = putative_enhanceosome,
      index = 999999, # Invalid index
      database = database,
      track_connection = track_connection
    )
  )

  # Test with invalid putative enhanceosome
  testthat::expect_error(
    epiRomics::plot_tracks(
      putative_enhanceosome = NULL,
      index = 1,
      database = database,
      track_connection = track_connection
    )
  )

  # Test with invalid database
  testthat::expect_error(
    epiRomics::plot_tracks(
      putative_enhanceosome = putative_enhanceosome,
      index = 1,
      database = NULL,
      track_connection = track_connection
    )
  )
})
