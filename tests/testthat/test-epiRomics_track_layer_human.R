library(testthat)
library(epiRomics)

# Build database and track connection using helpers
epiRomics_dB <- build_test_dB()
track_connection <- make_track_connection()

# Generate putative enhanceosome using real data
histone_marks <- epiRomics_dB@meta$name[epiRomics_dB@meta$type == "histone"]
tf_marks <- epiRomics_dB@meta$name[epiRomics_dB@meta$type == "chip"]
# Use first two histone marks for enhancers
epiRomics_enhancers <- epiRomics::epiRomics_enhancers(epiRomics_dB, histone_marks[1], histone_marks[2])
putative_enhanceosome <- epiRomics::epiRomics_enhanceosome(epiRomics_enhancers, epiRomics_dB)

# Test: Function works with valid input
testthat::test_that("epiRomics_track_layer_human works with valid input", {
  skip_if_no_extdata()
  result <- epiRomics::epiRomics_track_layer_human(
    epiRomics_putative_enhanceosome = putative_enhanceosome,
    epiRomics_index = 1,
    epiRomics_dB = epiRomics_dB,
    epiRomics_track_connection = track_connection,
    epiRomics_keep_epitracks = TRUE
  )
  # Check that function runs without errors
  testthat::expect_true(!inherits(result, "try-error"))
})

# Test: Function works with epitracks disabled
testthat::test_that("epiRomics_track_layer_human works with epitracks disabled", {
  skip_if_no_extdata()
  result <- epiRomics::epiRomics_track_layer_human(
    epiRomics_putative_enhanceosome = putative_enhanceosome,
    epiRomics_index = 1,
    epiRomics_dB = epiRomics_dB,
    epiRomics_track_connection = track_connection,
    epiRomics_keep_epitracks = FALSE
  )
  # Check that function runs without errors
  testthat::expect_true(!inherits(result, "try-error"))
})

# Test: Function performance is acceptable
testthat::test_that("epiRomics_track_layer_human performance is acceptable", {
  skip_if_no_extdata()
  start_time <- Sys.time()
  result <- epiRomics::epiRomics_track_layer_human(
    epiRomics_putative_enhanceosome = putative_enhanceosome,
    epiRomics_index = 1,
    epiRomics_dB = epiRomics_dB,
    epiRomics_track_connection = track_connection,
    epiRomics_keep_epitracks = TRUE
  )
  end_time <- Sys.time()
  execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Performance threshold: should complete within 30 seconds
  testthat::expect_true(execution_time < 30)
  testthat::expect_true(!inherits(result, "try-error"))
})
