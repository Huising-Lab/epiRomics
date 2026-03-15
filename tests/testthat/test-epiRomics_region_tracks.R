library(testthat)
library(epiRomics)

# Setup test data using helpers
epiRomics_dB <- build_test_dB()
track_connection <- make_track_connection()

test_region <- GenomicRanges::GRanges(
  seqnames = "chr11",
  ranges = IRanges::IRanges(start = 2150000, end = 2155000),
  strand = "+",
  SYMBOL = "INS"
)

# Test 1: Function runs without error
test_that("epiRomics_region_tracks runs without error", {
  skip_if_no_extdata()
  # The function produces a plot and returns invisible NULL
  expect_no_error(
    epiRomics::epiRomics_region_tracks(
      epiRomics_region = test_region,
      epiRomics_dB = epiRomics_dB,
      epiRomics_track_connection = track_connection
    )
  )
})

# Test 2: Function performance is acceptable
test_that("epiRomics_region_tracks performance is acceptable", {
  skip_if_no_extdata()
  start_time <- Sys.time()
  epiRomics::epiRomics_region_tracks(
    epiRomics_region = test_region,
    epiRomics_dB = epiRomics_dB,
    epiRomics_track_connection = track_connection
  )
  end_time <- Sys.time()
  execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  expect_true(execution_time < 30)
})

# Test 3: Function validates parameters correctly
test_that("epiRomics_region_tracks validates parameters correctly", {
  skip_if_no_extdata()
  # Test with invalid region (NULL)
  expect_error(
    epiRomics::epiRomics_region_tracks(
      epiRomics_region = NULL,
      epiRomics_dB = epiRomics_dB,
      epiRomics_track_connection = track_connection
    )
  )

  # Test with invalid database (string)
  expect_error(
    epiRomics::epiRomics_region_tracks(
      epiRomics_region = test_region,
      epiRomics_dB = "invalid",
      epiRomics_track_connection = track_connection
    )
  )
})
