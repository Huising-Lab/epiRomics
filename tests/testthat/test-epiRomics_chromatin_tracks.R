# Test file for epiRomics_chromatin_tracks function
# Tests the function in isolation with real data

# Load required annotation and txdb packages
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Setup test data using helpers
library(epiRomics)
epiRomics_dB <- build_test_dB()
track_connection <- make_track_connection()

# Test 1: Function structure
testthat::test_that("epiRomics_chromatin_tracks has correct function structure", {
  testthat::expect_true(is.function(epiRomics_chromatin_tracks))
  testthat::expect_equal(length(formals(epiRomics_chromatin_tracks)), 3)
  testthat::expect_equal(
    names(formals(epiRomics_chromatin_tracks)),
    c("epiRomics_gene_name", "epiRomics_dB", "epiRomics_track_connection")
  )
})

# Test 2: Function works with valid input
testthat::test_that("epiRomics_chromatin_tracks works with valid gene name", {
  skip_if_no_extdata()
  result <- epiRomics::epiRomics_chromatin_tracks(
    epiRomics_gene_name = "INS",
    epiRomics_dB = epiRomics_dB,
    epiRomics_track_connection = track_connection
  )
  print(result)
  testthat::expect_true(is.list(result))
  testthat::expect_true(all(c("chr11", "Axis", "INS", "Beta_ATAC", "Beta_RNA") %in% names(result)))
})

# Test 3: Function handles invalid input
testthat::test_that("epiRomics_chromatin_tracks handles invalid input", {
  skip_if_no_extdata()
  # Test with non-existent gene name
  testthat::expect_error(
    epiRomics_chromatin_tracks(
      epiRomics_gene_name = "NONEXISTENT_GENE",
      epiRomics_dB = epiRomics_dB,
      epiRomics_track_connection = track_connection
    )
  )

  # Test with invalid database object
  testthat::expect_error(
    epiRomics_chromatin_tracks(
      epiRomics_gene_name = "INS",
      epiRomics_dB = "not_a_database",
      epiRomics_track_connection = track_connection
    )
  )
})

# Test 4: Function performance
testthat::test_that("epiRomics_chromatin_tracks performance is acceptable", {
  skip_if_no_extdata()
  start_time <- Sys.time()
  result <- epiRomics::epiRomics_chromatin_tracks(
    epiRomics_gene_name = "INS",
    epiRomics_dB = epiRomics_dB,
    epiRomics_track_connection = track_connection
  )
  end_time <- Sys.time()
  execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  testthat::expect_true(execution_time < 30) # Should complete within 30 seconds
})

# Test 5: Function parameter validation
testthat::test_that("epiRomics_chromatin_tracks validates parameters correctly", {
  skip_if_no_extdata()
  # Test with NULL gene name
  testthat::expect_error(
    epiRomics_chromatin_tracks(
      epiRomics_gene_name = NULL,
      epiRomics_dB = epiRomics_dB,
      epiRomics_track_connection = track_connection
    )
  )

  # Test with empty track connection
  testthat::expect_error(
    epiRomics_chromatin_tracks(
      epiRomics_gene_name = "INS",
      epiRomics_dB = epiRomics_dB,
      epiRomics_track_connection = data.frame()
    )
  )
})
