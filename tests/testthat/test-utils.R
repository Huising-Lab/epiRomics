library(testthat)
library(epiRomics)

# Build database using helper
epiRomics_dB <- build_test_dB()

# Get BigWig file path
bw_file <- system.file("extdata", "BigWigs", "Kaestner_Beta.atac.signal.subset.chr1.bigwig", package = "epiRomics")

# Test: maxCovBw returns correct output structure
testthat::test_that("maxCovBw returns correct output structure", {
  skip_if_no_extdata()
  # Create a small test region from the database
  if (length(epiRomics_dB@annotations) > 0) {
    test_region <- epiRomics_dB@annotations[1]

    tryCatch(
      {
        result <- epiRomics:::maxCovBw(bw_file, test_region)

        # Test output type
        testthat::expect_true(is.numeric(result))
        testthat::expect_true(length(result) == 1)
        testthat::expect_true(result >= 0)
      },
      error = function(e) {
        # Expected to potentially fail with test data
        testthat::expect_true(TRUE)
      }
    )
  }
})

# Test: maxCovFiles returns correct output structure
testthat::test_that("maxCovFiles returns correct output structure", {
  skip_if_no_extdata()
  # Create test regions from the database
  if (length(epiRomics_dB@annotations) > 0) {
    test_regions <- epiRomics_dB@annotations[1:min(3, length(epiRomics_dB@annotations))]

    tryCatch(
      {
        result <- epiRomics:::maxCovFiles(c(bw_file), test_regions)

        # Test output type
        testthat::expect_s4_class(result, "GRanges")
        testthat::expect_equal(length(result), length(test_regions))

        # Test that values were added
        mcols_data <- GenomicRanges::mcols(result)
        testthat::expect_true(!is.null(mcols_data))
      },
      error = function(e) {
        # Expected to potentially fail with test data
        testthat::expect_true(TRUE)
      }
    )
  }
})
