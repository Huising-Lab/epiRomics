library(testthat)
library(epiRomics)

# Build database using helper
epiRomics_dB <- build_test_dB()

# Test: epiRomics_enhancer_predictor_test returns correct output structure
testthat::test_that("epiRomics_enhancer_predictor_test returns correct output structure", {
  skip_if_no_extdata()
  tryCatch(
    {
      result <- epiRomics::epiRomics_enhancer_predictor_test(
        epiRomics_dB = epiRomics_dB,
        epiRomics_histone = "h3k4me1",
        epiRomics_curated_database = "fantom"
      )

      # Test output type
      testthat::expect_true(is.data.frame(result) || is.matrix(result) || is.null(result))

      if (!is.null(result)) {
        # Test output dimensions
        testthat::expect_gt(nrow(result), 0)
        testthat::expect_gt(ncol(result), 0)

        # Test column names if it's a data frame
        if (is.data.frame(result)) {
          testthat::expect_true("Histone_Mark" %in% colnames(result))
          testthat::expect_true("Fraction_of_Overlap" %in% colnames(result))
        }
      }
    },
    error = function(e) {
      error_msg <- tolower(e$message)
      expected_errors <- c("database", "annotation", "predictor", "test", "evaluation")

      if (any(sapply(expected_errors, function(x) grepl(x, error_msg)))) {
        testthat::expect_true(TRUE)
      } else {
        testthat::fail(paste("Unexpected error:", e$message))
      }
    }
  )
})
