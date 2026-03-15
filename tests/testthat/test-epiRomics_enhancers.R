library(testthat)
library(epiRomics)

# Build database using helper
epiRomics_dB <- build_test_dB()

# Get histone marks from meta
ehistone_marks <- epiRomics_dB@meta$name[epiRomics_dB@meta$type == "histone"]

# Test: epiRomics_enhancers returns correct S4 class and expected output
testthat::test_that("epiRomics_enhancers returns correct S4 class and expected output", {
  skip_if_no_extdata()
  result <- epiRomics::epiRomics_enhancers(epiRomics_dB, ehistone_marks[1], ehistone_marks[2])
  testthat::expect_s4_class(result, "epiRomicsS4")
  testthat::expect_true(length(result@annotations) >= 0)
  testthat::expect_lte(length(result@annotations), length(epiRomics_dB@annotations))
})
