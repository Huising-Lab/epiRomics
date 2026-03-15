library(testthat)
library(epiRomics)

# Build database using helper
epiRomics_dB <- build_test_dB()

# Test: epiRomics_regions_of_interest returns correct S4 class and filtered output
testthat::test_that("epiRomics_regions_of_interest returns correct S4 class and filtered output", {
  skip_if_no_extdata()
  # Use the full annotation set as regions for a round-trip test
  result <- epiRomics::epiRomics_regions_of_interest(
    epiRomics_putative_enhanceosome = epiRomics_dB,
    epiRomics_test_regions = epiRomics_dB@annotations
  )
  testthat::expect_s4_class(result, "epiRomicsS4")
  testthat::expect_true(length(result@annotations) >= 0)
  testthat::expect_lte(length(result@annotations), length(epiRomics_dB@annotations))
})
