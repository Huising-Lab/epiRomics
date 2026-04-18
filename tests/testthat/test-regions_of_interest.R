library(testthat)
library(epiRomics)

# Test: get_regions_of_interest returns correct S4 class and filtered output
testthat::test_that("get_regions_of_interest returns correct S4 class and filtered output", {
  require_extdata()
  database <- build_test_dB()
  # Use the full annotation set as regions for a round-trip test
  result <- epiRomics::get_regions_of_interest(
    putative_enhanceosome = database,
    test_regions = database@annotations
  )
  testthat::expect_s4_class(result, "epiRomicsS4")
  testthat::expect_true(length(result@annotations) >= 0)
  testthat::expect_lte(length(result@annotations), length(database@annotations))
})
