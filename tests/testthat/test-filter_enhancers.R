library(testthat)
library(epiRomics)

# Test: filter_enhancers returns correct S4 class and filtered output
testthat::test_that("filter_enhancers returns correct S4 class and filtered output", {
  require_extdata()
  database <- build_test_dB()
  enhancers <- epiRomics::find_enhancers_by_comarks(
    database, "h3k4me1", "h3k27ac")
  enhanceosome <- epiRomics::find_enhanceosomes(enhancers, database)
  result <- epiRomics::filter_enhancers(
    putative_enhancers = enhanceosome,
    database = database,
    type = "hg38_custom_fantom"
  )
  testthat::expect_s4_class(result, "epiRomicsS4")
  testthat::expect_true(length(result@annotations) >= 0)
  testthat::expect_lte(length(result@annotations), length(enhanceosome@annotations))
})
