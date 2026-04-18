library(testthat)
library(epiRomics)

# Test: find_enhanceosomes returns correct S4 class and expected output
testthat::test_that("find_enhanceosomes returns correct S4 class and expected output", {
  require_extdata()
  database <- build_test_dB()
  histone_marks <- database@meta$name[database@meta$type == "histone"]
  epiRomics_enhancers <- epiRomics::find_enhancers_by_comarks(
    database, histone_marks[1], histone_marks[2])
  result <- epiRomics::find_enhanceosomes(epiRomics_enhancers, database)
  testthat::expect_s4_class(result, "epiRomicsS4")
  testthat::expect_true(length(result@annotations) >= 0)
  testthat::expect_lte(length(result@annotations), length(epiRomics_enhancers@annotations))
})
