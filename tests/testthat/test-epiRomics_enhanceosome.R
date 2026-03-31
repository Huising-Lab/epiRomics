library(testthat)
library(epiRomics)

# Build database using helper
epiRomics_dB <- build_test_dB()

# Generate putative enhancers using real data (guard against NULL when extdata missing)
epiRomics_enhancers <- NULL
if (!is.null(epiRomics_dB)) {
  histone_marks <- epiRomics_dB@meta$name[epiRomics_dB@meta$type == "histone"]
  epiRomics_enhancers <- epiRomics::epiRomics_enhancers_co_marks(epiRomics_dB, histone_marks[1], histone_marks[2])
}

# Test: epiRomics_enhanceosome returns correct S4 class and expected output
testthat::test_that("epiRomics_enhanceosome returns correct S4 class and expected output", {
  skip_if_no_extdata()
  result <- epiRomics::epiRomics_enhanceosome(epiRomics_enhancers, epiRomics_dB)
  testthat::expect_s4_class(result, "epiRomicsS4")
  testthat::expect_true(length(result@annotations) >= 0)
  testthat::expect_lte(length(result@annotations), length(epiRomics_enhancers@annotations))
})
