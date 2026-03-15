library(testthat)
library(epiRomics)

# Build database using helper
epiRomics_dB <- build_test_dB()

# Get putative enhancers and enhanceosome
enhancers <- epiRomics::epiRomics_enhancers(epiRomics_dB, "h3k4me1", "h3k27ac")
enhanceosome <- epiRomics::epiRomics_enhanceosome(enhancers, epiRomics_dB)

# Test: epiRomics_enhancers_filter returns correct S4 class and filtered output

testthat::test_that("epiRomics_enhancers_filter returns correct S4 class and filtered output", {
  skip_if_no_extdata()
  result <- epiRomics::epiRomics_enhancers_filter(
    epiRomics_putative_enhancers = enhanceosome,
    epiRomics_dB = epiRomics_dB,
    epiRomics_type = "hg38_custom_fantom"
  )
  testthat::expect_s4_class(result, "epiRomicsS4")
  testthat::expect_true(length(result@annotations) >= 0)
  testthat::expect_lte(length(result@annotations), length(enhanceosome@annotations))
})
