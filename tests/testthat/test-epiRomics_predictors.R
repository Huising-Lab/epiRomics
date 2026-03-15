# Test file for epiRomics_predictors function (deprecated)
# The function now delegates to epiRomics_tf_cobinding and returns NULL
# with a deprecation warning.

library(testthat)
library(epiRomics)

# Setup test data
test_that("epiRomics_predictors setup", {
  enhanceosome_path <- "expected_outputs/epiRomics_putative_enhanceosome_fantom.rds"
  skip_if(!file.exists(enhanceosome_path), "Expected output not available")
  enhanceosome <<- readRDS(enhanceosome_path)

  expect_true(inherits(enhanceosome, "epiRomicsS4"))
  expect_true(length(enhanceosome@annotations) > 0)
  expect_true(exists("epiRomics_predictors"))
})

# Test 1: Function structure
test_that("epiRomics_predictors has correct function structure", {
  expect_true(is.function(epiRomics_predictors))
  expect_equal(length(formals(epiRomics_predictors)), 1)
  expect_equal(names(formals(epiRomics_predictors)),
    "epiRomics_putative_enhanceosome")
})

# Test 2: Function is deprecated and returns NULL
test_that("epiRomics_predictors returns NULL with deprecation warning", {
  enhanceosome_path <- "expected_outputs/epiRomics_putative_enhanceosome_fantom.rds"
  skip_if(!file.exists(enhanceosome_path), "Expected output not available")

  expect_warning(
    result <- epiRomics_predictors(
      epiRomics_putative_enhanceosome = enhanceosome),
    "deprecated"
  )
  expect_null(result)
})

# Test 3: Function handles invalid input (deprecated, returns NULL with warnings)
test_that("epiRomics_predictors handles invalid input with warning", {
  # The deprecated function issues warnings but may not error
  # on all invalid inputs since it catches errors internally
  expect_warning(
    result <- epiRomics_predictors(
      epiRomics_putative_enhanceosome = "not_an_epiRomics_object")
  )
  expect_null(result)
})
