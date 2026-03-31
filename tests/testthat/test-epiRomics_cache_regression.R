# ============================================================================
# REGRESSION BASELINE TESTS for epiRomics cache function contracts
# Tests function signatures (formals) as regression baselines.
# Any refactoring that changes the function signature will be caught.
# Baselines saved as .rds in testdata/ via expect_regression_match().
# ============================================================================

library(testthat)
library(epiRomics)

test_that("epiRomics_has_cache regression: contract", {
  contract <- base::list(
    formals = base::as.list(base::formals(epiRomics_has_cache)),
    returns_logical = base::is.logical(epiRomics_has_cache())
  )

  expect_regression_match(contract, "has_cache_contract")
})

test_that("epiRomics_cache_path regression: contract", {
  contract <- base::list(
    formals = base::as.list(base::formals(epiRomics_cache_path)),
    return_class = base::class(epiRomics_cache_path())
  )

  expect_regression_match(contract, "cache_path_contract")
})

test_that("epiRomics_cache_data regression: contract", {
  # Do NOT call epiRomics_cache_data() -- it downloads from Zenodo.
  # Only test the function signature as a contract baseline.
  contract <- base::list(
    formals = base::as.list(base::formals(epiRomics_cache_data))
  )

  expect_regression_match(contract, "cache_data_contract")
})
