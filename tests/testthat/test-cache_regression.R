# ============================================================================
# REGRESSION BASELINE TESTS for epiRomics cache function contracts
# Tests function signatures (formals) as regression baselines.
# Any refactoring that changes the function signature will be caught.
# Baselines saved as .rds in testdata/ via expect_regression_match().
# ============================================================================

library(testthat)
library(epiRomics)

test_that("has_cache regression: contract", {
  contract <- base::list(
    formals = base::as.list(base::formals(has_cache)),
    returns_logical = base::is.logical(has_cache())
  )

  expect_regression_match(contract, "has_cache_contract")
})

test_that("get_cache_path regression: contract", {
  # Only assert the function signature (formals) — the return class is
  # environment-dependent: "NULL" when no cache, "character" when cached.
  contract <- base::list(
    formals = base::as.list(base::formals(get_cache_path))
  )

  expect_regression_match(contract, "cache_path_contract")

  # Separately assert the runtime return class is one of the two valid
  # options, independent of baseline (CI may or may not have cache).
  rc <- base::class(get_cache_path())
  expect_true(rc %in% c("NULL", "character"))
})

test_that("cache_data regression: contract", {
  # Do NOT call cache_data() -- it downloads from Zenodo.
  # Only test the function signature as a contract baseline.
  contract <- base::list(
    formals = base::as.list(base::formals(cache_data))
  )

  expect_regression_match(contract, "cache_data_contract")
})
