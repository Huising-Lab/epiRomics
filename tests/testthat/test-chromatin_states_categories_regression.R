# ============================================================================
# REGRESSION BASELINE TESTS for chromatin_state_categories
#
# Per D-04/D-05/D-06: Run function on synthetic data, save output as .rds,
# compare on future runs. Any code change that alters output will be caught.
#
# Uses shared builders from helper-synthetic-builders.R:
#   - make_synthetic_dB_histone_only() for pre-built 3-mark database
#   - expect_regression_match() for .rds save/load/compare
# ============================================================================

library(testthat)
library(epiRomics)
library(GenomicRanges)
library(IRanges)

# ============================================================================
# Test 1: chromatin_states_categories regression: default
# ============================================================================
test_that("chromatin_states_categories regression: default", {
  dB <- make_synthetic_dB_histone_only()
  result <- chromatin_state_categories(dB, refine_by_tss = FALSE)
  expect_regression_match(result, "chromatin_states_categories_default")
})
