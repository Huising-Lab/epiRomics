# ============================================================================
# REGRESSION BASELINE TESTS for find_enhancers_by_comarks
#
# Per D-04/D-05/D-06: Run function on synthetic data, save output as .rds,
# compare on future runs. Any code change that alters output will be caught.
#
# Uses shared builders from helper-synthetic-builders.R:
#   - make_synthetic_dB_histone_only() for pre-built histone-mark database
#   - expect_regression_match() for .rds save/load/compare
#
# Note: find_enhancers_by_comarks returns an epiRomicsS4 object.
# We extract @annotations as data.frame for stable regression comparison
# (S4 objects may have environment-dependent slots).
# ============================================================================

library(testthat)
library(epiRomics)
library(GenomicRanges)
library(IRanges)

# ============================================================================
# Test 1: enhancers_co_marks regression: default marks (h3k4me1 + h3k27ac)
# ============================================================================
test_that("enhancers_co_marks regression: default marks", {
  dB <- make_synthetic_dB_histone_only()
  result <- find_enhancers_by_comarks(dB)
  result_df <- base::as.data.frame(result@annotations)
  expect_regression_match(result_df, "enhancers_co_marks_default")
})

# ============================================================================
# Test 2: enhancers_co_marks regression: custom marks (h3k4me1 + h3k27me3)
# ============================================================================
test_that("enhancers_co_marks regression: custom marks", {
  dB <- make_synthetic_dB_histone_only()
  result <- find_enhancers_by_comarks(dB,
    histone_mark_1 = "h3k4me1",
    histone_mark_2 = "h3k27me3")
  result_df <- base::as.data.frame(result@annotations)
  expect_regression_match(result_df, "enhancers_co_marks_custom")
})
