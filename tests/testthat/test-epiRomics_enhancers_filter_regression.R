# ============================================================================
# REGRESSION BASELINE TESTS for epiRomics_enhancers_filter
#
# Per D-04/D-05/D-06: Run function on synthetic data, save output as .rds,
# compare on future runs. Any code change that alters output will be caught.
#
# Uses shared builders from helper-synthetic-builders.R:
#   - make_synthetic_dB_full() for database with histones + functional
#   - expect_regression_match() for .rds save/load/compare
#
# epiRomics_enhancers_filter requires:
#   1. A "putative enhancers" epiRomicsS4 (histone co-mark intersections)
#   2. A full epiRomicsS4 containing functional annotations (e.g., fantom)
#   3. A filter type matching the functional annotation naming convention
#
# Note: Returns epiRomicsS4. We extract @annotations as data.frame.
# ============================================================================

library(testthat)
library(epiRomics)
library(GenomicRanges)
library(IRanges)

# ============================================================================
# Test 1: enhancers_filter regression: default filter via fantom
# ============================================================================
test_that("enhancers_filter regression: default filter", {
  dB <- make_synthetic_dB_full()
  # Create putative enhancers by co-mark intersection (h3k4me1 + h3k27ac)
  pe_dB <- epiRomics_enhancers_co_marks(dB)
  # Filter against fantom reference; genome is "hg38" so type is "hg38_custom_fantom"
  result <- epiRomics_enhancers_filter(pe_dB, dB,
    epiRomics_type = "hg38_custom_fantom")
  result_df <- base::as.data.frame(result@annotations)
  expect_regression_match(result_df, "enhancers_filter_default")
})
