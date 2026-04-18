# ============================================================================
# REGRESSION BASELINE TESTS for classify_chromatin_states
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
# Test 1: chromatin_states regression: default params (refine_by_tss = FALSE)
# ============================================================================
test_that("chromatin_states regression: default params", {
  dB <- make_synthetic_dB_histone_only()
  result <- classify_chromatin_states(dB, refine_by_tss = FALSE)
  expect_regression_match(result, "chromatin_states_default")
})

# ============================================================================
# Test 2: chromatin_states regression: with TSS refinement
# ============================================================================
test_that("chromatin_states regression: with TSS refinement", {
  dB <- make_synthetic_dB_histone_only()
  result <- classify_chromatin_states(dB, refine_by_tss = TRUE)
  expect_regression_match(result, "chromatin_states_tss")
})

# ============================================================================
# Test 3: chromatin_states regression: custom regions
# ============================================================================
test_that("chromatin_states regression: custom regions", {
  dB <- make_synthetic_dB_histone_only()
  custom_regions <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = base::c(900L, 4500L, 9500L),
      end = base::c(2100L, 6500L, 11500L)
    )
  )
  result <- classify_chromatin_states(dB, regions = custom_regions,
                                        refine_by_tss = FALSE)
  expect_regression_match(result, "chromatin_states_custom_regions")
})
