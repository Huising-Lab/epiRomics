# ============================================================================
# REGRESSION BASELINE TESTS for epiRomics_enhancer_predictor_to_ref
#
# Per D-04/D-05/D-06: Run function on synthetic data, save output as .rds,
# compare on future runs. Any code change that alters output will be caught.
#
# Uses shared builders from helper-synthetic-builders.R:
#   - make_synthetic_dB_full() for database with histones + ChIP TFs
#   - expect_regression_match() for .rds save/load/compare
#
# epiRomics_enhancer_predictor_to_ref computes overlap fractions between
# histone marks and a curated reference database. The curated_database
# parameter must name an entry whose meta$type is "histone" or "chip".
# In make_synthetic_dB_full(), TF1 and TF2 are type "chip", so they
# can serve as curated reference databases for regression testing.
#
# Returns data.frame directly -- no S4 extraction needed.
# ============================================================================

library(testthat)
library(epiRomics)
library(GenomicRanges)
library(IRanges)

# ============================================================================
# Test 1: enhancer_predictor_to_ref regression: default (h3k4me1 vs TF1)
# ============================================================================
test_that("enhancer_predictor_to_ref regression: default", {
  dB <- make_synthetic_dB_full()
  result <- epiRomics_enhancer_predictor_to_ref(dB,
    epiRomics_histone = "h3k4me1",
    epiRomics_curated_database = "TF1")
  expect_regression_match(result, "enhancer_predictor_to_ref_default")
})
