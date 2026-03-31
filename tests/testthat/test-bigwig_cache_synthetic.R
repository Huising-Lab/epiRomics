# ============================================================
# SYNTHETIC DATA TESTS for BigWig coverage utilities
# Tests: maxCovBwCached validation
# ============================================================

library(testthat)
library(epiRomics)
library(GenomicRanges)
library(IRanges)

# ============================================================
# maxCovBwCached validation tests
# ============================================================

test_that("maxCovBwCached errors on non-existent BigWig", {
  gr <- GRanges("chr1", IRanges(100, 200))
  expect_error(
    maxCovBwCached("/nonexistent/file.bw", gr),
    "do not exist")
})

test_that("maxCovBwCached errors on invalid GRanges", {
  expect_error(
    epiRomics:::validate_genomic_ranges("not_granges"),
    "must be a GRanges object")
})
