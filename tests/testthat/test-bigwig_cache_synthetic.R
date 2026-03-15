# ============================================================================
# SYNTHETIC DATA TESTS for BigWig caching utilities
# Tests: cache_bigwig (mock), load_cached_bigwig, maxCovBwCached logic
# ============================================================================

library(testthat)
library(epiRomics)
library(GenomicRanges)
library(IRanges)

# ============================================================================
# load_cached_bigwig tests (with synthetic RDS)
# ============================================================================

test_that("load_cached_bigwig loads full GRanges from RDS", {
  bw_data <- GRanges("chr1", IRanges(c(100, 500), c(200, 600)),
    score = c(10, 20))
  cache_file <- tempfile(fileext = ".rds")
  saveRDS(bw_data, cache_file)
  on.exit(unlink(cache_file))

  loaded <- epiRomics:::load_cached_bigwig(cache_file)
  expect_s4_class(loaded, "GRanges")
  expect_equal(length(loaded), 2)
  expect_equal(loaded$score, c(10, 20))
})

test_that("load_cached_bigwig subsets by gr when provided", {
  bw_data <- GRanges("chr1", IRanges(c(100, 500, 1000), c(200, 600, 1100)),
    score = c(10, 20, 30))
  cache_file <- tempfile(fileext = ".rds")
  saveRDS(bw_data, cache_file)
  on.exit(unlink(cache_file))

  # Only query the first two regions
  gr <- GRanges("chr1", IRanges(50, 700))
  loaded <- epiRomics:::load_cached_bigwig(cache_file, gr = gr)
  expect_equal(length(loaded), 2)
  expect_equal(loaded$score, c(10, 20))
})

test_that("load_cached_bigwig errors on missing cache file", {
  expect_error(
    epiRomics:::load_cached_bigwig("/nonexistent/cache.rds"),
    "do not exist")
})

test_that("load_cached_bigwig returns empty for non-overlapping gr", {
  bw_data <- GRanges("chr1", IRanges(100, 200), score = 10)
  cache_file <- tempfile(fileext = ".rds")
  saveRDS(bw_data, cache_file)
  on.exit(unlink(cache_file))

  gr <- GRanges("chr1", IRanges(500, 600))
  loaded <- epiRomics:::load_cached_bigwig(cache_file, gr = gr)
  expect_equal(length(loaded), 0)
})

# ============================================================================
# cache_bigwig validation tests (without actual BigWig)
# ============================================================================

test_that("cache_bigwig errors on non-existent BigWig path", {
  expect_error(
    epiRomics:::cache_bigwig("/nonexistent/file.bw"),
    "do not exist")
})

# ============================================================================
# maxCovBwCached validation tests
# ============================================================================

test_that("maxCovBwCached errors on non-existent BigWig", {
  gr <- GRanges("chr1", IRanges(100, 200))
  expect_error(
    maxCovBwCached("/nonexistent/file.bw", gr),
    "do not exist")
})

test_that("maxCovBwCached errors on invalid GRanges", {
  # It will error on file validation first; test the gr validation separately
  # by using validate_genomic_ranges directly
  expect_error(
    epiRomics:::validate_genomic_ranges("not_granges"),
    "must be a GRanges object")
})
