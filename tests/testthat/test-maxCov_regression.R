# ============================================================================
# REGRESSION BASELINE TESTS for maxCovBwCached and maxCovFilesCached
# Uses synthetic BigWig files for deterministic testing.
# Baselines saved as .rds in testdata/ via expect_regression_match().
# ============================================================================

library(testthat)
library(epiRomics)
library(GenomicRanges)
library(IRanges)

test_that("maxCovBwCached regression: single BigWig", {
  testthat::skip_if_not_installed("rtracklayer")

  # Create synthetic BigWig with known scores
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = base::seq(1L, 9001L, 1000L),
      end = base::seq(1000L, 10000L, 1000L)
    )
  )
  scores <- base::c(5, 10, 15, 80, 90, 5, 10, 95, 5, 10)

  bw_path <- make_synthetic_bigwig(gr, scores)
  base::on.exit(base::file.remove(bw_path), add = TRUE)

  # Query region covering all bins
  query_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1L, end = 10000L)
  )

  result <- maxCovBwCached(bw_path, query_gr)

  expect_regression_match(result, "maxCovBwCached_default", tolerance = 1e-6)
})

test_that("maxCovFilesCached regression: multiple BigWigs", {
  testthat::skip_if_not_installed("rtracklayer")

  # Create two synthetic BigWig files with different max scores
  gr1 <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = base::seq(1L, 4001L, 1000L),
      end = base::seq(1000L, 5000L, 1000L)
    )
  )
  scores1 <- base::c(10, 20, 30, 40, 50)

  gr2 <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = base::seq(1L, 4001L, 1000L),
      end = base::seq(1000L, 5000L, 1000L)
    )
  )
  scores2 <- base::c(60, 70, 80, 90, 100)

  bw1 <- make_synthetic_bigwig(gr1, scores1)
  bw2 <- make_synthetic_bigwig(gr2, scores2)
  base::on.exit({
    base::file.remove(bw1)
    base::file.remove(bw2)
  }, add = TRUE)

  # Query region covering all bins
  query_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1L, end = 5000L)
  )

  result <- maxCovFilesCached(base::c(bw1, bw2), query_gr)

  # Convert GRanges result to data.frame for stable serialization
  result_df <- base::as.data.frame(result)

  expect_regression_match(result_df, "maxCovFilesCached_default", tolerance = 1e-6)
})
