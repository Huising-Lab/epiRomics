# ============================================================================
# REGRESSION BASELINE TESTS for call_accessible_regions
# Uses synthetic BigWig files for deterministic testing.
# Baselines saved as .rds in testdata/ via expect_regression_match().
# ============================================================================

library(testthat)
library(epiRomics)
library(GenomicRanges)
library(IRanges)

test_that("call_accessible_regions regression: z-score method", {
  # rtracklayer is an Imports of epiRomics; skip gate removed.


  # Create synthetic BigWig with known score distribution:
  # scores 5,10,15,80,90,5,10,95,5,10 -- some high, some low
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

  # Query regions covering the full extent
  regions <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = base::seq(1L, 9001L, 1000L),
      end = base::seq(1000L, 10000L, 1000L)
    )
  )

  result <- call_accessible_regions(
    bw_path, regions,
    z_threshold = 1.0,
    return_scores = TRUE
  )

  # Save the full result list as the regression baseline
  expect_regression_match(result, "call_accessible_default", tolerance = 1e-6)
})

test_that("call_accessible_regions regression: auto threshold", {
  # rtracklayer is an Imports of epiRomics; skip gate removed.

  # Same synthetic BigWig as above

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

  regions <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = base::seq(1L, 9001L, 1000L),
      end = base::seq(1000L, 10000L, 1000L)
    )
  )

  result <- call_accessible_regions(
    bw_path, regions,
    auto_threshold = TRUE,
    return_scores = TRUE
  )

  expect_regression_match(result, "call_accessible_auto", tolerance = 1e-6)
})
