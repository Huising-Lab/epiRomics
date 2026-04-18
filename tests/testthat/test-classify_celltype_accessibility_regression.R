# ============================================================================
# REGRESSION BASELINE TESTS for classify_celltype_accessibility
# Uses synthetic BigWig files for deterministic testing.
# Baselines saved as .rds in testdata/ via expect_regression_match().
# ============================================================================

library(testthat)
library(epiRomics)
library(GenomicRanges)
library(IRanges)

test_that("classify_celltype_accessibility regression: two celltypes", {
  # rtracklayer is an Imports of epiRomics; skip gate removed.

  # Alpha: high signal in first half, low in second
  gr_alpha <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = base::seq(1L, 9001L, 1000L),
      end = base::seq(1000L, 10000L, 1000L)
    )
  )
  scores_alpha <- base::c(90, 85, 80, 75, 70, 5, 5, 5, 5, 5)

  # Beta: low signal in first half, high in second
  gr_beta <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = base::seq(1L, 9001L, 1000L),
      end = base::seq(1000L, 10000L, 1000L)
    )
  )
  scores_beta <- base::c(5, 5, 5, 5, 5, 70, 75, 80, 85, 90)

  alpha_bw <- make_synthetic_bigwig(gr_alpha, scores_alpha)
  beta_bw <- make_synthetic_bigwig(gr_beta, scores_beta)
  base::on.exit({
    base::file.remove(alpha_bw)
    base::file.remove(beta_bw)
  }, add = TRUE)

  # Query regions covering all bins
  regions <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = base::seq(1L, 9001L, 1000L),
      end = base::seq(1000L, 10000L, 1000L)
    )
  )

  bw_paths <- base::c(Alpha = alpha_bw, Beta = beta_bw)

  result <- classify_celltype_accessibility(
    bw_paths = bw_paths,
    regions = regions,
    z_threshold = 1.0
  )

  # Remove the accessibility_matrix attribute for stable serialization
  # (matrix attributes can be fragile across R versions)
  result_clean <- result
  base::attr(result_clean, "accessibility_matrix") <- NULL

  expect_regression_match(result_clean, "classify_accessibility_default")
})
