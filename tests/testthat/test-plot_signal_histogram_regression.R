# ============================================================================
# REGRESSION BASELINE TESTS for plot_signal_histogram
# Uses synthetic BigWig files for deterministic testing.
# Baselines saved as .rds in testdata/ via expect_regression_match().
# ============================================================================

library(testthat)
library(epiRomics)
library(GenomicRanges)
library(IRanges)

test_that("plot_signal_histogram regression: default", {
  testthat::skip_if_not_installed("rtracklayer")

  # Create synthetic BigWig with known score distribution
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = base::seq(1L, 9001L, 1000L),
      end = base::seq(1000L, 10000L, 1000L)
    )
  )
  scores <- base::c(5, 10, 15, 80, 90, 5, 10, 95, 5, 10)

  bw_path <- make_synthetic_bigwig(gr, scores)
  pdf_path <- base::tempfile(fileext = ".pdf")
  base::on.exit({
    grDevices::dev.off()
    base::file.remove(bw_path)
    if (base::file.exists(pdf_path)) base::file.remove(pdf_path)
  }, add = TRUE)

  # Query regions covering all bins
  regions <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = base::seq(1L, 9001L, 1000L),
      end = base::seq(1000L, 10000L, 1000L)
    )
  )

  # Suppress plot output by directing to PDF null device
  grDevices::pdf(pdf_path)

  result <- plot_signal_histogram(
    bw_paths = base::c(Synthetic = bw_path),
    regions = regions
  )

  expect_regression_match(result, "plot_signal_histogram_default", tolerance = 1e-6)
})
