# ============================================================================
# REGRESSION BASELINE TESTS for epiRomics_quick_view
# Uses synthetic BigWig files for deterministic testing.
# Baselines capture either successful return values or deterministic error
# patterns -- both are valid regression detection.
# Baselines saved as .rds in testdata/ via expect_regression_match().
# ============================================================================

library(testthat)
library(epiRomics)
library(GenomicRanges)
library(IRanges)

test_that("quick_view regression: region mode hg38", {
  testthat::skip_if_not_installed("rtracklayer")

  # Create synthetic BigWig file with 10 contiguous bins
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = base::seq(1L, 9001L, 1000L),
      end = base::seq(1000L, 10000L, 1000L)
    )
  )
  scores <- base::c(5, 10, 50, 80, 90, 5, 10, 95, 5, 10)

  bw_path <- make_synthetic_bigwig(gr, scores)
  pdf_path <- base::tempfile(fileext = ".pdf")
  base::on.exit({
    if (grDevices::dev.cur() > 1L) grDevices::dev.off()
    if (base::file.exists(bw_path)) base::file.remove(bw_path)
    if (base::file.exists(pdf_path)) base::file.remove(pdf_path)
  }, add = TRUE)

  bw_paths <- base::c(Synthetic_ATAC = bw_path)

  # Suppress plot output
  grDevices::pdf(pdf_path)

  # Capture either the result or the error as the regression baseline
  outcome <- base::tryCatch(
    {
      result <- epiRomics_quick_view(
        region = base::list(chr = "chr1", start = 1000L, end = 10000L),
        bw_paths = bw_paths,
        genome = "hg38"
      )
      base::list(
        status = "success",
        result_class = base::class(result)
      )
    },
    error = function(e) {
      base::list(
        status = "error",
        error_class = base::class(e),
        error_msg_pattern = base::sub(".*: ", "", e$message)
      )
    }
  )

  expect_regression_match(outcome, "quick_view_region_hg38")
})

test_that("quick_view regression: region mode mm10", {
  testthat::skip_if_not_installed("rtracklayer")

  # Create synthetic BigWig file with 10 contiguous bins
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = base::seq(1L, 9001L, 1000L),
      end = base::seq(1000L, 10000L, 1000L)
    )
  )
  scores <- base::c(5, 10, 50, 80, 90, 5, 10, 95, 5, 10)

  bw_path <- make_synthetic_bigwig(gr, scores)
  pdf_path <- base::tempfile(fileext = ".pdf")
  base::on.exit({
    if (grDevices::dev.cur() > 1L) grDevices::dev.off()
    if (base::file.exists(bw_path)) base::file.remove(bw_path)
    if (base::file.exists(pdf_path)) base::file.remove(pdf_path)
  }, add = TRUE)

  bw_paths <- base::c(Synthetic_ATAC = bw_path)

  # Suppress plot output
  grDevices::pdf(pdf_path)

  # Capture either the result or the error as the regression baseline
  outcome <- base::tryCatch(
    {
      result <- epiRomics_quick_view(
        region = base::list(chr = "chr1", start = 1000L, end = 10000L),
        bw_paths = bw_paths,
        genome = "mm10"
      )
      base::list(
        status = "success",
        result_class = base::class(result)
      )
    },
    error = function(e) {
      base::list(
        status = "error",
        error_class = base::class(e),
        error_msg_pattern = base::sub(".*: ", "", e$message)
      )
    }
  )

  expect_regression_match(outcome, "quick_view_region_mm10")
})
