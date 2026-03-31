## Regression baseline tests for epiRomics_tf_overlap
##
## Pattern: Run function on synthetic data, save output as .rds baseline,
## compare on future runs. First run creates the baseline; subsequent runs
## verify output has not changed.

library(testthat)
library(epiRomics)
library(GenomicRanges)
library(IRanges)

# ---------------------------------------------------------------------------
# Test 1: tf_overlap regression: default (self-overlap)
# ---------------------------------------------------------------------------
test_that("tf_overlap regression: default", {
  dB <- make_synthetic_dB_full()

  # epiRomics_tf_overlap takes two epiRomicsS4 dBs and computes pairwise

  # overlap between their ChIP TF annotations. Using the same dB for both
  # gives deterministic self-overlap results.
  result <- base::suppressWarnings(base::suppressMessages(
    epiRomics_tf_overlap(dB, dB)
  ))

  # Sanity: result should be a list with expected components
  testthat::expect_type(result, "list")
  testthat::expect_true("overlap_matrix" %in% base::names(result))
  testthat::expect_true("jaccard_matrix" %in% base::names(result))
  testthat::expect_true("summary" %in% base::names(result))
  testthat::expect_true("tf_names" %in% base::names(result))
  testthat::expect_equal(base::length(result$tf_names), 2L)

  # Overlap matrix should be 2x2
  testthat::expect_equal(base::dim(result$overlap_matrix), base::c(2L, 2L))

  # Self-overlap diagonal should be 1.0
  testthat::expect_equal(
    base::as.numeric(base::diag(result$overlap_matrix)),
    base::c(1.0, 1.0)
  )

  # Regression baseline comparison
  expect_regression_match(result, "tf_overlap_default")
})

# ---------------------------------------------------------------------------
# Test 2: tf_overlap regression: with custom regions
# ---------------------------------------------------------------------------
test_that("tf_overlap regression: with custom regions", {
  dB <- make_synthetic_dB_full()

  # Create a custom regions GRanges restricting analysis to a subset
  custom_regions <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = base::c(1000L, 5000L, 50000L),
      end = base::c(2000L, 6000L, 51000L)
    )
  )

  result <- base::suppressWarnings(base::suppressMessages(
    epiRomics_tf_overlap(dB, dB, regions = custom_regions)
  ))

  # Sanity checks
  testthat::expect_type(result, "list")
  testthat::expect_true("overlap_matrix" %in% base::names(result))
  testthat::expect_equal(base::length(result$tf_names), 2L)

  # Regression baseline comparison
  expect_regression_match(result, "tf_overlap_regions")
})
