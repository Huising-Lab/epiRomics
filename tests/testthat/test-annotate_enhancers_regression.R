## Regression baseline tests for annotate_enhancers
##
## Pattern: Run function on synthetic data, save output as .rds baseline,
## compare on future runs. First run creates the baseline; subsequent runs
## verify output has not changed.

library(testthat)
library(epiRomics)
library(GenomicRanges)
library(IRanges)

# ---------------------------------------------------------------------------
# Test 1: annotate_putative default (functional database overlap)
# ---------------------------------------------------------------------------
test_that("annotate_putative regression: default", {
  dB <- make_synthetic_dB_full()

  # Create putative enhancers data.frame using helper
  pe_df <- make_synthetic_putative_enhancers(n = 5L)

  # annotate_enhancers checks for functional databases in dB.
  # Our synthetic dB has a "fantom" functional entry, so it will compute

  # overlap columns for that database.
  # Expected benign noise: annotate_enhancers emits progress messages and
  # countOverlaps can warn about seqlevel style mismatches when the synthetic
  # fantom entry uses UCSC-style "chr1". Assertion is on returned columns,
  # not on message/warning text.
  result <- base::suppressWarnings(base::suppressMessages(
    annotate_enhancers(pe_df, dB)
  ))

  # Sanity: must return a data.frame with additional annotation columns
  testthat::expect_true(base::is.data.frame(result))
  testthat::expect_true(base::nrow(result) > 0L)
  # Should have added overlap columns + n_databases + novel
  testthat::expect_true("n_databases" %in% base::colnames(result))
  testthat::expect_true("novel" %in% base::colnames(result))
  # Should have at least 1 overlap column for "fantom"
  overlap_cols <- base::grep("^overlap_", base::colnames(result), value = TRUE)
  testthat::expect_true(base::length(overlap_cols) >= 1L)

  # Regression baseline comparison
  expect_regression_match(result, "annotate_putative_default")
})
