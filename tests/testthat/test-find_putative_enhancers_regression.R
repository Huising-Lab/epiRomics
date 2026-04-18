## Regression baseline tests for find_putative_enhancers
##
## Pattern: Run function on synthetic data, save output as .rds baseline,
## compare on future runs. First run creates the baseline; subsequent runs
## verify output has not changed.

library(testthat)
library(epiRomics)
library(GenomicRanges)
library(IRanges)

# ---------------------------------------------------------------------------
# Test 1: putative_enhancers default (auto-compute chromatin states)
# ---------------------------------------------------------------------------
test_that("putative_enhancers regression: default", {
  dB <- make_synthetic_dB_full()

  # Expected benign noise from find_putative_enhancers against synthetic data:
  #   * classify_chromatin_states progress messages.
  #   * ChIPseeker "peak has no overlap with gene" for synthetic peaks that
  #     fall outside real hg38 gene bodies.
  # Assertion is on expected_cols presence, not warning text.
  result <- base::suppressWarnings(base::suppressMessages(
    find_putative_enhancers(dB)
  ))

  # Sanity: must return a data.frame with expected columns
  testthat::expect_true(base::is.data.frame(result))
  expected_cols <- base::c(
    "putative_id", "chr", "start", "end", "width", "source",
    "chromatin_state", "chromatin_state_detail", "histone_marks",
    "n_histone_marks", "h2az", "tf_names", "n_tfs"
  )
  testthat::expect_true(
    base::all(expected_cols %in% base::colnames(result))
  )
  testthat::expect_true(base::nrow(result) > 0L)

  # Regression baseline comparison
  expect_regression_match(result, "putative_enhancers_default")
})

# ---------------------------------------------------------------------------
# Test 2: putative_enhancers with pre-computed chromatin_states
# ---------------------------------------------------------------------------
test_that("putative_enhancers regression: with pre-computed chromatin_states", {
  dB <- make_synthetic_dB_full()

  # Expected benign noise: classify_chromatin_states progress messages when
  # building a chromatin-states fixture for the subsequent test step. Not
  # under test here — the assertion is on find_putative_enhancers' output.
  cs <- base::suppressWarnings(base::suppressMessages(
    classify_chromatin_states(dB, refine_by_tss = FALSE)
  ))

  # Expected benign noise: find_putative_enhancers progress messages plus
  # ChIPseeker "peak has no overlap with gene" on synthetic chr1 peaks.
  result <- base::suppressWarnings(base::suppressMessages(
    find_putative_enhancers(dB, chromatin_states = cs)
  ))

  testthat::expect_true(base::is.data.frame(result))
  testthat::expect_true(base::nrow(result) > 0L)

  # Regression baseline comparison
  expect_regression_match(result, "putative_enhancers_precomputed_cs")
})
