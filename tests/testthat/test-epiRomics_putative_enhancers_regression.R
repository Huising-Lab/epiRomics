## Regression baseline tests for epiRomics_putative_enhancers
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

  result <- base::suppressWarnings(base::suppressMessages(
    epiRomics_putative_enhancers(dB)
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

  cs <- base::suppressWarnings(base::suppressMessages(
    epiRomics_chromatin_states(dB, refine_by_tss = FALSE)
  ))

  result <- base::suppressWarnings(base::suppressMessages(
    epiRomics_putative_enhancers(dB, chromatin_states = cs)
  ))

  testthat::expect_true(base::is.data.frame(result))
  testthat::expect_true(base::nrow(result) > 0L)

  # Regression baseline comparison
  expect_regression_match(result, "putative_enhancers_precomputed_cs")
})
