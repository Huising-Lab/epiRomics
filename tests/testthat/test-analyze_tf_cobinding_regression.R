## Regression baseline tests for analyze_tf_cobinding
##
## Pattern: Run function on synthetic data, save output as .rds baseline,
## compare on future runs. First run creates the baseline; subsequent runs
## verify output has not changed.

library(testthat)
library(epiRomics)
library(GenomicRanges)
library(IRanges)

# ---------------------------------------------------------------------------
# Test 1: tf_cobinding regression: fisher method
# ---------------------------------------------------------------------------
test_that("tf_cobinding regression: fisher method", {
  dB <- make_synthetic_dB_full()

  # Build an enhanceosome-like epiRomicsS4 with TF overlap columns.
  # We run find_enhanceosomes first to get proper input.
  pe_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = base::c(1000L, 5000L, 20000L, 50000L),
      end = base::c(2000L, 6000L, 21000L, 51000L)
    )
  )
  pe_gr$type <- base::paste0("hg38_custom_putative")

  pe_dB <- methods::new("epiRomicsS4",
    annotations = pe_gr,
    meta = base::data.frame(
      name = "putative",
      type = "histone",
      stringsAsFactors = FALSE
    ),
    genome = "hg38",
    txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene",
    organism = "org.Hs.eg.db"
  )

  # Expected benign noise from find_enhanceosomes against synthetic peaks:
  # ChIPseeker "peak has no overlap with gene" for synthetic chr1 peaks
  # outside real hg38 gene bodies, and an annotatPeak progress message.
  # This regression test asserts the downstream result shape, not the
  # construction warnings.
  enhanceosome <- base::suppressWarnings(base::suppressMessages(
    find_enhanceosomes(pe_dB, dB)
  ))

  # analyze_tf_cobinding requires at least 2 ChIP TFs in the dB
  # (which make_synthetic_dB_full provides: TF1, TF2)
  # Expected benign noise: progress messages from analyze_tf_cobinding and
  # Fisher-test "approximation may be incorrect" warnings on small synthetic
  # 4-region matrices. Assertion is on $method and pairwise shape.
  result <- base::suppressWarnings(base::suppressMessages(
    analyze_tf_cobinding(enhanceosome, dB, method = "fisher")
  ))

  # Sanity: result should be a list with expected components
  testthat::expect_type(result, "list")
  testthat::expect_true("pairwise" %in% base::names(result))
  testthat::expect_true("presence_matrix" %in% base::names(result))
  testthat::expect_true("clustering" %in% base::names(result))
  testthat::expect_true("tf_names" %in% base::names(result))
  testthat::expect_true("n_regions" %in% base::names(result))
  testthat::expect_true("method" %in% base::names(result))
  testthat::expect_equal(result$method, "fisher")

  # For stable regression comparison, save pairwise results and
  # presence_matrix (hclust objects may vary). Extract stable components.
  stable_result <- base::list(
    pairwise = result$pairwise,
    presence_matrix = result$presence_matrix,
    tf_names = result$tf_names,
    n_regions = result$n_regions,
    method = result$method
  )

  # Regression baseline comparison
  expect_regression_match(stable_result, "tf_cobinding_fisher")
})
