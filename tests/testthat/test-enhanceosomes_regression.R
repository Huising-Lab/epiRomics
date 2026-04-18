## Regression baseline tests for find_enhanceosomes
##
## Pattern: Run function on synthetic data, save output as .rds baseline,
## compare on future runs. First run creates the baseline; subsequent runs
## verify output has not changed.

library(testthat)
library(epiRomics)
library(GenomicRanges)
library(IRanges)

# ---------------------------------------------------------------------------
# Test 1: enhanceosome regression: default
# ---------------------------------------------------------------------------
test_that("enhanceosome regression: default", {
  dB <- make_synthetic_dB_full()

  # Create a "putative enhancers" epiRomicsS4 object.
  # Use regions that overlap the chip TF annotations (chr1:1000-2000,
  # chr1:5000-6000, chr1:50000-51000) so countOverlaps finds actual overlaps.
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

  # find_enhanceosomes adds TF overlap counts as mcols and calls
  # ChIPseeker::annotatePeak which needs a TxDb. This should work with
  # synthetic data on hg38 coordinates.
  # Expected benign noise:
  #   * ChIPseeker progress messages + "peak has no overlap with gene" for
  #     synthetic peaks outside real hg38 gene bodies.
  # Assertion is on TF1/TF2/ChIP_Hits columns, not warning text.
  result <- base::suppressWarnings(base::suppressMessages(
    find_enhanceosomes(pe_dB, dB)
  ))

  # Sanity: result should be an epiRomicsS4 object
  testthat::expect_true(methods::is(result, "epiRomicsS4"))
  testthat::expect_true(base::length(result@annotations) > 0L)

  # Extract annotations as data.frame for stable .rds comparison
  result_df <- base::as.data.frame(result@annotations)

  # Verify TF columns were added
  testthat::expect_true("TF1" %in% base::colnames(result_df))
  testthat::expect_true("TF2" %in% base::colnames(result_df))
  testthat::expect_true("ChIP_Hits" %in% base::colnames(result_df))

  # Regression baseline comparison
  expect_regression_match(result_df, "enhanceosome_default")
})
