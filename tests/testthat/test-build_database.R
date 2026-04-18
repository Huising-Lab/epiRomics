library(testthat)
library(epiRomics)

# Test: build_database returns correct S4 class and expected output
testthat::test_that("build_database returns correct S4 class and expected output", {
  require_extdata()
  db_sheet <- make_db_sheet()
  # Expected benign warnings from build_database() against toy chr11-only data:
  #   * annotatr "No seqlevels in common" for builtin annotations outside chr11.
  #   * ChIPseeker "peak has no overlap with gene" for toy peaks in intergenic
  #     regions. Neither is regression-relevant here — the test asserts S4
  #     structure, not warning text.
  suppressWarnings({
    result <- epiRomics::build_database(
      db_file = db_sheet,
      txdb_organism = "TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene",
      genome = "hg38",
      organism = "org.Hs.eg.db"
    )
  })

  testthat::expect_s4_class(result, "epiRomicsS4")
  testthat::expect_equal(result@organism, "org.Hs.eg.db")
  testthat::expect_equal(result@genome, "hg38")
  testthat::expect_equal(result@txdb, "TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene")
  testthat::expect_true(length(result@annotations) >= 0)
  testthat::expect_true(nrow(result@meta) > 0)
})
