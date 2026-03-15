library(testthat)
library(epiRomics)

# Test: epiRomics_build_dB_2 returns correct S4 class and expected output
testthat::test_that("epiRomics_build_dB_2 returns correct S4 class and expected output", {
  skip_if_no_extdata()
  db_sheet <- make_db_sheet()
  suppressWarnings({
    result <- epiRomics::epiRomics_build_dB_2(
      epiRomics_db_file = db_sheet,
      txdb_organism = "TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene",
      epiRomics_genome = "hg38",
      epiRomics_organism = "org.Hs.eg.db"
    )
  })

  testthat::expect_s4_class(result, "epiRomicsS4")
  testthat::expect_equal(result@organism, "org.Hs.eg.db")
  testthat::expect_equal(result@genome, "hg38")
  testthat::expect_equal(result@txdb, "TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene")
  testthat::expect_true(length(result@annotations) >= 0)
  testthat::expect_true(nrow(result@meta) > 0)
})
