testthat::test_that("errors", {
  testthat::expect_s4_class(
    epiRomics_build_dB(
      epiRomics_db_file = system.file(
        "extdata",
        "example_epiRomics_Db_sheet_user_paths.csv",
        package = "epiRomics"
      ),
      txdb_organism = "TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene",
      epiRomics_genome = "hg38",
      epiRomics_organism = "org.Hs.eg.db"
    ),
    "epiRomicsS4"
  )
})