library(testthat)
library(epiRomics)

# The validator must accept genome strings that match the TxDb's
# recorded assembly and reject mismatches. These tests exercise
# the internal helper via a direct :::-qualified call because
# validate_genome_matches_txdb is @noRd (not exported).

testthat::test_that("validate_genome_matches_txdb accepts matching hg38 TxDb", {
  # TxDb.Hsapiens.UCSC.hg38.knownGene is a Suggests; keep that skip.
  # GenomeInfoDb is an Imports — no skip.
  testthat::skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  result <- epiRomics:::validate_genome_matches_txdb(
    genome = "hg38",
    txdb_string = paste0(
      "TxDb.Hsapiens.UCSC.hg38.knownGene::",
      "TxDb.Hsapiens.UCSC.hg38.knownGene"
    ),
    function_name = "test"
  )
  testthat::expect_true(result)
})

testthat::test_that("validate_genome_matches_txdb rejects mismatched genome", {
  # TxDb.Hsapiens.UCSC.hg38.knownGene is a Suggests; keep that skip.
  # GenomeInfoDb is an Imports — no skip.
  testthat::skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  testthat::expect_error(
    epiRomics:::validate_genome_matches_txdb(
      genome = "mm10",
      txdb_string = paste0(
        "TxDb.Hsapiens.UCSC.hg38.knownGene::",
        "TxDb.Hsapiens.UCSC.hg38.knownGene"
      ),
      function_name = "test"
    ),
    regexp = "does not match TxDb genome"
  )
})

testthat::test_that("validate_genome_matches_txdb validates inputs", {
  testthat::expect_error(
    epiRomics:::validate_genome_matches_txdb(
      genome = NULL,
      txdb_string = "a::b",
      function_name = "test"
    ),
    regexp = "genome must be a character string"
  )
  testthat::expect_error(
    epiRomics:::validate_genome_matches_txdb(
      genome = "hg38",
      txdb_string = "",
      function_name = "test"
    ),
    regexp = "txdb"
  )
})

testthat::test_that("validate_genome_matches_txdb rejects malformed TxDb string", {
  testthat::expect_error(
    epiRomics:::validate_genome_matches_txdb(
      genome = "hg38",
      txdb_string = "not_a_valid_spec",
      function_name = "test"
    ),
    regexp = "Invalid TxDb format"
  )
})
