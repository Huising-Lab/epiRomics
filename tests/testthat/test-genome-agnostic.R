library(testthat)
library(epiRomics)

# Prove organism-agnosticism.
#
# epiRomics must accept ANY genome string that is consistent with the
# TxDb the caller supplies. hg38/mm10 remain only as convenience
# defaults in plot_quick_view; the core validator is generic.
#
# Every test that needs an optional TxDb uses skip_if_not_installed()
# so these never fail the Bioconductor Single Package Builder when a
# non-default TxDb is not present on the build machine.

testthat::test_that(
  "validator accepts rn6 (rat) TxDb — proves organism-agnosticism", {
    # TxDb.Rnorvegicus is an optional organism package (not a Suggests of
    # this package), so keep the skip. GenomeInfoDb is an Imports — no skip.
    testthat::skip_if_not_installed("TxDb.Rnorvegicus.UCSC.rn6.refGene")
    result <- epiRomics:::validate_genome_matches_txdb(
      genome = "rn6",
      txdb_string = base::paste0(
        "TxDb.Rnorvegicus.UCSC.rn6.refGene::",
        "TxDb.Rnorvegicus.UCSC.rn6.refGene"
      ),
      function_name = "test_agnostic"
    )
    testthat::expect_true(result)
  }
)

testthat::test_that(
  "validator accepts dm6 (fruit fly) TxDb — proves non-mammalian support", {
    # TxDb.Dmelanogaster is an optional organism package; keep the skip.
    # GenomeInfoDb is an Imports — no skip.
    testthat::skip_if_not_installed("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
    result <- epiRomics:::validate_genome_matches_txdb(
      genome = "dm6",
      txdb_string = base::paste0(
        "TxDb.Dmelanogaster.UCSC.dm6.ensGene::",
        "TxDb.Dmelanogaster.UCSC.dm6.ensGene"
      ),
      function_name = "test_agnostic"
    )
    testthat::expect_true(result)
  }
)

testthat::test_that(
  "validator rejects hg38 genome with rn6 TxDb — cross-organism mismatch", {
    # TxDb.Rnorvegicus is an optional organism package; keep the skip.
    # GenomeInfoDb is an Imports — no skip.
    testthat::skip_if_not_installed("TxDb.Rnorvegicus.UCSC.rn6.refGene")
    testthat::expect_error(
      epiRomics:::validate_genome_matches_txdb(
        genome = "hg38",
        txdb_string = base::paste0(
          "TxDb.Rnorvegicus.UCSC.rn6.refGene::",
          "TxDb.Rnorvegicus.UCSC.rn6.refGene"
        ),
        function_name = "test_agnostic"
      ),
      regexp = "does not match TxDb genome"
    )
  }
)

testthat::test_that(
  "validator source has no hardcoded organism whitelist", {
    # This introspective test needs no external TxDb and runs everywhere.
    # It guards against regression: a future refactor that hardcodes a
    # genome list would fail here.
    validator_src <- base::paste(
      base::deparse(epiRomics:::validate_genome_matches_txdb),
      collapse = "\n"
    )
    # Must use GenomeInfoDb::genome() for dynamic assembly lookup
    testthat::expect_true(
      base::grepl("GenomeInfoDb::genome", validator_src, fixed = TRUE)
    )
    # Must NOT use match.arg() against a fixed genome vector
    testthat::expect_false(
      base::grepl("match.arg", validator_src, fixed = TRUE)
    )
  }
)

testthat::test_that(
  "plot_quick_view no longer whitelists genome via match.arg", {
    # plot_quick_view previously used match.arg(genome, c("hg38","mm10")).
    # It now uses a switch() convenience default that errors informatively
    # when the user supplies a non-default genome without a txdb. Confirm
    # the whitelist is gone.
    qv_src <- base::paste(
      base::deparse(plot_quick_view),
      collapse = "\n"
    )
    testthat::expect_false(
      base::grepl(
        "match.arg(genome", qv_src, fixed = TRUE
      )
    )
  }
)
