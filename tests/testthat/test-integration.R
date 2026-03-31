library(testthat)
library(epiRomics)

test_that("epiRomics vignette workflow integration test", {
  skip_if_no_extdata()
  # Build the database using helper with resolved paths
  db_file <- make_db_sheet()

  # Build the database
  epiRomics_dB <- epiRomics_build_dB(
    epiRomics_db_file = db_file,
    txdb_organism = "TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene",
    epiRomics_genome = "hg38",
    epiRomics_organism = "org.Hs.eg.db"
  )
  expect_s4_class(epiRomics_dB, "epiRomicsS4")
  expect_true(length(epiRomics_dB@annotations) > 0)
  expect_true(nrow(epiRomics_dB@meta) > 0)

  # Identify putative enhancers
  enhancers <- epiRomics_enhancers_co_marks(
    epiRomics_dB = epiRomics_dB,
    epiRomics_histone_mark_1 = "h3k4me1",
    epiRomics_histone_mark_2 = "h3k27ac"
  )
  expect_s4_class(enhancers, "epiRomicsS4")
  expect_true(length(enhancers@annotations) > 0)

  # Filter enhancers with FANTOM
  enhancers_fantom <- epiRomics_enhancers_filter(
    epiRomics_putative_enhancers = enhancers,
    epiRomics_dB = epiRomics_dB,
    epiRomics_type = "hg38_custom_fantom"
  )
  expect_s4_class(enhancers_fantom, "epiRomicsS4")
  expect_true(length(enhancers_fantom@annotations) > 0)

  # Filter enhancers with regulome_active
  enhancers_regulome_active <- epiRomics_enhancers_filter(
    epiRomics_putative_enhancers = enhancers,
    epiRomics_dB = epiRomics_dB,
    epiRomics_type = "hg38_custom_regulome_active"
  )
  expect_s4_class(enhancers_regulome_active, "epiRomicsS4")
  expect_true(length(enhancers_regulome_active@annotations) > 0)

  # Filter enhancers with regulome_super
  enhancers_regulome_super <- epiRomics_enhancers_filter(
    epiRomics_putative_enhancers = enhancers,
    epiRomics_dB = epiRomics_dB,
    epiRomics_type = "hg38_custom_regulome_super"
  )
  expect_s4_class(enhancers_regulome_super, "epiRomicsS4")
  expect_true(length(enhancers_regulome_super@annotations) > 0)
})
