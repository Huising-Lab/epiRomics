library(testthat)
library(epiRomics)

test_that("epiRomics vignette workflow integration test", {
  require_extdata()
  # Build the database using helper with resolved paths
  db_file <- make_db_sheet()

  # Build the database
  database <- build_database(
    db_file = db_file,
    txdb_organism = "TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene",
    genome = "hg38",
    organism = "org.Hs.eg.db"
  )
  expect_s4_class(database, "epiRomicsS4")
  expect_true(length(database@annotations) > 0)
  expect_true(nrow(database@meta) > 0)

  # Identify putative enhancers
  enhancers <- find_enhancers_by_comarks(
    database = database,
    histone_mark_1 = "h3k4me1",
    histone_mark_2 = "h3k27ac"
  )
  expect_s4_class(enhancers, "epiRomicsS4")
  expect_true(length(enhancers@annotations) > 0)

  # Filter enhancers with FANTOM
  enhancers_fantom <- filter_enhancers(
    putative_enhancers = enhancers,
    database = database,
    type = "hg38_custom_fantom"
  )
  expect_s4_class(enhancers_fantom, "epiRomicsS4")
  expect_true(length(enhancers_fantom@annotations) > 0)

  # Filter enhancers with regulome_active
  enhancers_regulome_active <- filter_enhancers(
    putative_enhancers = enhancers,
    database = database,
    type = "hg38_custom_regulome_active"
  )
  expect_s4_class(enhancers_regulome_active, "epiRomicsS4")
  expect_true(length(enhancers_regulome_active@annotations) > 0)

  # Filter enhancers with regulome_super
  enhancers_regulome_super <- filter_enhancers(
    putative_enhancers = enhancers,
    database = database,
    type = "hg38_custom_regulome_super"
  )
  expect_s4_class(enhancers_regulome_super, "epiRomicsS4")
  expect_true(length(enhancers_regulome_super@annotations) > 0)
})
