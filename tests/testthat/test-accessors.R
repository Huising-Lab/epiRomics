# =============================================================================
# Accessor tests for epiRomicsS4 public getters and setters
#
# Covers every getter/setter pair (annotations, meta, txdb, organism, genome),
# validity enforcement on empty-string assignment, round-trip behavior, and
# generic-extension behavior (organism/genome extend BiocGenerics and
# GenomeInfoDb respectively).
#
# Fixtures: prefer make_example_database() over hand-built objects.
# =============================================================================

library(testthat)
library(epiRomics)
library(GenomicRanges)
library(IRanges)

# ============================================================================
# Getters: annotations()
# ============================================================================
test_that("annotations() returns the populated GRanges from example hg38 db", {
  db <- make_example_database("hg38")

  result <- annotations(db)

  testthat::expect_s4_class(result, "GRanges")
  testthat::expect_gt(base::length(result), 0L)
  testthat::expect_identical(result, db@annotations)
})

test_that("annotations() returns empty GRanges on a fresh epiRomicsS4", {
  obj <- methods::new("epiRomicsS4")

  result <- annotations(obj)

  testthat::expect_s4_class(result, "GRanges")
  testthat::expect_equal(base::length(result), 0L)
})

# ============================================================================
# Getters: meta()
# ============================================================================
test_that("meta() returns data.frame with name and type columns", {
  db <- make_example_database("hg38")

  result <- meta(db)

  testthat::expect_s3_class(result, "data.frame")
  testthat::expect_true(base::all(base::c("name", "type") %in% base::colnames(result)))
  testthat::expect_gt(base::nrow(result), 0L)
  testthat::expect_identical(result, db@meta)
})

test_that("meta() returns empty data.frame on a fresh epiRomicsS4", {
  obj <- methods::new("epiRomicsS4")

  result <- meta(obj)

  testthat::expect_s3_class(result, "data.frame")
  testthat::expect_equal(base::nrow(result), 0L)
})

# ============================================================================
# Getters: txdb()
# ============================================================================
test_that("txdb() returns the TxDb identifier string from example hg38 db", {
  db <- make_example_database("hg38")

  result <- txdb(db)

  testthat::expect_type(result, "character")
  testthat::expect_equal(
    result,
    "TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene"
  )
})

test_that("txdb() returns empty character on a fresh epiRomicsS4", {
  obj <- methods::new("epiRomicsS4")

  result <- txdb(obj)

  testthat::expect_type(result, "character")
  testthat::expect_equal(base::length(result), 0L)
})

# ============================================================================
# Getters: organism()
# ============================================================================
test_that("organism() returns org.Hs.eg.db from example hg38 db", {
  db <- make_example_database("hg38")

  result <- organism(db)

  testthat::expect_type(result, "character")
  testthat::expect_equal(result, "org.Hs.eg.db")
})

test_that("organism() returns empty character on a fresh epiRomicsS4", {
  obj <- methods::new("epiRomicsS4")

  result <- organism(obj)

  testthat::expect_type(result, "character")
  testthat::expect_equal(base::length(result), 0L)
})

# ============================================================================
# Getters: genome()
# ============================================================================
test_that("genome() returns 'hg38' from hg38 example db", {
  db <- make_example_database("hg38")

  result <- genome(db)

  testthat::expect_type(result, "character")
  testthat::expect_equal(result, "hg38")
})

test_that("genome() returns 'mm10' from mm10 example db", {
  db <- make_example_database("mm10")

  result <- genome(db)

  testthat::expect_type(result, "character")
  testthat::expect_equal(result, "mm10")
})

test_that("genome() returns empty character on a fresh epiRomicsS4", {
  obj <- methods::new("epiRomicsS4")

  result <- genome(obj)

  testthat::expect_type(result, "character")
  testthat::expect_equal(base::length(result), 0L)
})

# ============================================================================
# Setters: annotations<-
# ============================================================================
test_that("annotations<- updates slot and returns epiRomicsS4", {
  db <- make_example_database("hg38")
  original_length <- base::length(annotations(db))

  new_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = base::c(100L, 200L), end = base::c(150L, 250L))
  )
  new_gr$type <- base::c("hg38_custom_test", "hg38_custom_test")

  db <- `annotations<-`(db, value = new_gr)

  testthat::expect_s4_class(db, "epiRomicsS4")
  testthat::expect_equal(base::length(annotations(db)), 2L)
  testthat::expect_false(base::length(annotations(db)) == original_length)
  testthat::expect_identical(annotations(db), new_gr)
})

# ============================================================================
# Setters: meta<-
# ============================================================================
test_that("meta<- updates slot with extra column and returns epiRomicsS4", {
  db <- make_example_database("hg38")
  new_meta <- meta(db)
  new_meta$extra <- "annotated"

  db <- `meta<-`(db, value = new_meta)

  testthat::expect_s4_class(db, "epiRomicsS4")
  testthat::expect_true("extra" %in% base::colnames(meta(db)))
  testthat::expect_equal(base::ncol(meta(db)), base::ncol(new_meta))
  testthat::expect_identical(meta(db), new_meta)
})

# ============================================================================
# Setters: txdb<-
# ============================================================================
test_that("txdb<- updates slot with a new TxDb identifier and round-trips", {
  db <- make_example_database("hg38")
  new_txdb <- "TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene"

  db <- `txdb<-`(db, value = new_txdb)

  testthat::expect_s4_class(db, "epiRomicsS4")
  testthat::expect_equal(txdb(db), new_txdb)
})

# ============================================================================
# Setters: organism<-
# ============================================================================
test_that("organism<- updates slot with a new org.db name and round-trips", {
  db <- make_example_database("hg38")

  organism(db) <- "org.Mm.eg.db"

  testthat::expect_s4_class(db, "epiRomicsS4")
  testthat::expect_equal(organism(db), "org.Mm.eg.db")
})

# ============================================================================
# Setters: genome<-
# ============================================================================
test_that("genome<- switches hg38 -> mm10 and round-trips", {
  db <- make_example_database("hg38")
  testthat::expect_equal(genome(db), "hg38")

  genome(db) <- "mm10"

  testthat::expect_s4_class(db, "epiRomicsS4")
  testthat::expect_equal(genome(db), "mm10")
})

# ============================================================================
# Validity: empty-string rejection on genome and txdb
# ============================================================================
test_that("genome<- rejects empty string via validity", {
  db <- make_example_database("hg38")

  testthat::expect_error(
    {
      genome(db) <- ""
    },
    "genome must be a non-empty character string"
  )
})

test_that("txdb<- rejects empty string via validity", {
  db <- make_example_database("hg38")

  testthat::expect_error(
    {
      txdb(db) <- ""
    },
    "txdb must be a non-empty character string"
  )
})

test_that("genome<- accepts valid assembly name without error", {
  db <- make_example_database("hg38")

  testthat::expect_silent({
    genome(db) <- "mm10"
  })
  testthat::expect_equal(genome(db), "mm10")
})

# ============================================================================
# Generic extension behavior: organism() and genome() extend upstream generics
# ============================================================================
test_that("organism and genome are S4 generics with methods on epiRomicsS4", {
  testthat::expect_true(methods::isGeneric("organism"))
  testthat::expect_true(methods::isGeneric("genome"))
  testthat::expect_true(methods::existsMethod("organism", "epiRomicsS4"))
  testthat::expect_true(methods::existsMethod("genome", "epiRomicsS4"))
  testthat::expect_true(methods::existsMethod("organism<-", "epiRomicsS4"))
  testthat::expect_true(methods::existsMethod("genome<-", "epiRomicsS4"))
})
