# Regression baseline tests for epiRomicsS4 S4 class constructor
# Verifies constructor output stability for hg38 and mm10 genome configurations
# Uses expect_regression_match() from helper-synthetic-builders.R

library(testthat)
library(epiRomics)
library(GenomicRanges)
library(IRanges)

# ============================================================================
# Test 1: epiRomicsS4 regression -- hg38 construction
# ============================================================================
test_that("epiRomicsS4 regression: hg38 construction", {
  # Build annotations with known structure
  gr <- GenomicRanges::GRanges(
    seqnames = base::c("chr1", "chr1", "chr2"),
    ranges = IRanges::IRanges(
      start = base::c(1000L, 5000L, 10000L),
      end = base::c(2000L, 6000L, 11000L)
    )
  )
  gr$type <- base::c("hg38_custom_h3k4me1", "hg38_custom_h3k27ac",
                      "hg38_custom_h3k4me1")
  gr$SYMBOL <- base::c("GENE1", "GENE2", "GENE3")

  meta <- base::data.frame(
    name = base::c("h3k4me1", "h3k27ac"),
    type = base::c("histone", "histone"),
    stringsAsFactors = FALSE
  )

  obj <- methods::new("epiRomicsS4",
    annotations = gr,
    meta = meta,
    genome = "hg38",
    txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene",
    organism = "org.Hs.eg.db"
  )

  # Convert to serializable list for stable comparison
  result <- base::list(
    annotations = base::as.data.frame(obj@annotations),
    meta = obj@meta,
    genome = obj@genome,
    txdb = obj@txdb,
    organism = obj@organism,
    n_annotations = base::length(obj@annotations)
  )

  expect_regression_match(result, "epiRomicsS4_hg38")
})

# ============================================================================
# Test 2: epiRomicsS4 regression -- mm10 construction
# ============================================================================
test_that("epiRomicsS4 regression: mm10 construction", {
  gr <- GenomicRanges::GRanges(
    seqnames = base::c("chr1", "chr5"),
    ranges = IRanges::IRanges(
      start = base::c(2000L, 8000L),
      end = base::c(3000L, 9000L)
    )
  )
  gr$type <- base::c("mm10_custom_h3k4me1", "mm10_custom_h3k27ac")
  gr$SYMBOL <- base::c("Ins1", "Gcg")

  meta <- base::data.frame(
    name = base::c("h3k4me1", "h3k27ac"),
    type = base::c("histone", "histone"),
    stringsAsFactors = FALSE
  )

  obj <- methods::new("epiRomicsS4",
    annotations = gr,
    meta = meta,
    genome = "mm10",
    txdb = "TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene",
    organism = "org.Mm.eg.db"
  )

  result <- base::list(
    annotations = base::as.data.frame(obj@annotations),
    meta = obj@meta,
    genome = obj@genome,
    txdb = obj@txdb,
    organism = obj@organism,
    n_annotations = base::length(obj@annotations)
  )

  expect_regression_match(result, "epiRomicsS4_mm10")
})
