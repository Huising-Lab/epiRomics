# Regression baseline tests for get_regions_of_interest()
# Verifies output stability across 3 input modes: granges, genelist, bed
# Uses expect_regression_match() from helper-synthetic-builders.R

library(testthat)
library(epiRomics)
library(GenomicRanges)
library(IRanges)

# --- Shared dB builder with SYMBOL metadata ---
.make_roi_dB <- function() {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = base::c(1000L, 6000L, 11000L, 16000L, 21000L,
                       26000L, 31000L, 36000L, 41000L, 46000L),
      width = 500L
    )
  )
  gr$SYMBOL <- base::paste0("GENE", base::seq_len(10L))
  gr$type <- base::paste0("hg38_custom_h3k4me1")

  meta <- base::data.frame(
    name = "h3k4me1",
    type = "histone",
    stringsAsFactors = FALSE
  )

  methods::new("epiRomicsS4",
    annotations = gr,
    meta = meta,
    txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene",
    organism = "org.Hs.eg.db",
    genome = "hg38"
  )
}

# ============================================================================
# Test 1: GRanges mode regression baseline
# ============================================================================
test_that("regions_of_interest regression: granges mode", {
  dB <- .make_roi_dB()

  # Query overlapping regions 1, 2, and 3 (start = 1000, 6000, 11000)
  query_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = base::c(800L, 5800L, 10800L),
      end = base::c(1200L, 6200L, 11200L)
    )
  )

  # Expected benign noise: get_regions_of_interest emits info-level progress
  # messages ("Overlap found for N regions"). The assertion is a snapshot of
  # the returned annotations data.frame.
  result <- base::suppressMessages(
    get_regions_of_interest(dB, query_gr)
  )

  # Convert to data.frame for stable serialization
  result_df <- base::as.data.frame(result@annotations)

  expect_regression_match(result_df, "roi_granges")
})

# ============================================================================
# Test 2: Genelist mode regression baseline
# ============================================================================
test_that("regions_of_interest regression: genelist mode", {
  dB <- .make_roi_dB()

  # Expected benign noise: get_regions_of_interest genelist-mode progress
  # messages ("Matched N genes"). Not under test — assertion is on snapshot.
  result <- base::suppressMessages(
    get_regions_of_interest(
      dB,
      input_type = "genelist",
      gene_list = base::c("GENE1", "GENE5", "GENE8")
    )
  )

  result_df <- base::as.data.frame(result@annotations)

  expect_regression_match(result_df, "roi_genelist")
})

# ============================================================================
# Test 3: BED mode regression baseline
# ============================================================================
test_that("regions_of_interest regression: bed mode", {
  dB <- .make_roi_dB()

  # Create a temp BED file covering regions 2 and 4
  bed_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = base::c(5800L, 15800L),
      end = base::c(6200L, 16200L)
    )
  )
  bed_file <- base::tempfile(fileext = ".bed")
  base::on.exit(base::file.remove(bed_file), add = TRUE)
  rtracklayer::export(bed_gr, bed_file, format = "BED")

  # Expected benign noise: get_regions_of_interest BED-mode progress
  # messages ("Loading BED", "Overlap found"). Not under test — assertion is
  # on the result snapshot.
  result <- base::suppressMessages(
    get_regions_of_interest(
      dB,
      input_type = "bed",
      bed_path = bed_file
    )
  )

  result_df <- base::as.data.frame(result@annotations)

  expect_regression_match(result_df, "roi_bed")
})
