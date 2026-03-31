# Regression baseline tests for epiRomics_build_dB()
# Verifies output stability with fully synthetic CSV + BED file input
# Uses make_synthetic_bed_csv() and expect_regression_match() from
# helper-synthetic-builders.R

library(testthat)
library(epiRomics)
library(GenomicRanges)
library(IRanges)

# ============================================================================
# Test 1: build_dB regression -- synthetic CSV + BED files
# ============================================================================
test_that("build_dB regression: synthetic CSV + BED", {
  # Create synthetic BED entries: 1 histone, 1 chip, 1 functional
  histone_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = base::c(1000L, 5000L, 10000L),
      end = base::c(2000L, 6000L, 11000L)
    )
  )

  chip_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = base::c(1000L, 5000L),
      end = base::c(2000L, 6000L)
    )
  )

  functional_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = base::c(1000L),
      end = base::c(2000L)
    )
  )

  bed_entries <- base::list(
    base::list(name = "h3k4me1", type = "histone", content = histone_gr),
    base::list(name = "TF1", type = "chip", content = chip_gr),
    base::list(name = "fantom", type = "functional", content = functional_gr)
  )

  csv_path <- make_synthetic_bed_csv(bed_entries, genome = "hg38")
  base::on.exit({
    # Clean up CSV and BED files
    manifest <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
    base::file.remove(csv_path)
    for (p in manifest$path) {
      if (base::file.exists(p)) base::file.remove(p)
    }
  }, add = TRUE)

  # Call epiRomics_build_dB with synthetic data
  # This will call annotatr::read_annotations + annotatr::build_annotations
  result <- base::suppressMessages(base::suppressWarnings(
    epiRomics_build_dB(
      epiRomics_db_file = csv_path,
      txdb_organism = "TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene",
      epiRomics_genome = "hg38",
      epiRomics_organism = "org.Hs.eg.db"
    )
  ))

  # Verify result is epiRomicsS4
  testthat::expect_s4_class(result, "epiRomicsS4")
  testthat::expect_true(base::length(result@annotations) > 0L)

  # Extract stable representation for regression comparison
  # Avoid saving the full GRanges (may include environment refs)
  # Instead save key structural properties
  annot_df <- base::as.data.frame(result@annotations)

  # Filter to only custom annotations (our synthetic BED data)
  # to avoid builtin annotation variability across annotatr versions
  custom_mask <- base::grepl("custom", annot_df$type)
  custom_df <- annot_df[custom_mask, , drop = FALSE]
  base::rownames(custom_df) <- NULL

  result_snapshot <- base::list(
    custom_annotations_df = custom_df,
    meta = result@meta[, base::c("name", "type", "genome", "format")],
    genome = result@genome,
    txdb = result@txdb,
    organism = result@organism,
    n_total_annotations = base::length(result@annotations),
    n_custom_annotations = base::nrow(custom_df),
    annotation_types = base::sort(base::unique(annot_df$type))
  )

  expect_regression_match(result_snapshot, "build_dB_synthetic")
})
