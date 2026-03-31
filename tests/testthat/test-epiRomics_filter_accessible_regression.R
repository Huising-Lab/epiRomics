# Regression baseline tests for epiRomics_filter_accessible()
# Verifies output stability across 3 modes: genelist, bed, combined
# Signal mode requires real BigWig files and is skipped here (covered in plan 05)
# Uses expect_regression_match() from helper-synthetic-builders.R

library(testthat)
library(epiRomics)
library(GenomicRanges)
library(IRanges)

# ============================================================================
# Test 1: Genelist mode regression baseline
# ============================================================================
test_that("filter_accessible regression: genelist mode", {
  pe <- make_synthetic_putative_enhancers(5L)

  result <- base::suppressMessages(
    epiRomics_filter_accessible(
      pe,
      mode = "genelist",
      gene_list = base::c("GENE1", "GENE2")
    )
  )

  expect_regression_match(result, "filter_accessible_genelist")
})

# ============================================================================
# Test 2: BED mode regression baseline
# ============================================================================
test_that("filter_accessible regression: bed mode", {
  pe <- make_synthetic_putative_enhancers(5L)

  # Create a temp BED overlapping PE rows 1 and 3
  # PE starts: 1000, 6000, 11000, 16000, 21000
  bed_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = base::c(900L, 10900L),
      end = base::c(1500L, 11500L)
    )
  )
  bed_file <- base::tempfile(fileext = ".bed")
  base::on.exit(base::file.remove(bed_file), add = TRUE)
  rtracklayer::export(bed_gr, bed_file, format = "BED")

  result <- base::suppressMessages(
    epiRomics_filter_accessible(
      pe,
      mode = "bed",
      bed_path = bed_file
    )
  )

  expect_regression_match(result, "filter_accessible_bed")
})

# ============================================================================
# Test 3: Combined mode regression baseline (genelist + bed)
# ============================================================================
test_that("filter_accessible regression: combined mode", {
  pe <- make_synthetic_putative_enhancers(5L)

  # BED overlaps PE row 1 (start=1000)
  bed_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = 900L,
      end = 1500L
    )
  )
  bed_file <- base::tempfile(fileext = ".bed")
  base::on.exit(base::file.remove(bed_file), add = TRUE)
  rtracklayer::export(bed_gr, bed_file, format = "BED")

  # Genelist matches PE row 4 (GENE4)
  result <- base::suppressMessages(
    epiRomics_filter_accessible(
      pe,
      mode = "combined",
      bed_path = bed_file,
      gene_list = base::c("GENE4")
    )
  )

  expect_regression_match(result, "filter_accessible_combined")
})
