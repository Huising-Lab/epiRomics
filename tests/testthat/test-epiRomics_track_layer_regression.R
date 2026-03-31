# ============================================================================
# REGRESSION BASELINE TESTS for epiRomics_track_layer, _track_layer_fast,
# and epiRomics_track_layer_gene
# Uses synthetic BigWig files and synthetic dB objects for testing.
# Baselines capture either successful return values or deterministic error
# patterns -- both are valid regression detection.
# Baselines saved as .rds in testdata/ via expect_regression_match().
# ============================================================================

library(testthat)
library(epiRomics)
library(GenomicRanges)
library(IRanges)

test_that("track_layer regression: synthetic data", {
  testthat::skip_if_not_installed("rtracklayer")

  # Create synthetic dB
  dB <- make_synthetic_dB_full()

  # Create synthetic BigWig files
  gr_bw <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = base::seq(1L, 9001L, 1000L),
      end = base::seq(1000L, 10000L, 1000L)
    )
  )
  scores1 <- base::c(5, 10, 50, 80, 90, 5, 10, 95, 5, 10)
  scores2 <- base::c(10, 5, 5, 5, 5, 50, 80, 90, 95, 5)

  bw1 <- make_synthetic_bigwig(gr_bw, scores1)
  bw2 <- make_synthetic_bigwig(gr_bw, scores2)

  pdf_path <- base::tempfile(fileext = ".pdf")
  base::on.exit({
    if (grDevices::dev.cur() > 1L) grDevices::dev.off()
    if (base::file.exists(bw1)) base::file.remove(bw1)
    if (base::file.exists(bw2)) base::file.remove(bw2)
    if (base::file.exists(pdf_path)) base::file.remove(pdf_path)
  }, add = TRUE)

  bw_paths <- base::c(Alpha_ATAC = bw1, Beta_ATAC = bw2)

  # Create track_connection data.frame
  tc <- base::data.frame(
    path = base::unname(bw_paths),
    name = base::names(bw_paths),
    color = base::c("#FF0000", "#0000FF"),
    type = base::c("atac", "atac"),
    stringsAsFactors = FALSE
  )

  # Create a synthetic enhanceosome-like object
  # Needs: annotations with SYMBOL, geneStart, geneEnd, geneChr
  synth_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 5000L, end = 6000L)
  )
  synth_gr$SYMBOL <- "SyntheticGene"
  synth_gr$geneStart <- 1000L
  synth_gr$geneEnd <- 10000L
  synth_gr$geneChr <- "1"
  base::names(synth_gr) <- "PE_1"

  pe_dB <- methods::new("epiRomicsS4")
  pe_dB@annotations <- synth_gr
  pe_dB@meta <- dB@meta
  pe_dB@txdb <- dB@txdb
  pe_dB@organism <- dB@organism
  pe_dB@genome <- dB@genome

  # Suppress plot output
  grDevices::pdf(pdf_path)

  outcome <- base::tryCatch(
    {
      result <- epiRomics_track_layer(
        epiRomics_putative_enhanceosome = pe_dB,
        epiRomics_index = 1L,
        epiRomics_dB = dB,
        epiRomics_track_connection = tc,
        chromatin_states = NULL,
        show_bigwig = TRUE,
        show_chromatin = FALSE,
        show_annotations = FALSE,
        show_gene_model = TRUE,
        show_enhancer_highlight = TRUE
      )
      # If we get here, function succeeded
      base::list(status = "success", result_class = base::class(result))
    },
    error = function(e) {
      base::list(
        status = "error",
        error_class = base::class(e),
        error_msg_pattern = base::sub(".*: ", "", e$message)
      )
    }
  )

  expect_regression_match(outcome, "track_layer_default")
})

test_that("track_layer_fast is alias for track_layer", {
  # Verify that epiRomics_track_layer_fast is identical to epiRomics_track_layer
  is_alias <- base::identical(epiRomics_track_layer_fast, epiRomics_track_layer)

  alias_result <- base::list(is_alias = is_alias)

  expect_regression_match(alias_result, "track_layer_fast_alias")
})

test_that("track_layer_gene regression: synthetic data", {
  testthat::skip_if_not_installed("rtracklayer")

  # Create synthetic dB
  dB <- make_synthetic_dB_full()

  # Create synthetic BigWig files
  gr_bw <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = base::seq(1L, 9001L, 1000L),
      end = base::seq(1000L, 10000L, 1000L)
    )
  )
  scores1 <- base::c(5, 10, 50, 80, 90, 5, 10, 95, 5, 10)

  bw1 <- make_synthetic_bigwig(gr_bw, scores1)
  pdf_path <- base::tempfile(fileext = ".pdf")
  base::on.exit({
    if (grDevices::dev.cur() > 1L) grDevices::dev.off()
    if (base::file.exists(bw1)) base::file.remove(bw1)
    if (base::file.exists(pdf_path)) base::file.remove(pdf_path)
  }, add = TRUE)

  bw_paths <- base::c(Alpha_ATAC = bw1)

  tc <- base::data.frame(
    path = base::unname(bw_paths),
    name = base::names(bw_paths),
    color = base::c("#FF0000"),
    type = base::c("atac"),
    stringsAsFactors = FALSE
  )

  # Suppress plot output
  grDevices::pdf(pdf_path)

  # Use a real gene symbol that exists in hg38 TxDb.
  # The synthetic BigWig only has chr1 data so warnings about chromosome
  # mismatch are expected and suppressed.
  outcome <- base::suppressWarnings(base::tryCatch(
    {
      result <- epiRomics_track_layer_gene(
        gene_symbol = "TP53",
        epiRomics_dB = dB,
        epiRomics_track_connection = tc,
        chromatin_states = NULL,
        show_bigwig = TRUE,
        show_chromatin = FALSE,
        show_annotations = FALSE,
        show_gene_model = TRUE,
        show_enhancer_highlight = FALSE
      )
      base::list(status = "success", result_class = base::class(result))
    },
    error = function(e) {
      base::list(
        status = "error",
        error_class = base::class(e),
        error_msg_pattern = base::sub(".*: ", "", e$message)
      )
    }
  ))

  expect_regression_match(outcome, "track_layer_gene_default")
})
