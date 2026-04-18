# ============================================================================
# Shared synthetic data builders for regression testing
# Auto-sourced by testthat before all tests (helper-* convention)
#
# Provides reusable constructors for:
#   - epiRomicsS4 objects with histone/chip/functional annotations
#   - Temporary BigWig files for coverage testing
#   - CSV manifest files for build_database testing
#   - Putative enhancer data.frames for annotation testing
#   - Regression baseline save/load/compare utilities
# ============================================================================

# ----------------------------------------------------------------------------
# 1. make_synthetic_dB: Generic epiRomicsS4 builder
# ----------------------------------------------------------------------------
#' Build a synthetic epiRomicsS4 object from named lists of GRanges
#'
#' @param marks_list Named list of GRanges for histone marks
#' @param genome Character genome assembly name (default "hg38")
#' @param chip_list Optional named list of GRanges for ChIP TF peaks
#' @param functional_list Optional named list of GRanges for functional annotations
#' @return A valid epiRomicsS4 object
make_synthetic_dB <- function(marks_list, genome = "hg38",
                              chip_list = NULL, functional_list = NULL) {
  all_gr <- GenomicRanges::GRanges()
  meta_rows <- base::list()

  # Process histone marks
  if (base::length(marks_list) > 0) {
    for (mark_name in base::names(marks_list)) {
      gr <- marks_list[[mark_name]]
      gr$type <- base::paste0(genome, "_custom_", mark_name)
      all_gr <- base::c(all_gr, gr)
      meta_rows[[mark_name]] <- base::data.frame(
        name = mark_name,
        type = "histone",
        stringsAsFactors = FALSE
      )
    }
  }

  # Process ChIP TF peaks
  if (!base::is.null(chip_list) && base::length(chip_list) > 0) {
    for (chip_name in base::names(chip_list)) {
      gr <- chip_list[[chip_name]]
      gr$type <- base::paste0(genome, "_custom_", chip_name)
      all_gr <- base::c(all_gr, gr)
      meta_rows[[chip_name]] <- base::data.frame(
        name = chip_name,
        type = "chip",
        stringsAsFactors = FALSE
      )
    }
  }

  # Process functional annotations
  if (!base::is.null(functional_list) && base::length(functional_list) > 0) {
    for (func_name in base::names(functional_list)) {
      gr <- functional_list[[func_name]]
      gr$type <- base::paste0(genome, "_custom_", func_name)
      all_gr <- base::c(all_gr, gr)
      meta_rows[[func_name]] <- base::data.frame(
        name = func_name,
        type = "functional",
        stringsAsFactors = FALSE
      )
    }
  }

  meta <- base::do.call(base::rbind, meta_rows)

  methods::new("epiRomicsS4",
    annotations = all_gr,
    meta = meta,
    genome = genome,
    txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene",
    organism = "org.Hs.eg.db"
  )
}

# ----------------------------------------------------------------------------
# 2. make_synthetic_dB_histone_only: Pre-built with 3 histone marks
# ----------------------------------------------------------------------------
#' Build a pre-configured epiRomicsS4 with 3 histone marks on chr1
#'
#' Marks: h3k4me1, h3k27ac, h3k27me3 with overlapping and non-overlapping regions.
#' Good for: chromatin_states, enhancers_co_marks, putative_enhancers
#'
#' @param genome Character genome assembly name (default "hg38")
#' @return A valid epiRomicsS4 object with 3 histone marks
make_synthetic_dB_histone_only <- function(genome = "hg38") {
  marks <- base::list(
    h3k4me1 = GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges = IRanges::IRanges(
        start = base::c(1000L, 5000L, 10000L, 20000L),
        end = base::c(2000L, 6000L, 11000L, 21000L)
      )
    ),
    h3k27ac = GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges = IRanges::IRanges(
        start = base::c(1000L, 5000L, 30000L),
        end = base::c(2000L, 6000L, 31000L)
      )
    ),
    h3k27me3 = GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges = IRanges::IRanges(
        start = base::c(10000L, 40000L),
        end = base::c(11000L, 41000L)
      )
    )
  )
  make_synthetic_dB(marks, genome = genome)
}

# ----------------------------------------------------------------------------
# 3. make_synthetic_dB_full: Pre-built with histones + ChIP TFs + functional
# ----------------------------------------------------------------------------
#' Build a full synthetic epiRomicsS4 with histones, ChIP TFs, and functional data
#'
#' Histones: h3k4me1, h3k27ac, h3k27me3, h3k4me3, h3k36me3
#' ChIP: TF1, TF2
#' Functional: fantom
#' Good for: enhanceosome, tf_cobinding, tf_overlap, enhancers_filter,
#'           enhancer_predictor_to_ref
#'
#' @param genome Character genome assembly name (default "hg38")
#' @return A valid epiRomicsS4 object with full annotation set
make_synthetic_dB_full <- function(genome = "hg38") {
  marks <- base::list(
    h3k4me1 = GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges = IRanges::IRanges(
        start = base::c(1000L, 5000L, 10000L, 20000L, 50000L),
        end = base::c(2000L, 6000L, 11000L, 21000L, 51000L)
      )
    ),
    h3k27ac = GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges = IRanges::IRanges(
        start = base::c(1000L, 5000L, 30000L, 50000L),
        end = base::c(2000L, 6000L, 31000L, 51000L)
      )
    ),
    h3k27me3 = GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges = IRanges::IRanges(
        start = base::c(10000L, 40000L),
        end = base::c(11000L, 41000L)
      )
    ),
    h3k4me3 = GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges = IRanges::IRanges(
        start = base::c(1000L, 50000L, 60000L),
        end = base::c(2000L, 51000L, 61000L)
      )
    ),
    h3k36me3 = GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges = IRanges::IRanges(
        start = base::c(70000L, 80000L),
        end = base::c(71000L, 81000L)
      )
    )
  )

  chip <- base::list(
    TF1 = GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges = IRanges::IRanges(
        start = base::c(1000L, 5000L, 50000L),
        end = base::c(2000L, 6000L, 51000L)
      )
    ),
    TF2 = GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges = IRanges::IRanges(
        start = base::c(1000L, 20000L, 50000L),
        end = base::c(2000L, 21000L, 51000L)
      )
    )
  )

  functional <- base::list(
    fantom = GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges = IRanges::IRanges(
        start = base::c(1000L, 50000L),
        end = base::c(2000L, 51000L)
      )
    )
  )

  make_synthetic_dB(marks, genome = genome, chip_list = chip,
                    functional_list = functional)
}

# ----------------------------------------------------------------------------
# 4. make_synthetic_bigwig: Create a temporary BigWig file
# ----------------------------------------------------------------------------
#' Create a temporary BigWig file from GRanges with scores
#'
#' Requires rtracklayer. Caller should guard with
#' testthat::skip_if_not_installed("rtracklayer").
#'
#' @param regions_gr GRanges with the regions to export
#' @param scores Numeric vector of scores (same length as regions_gr)
#' @param path Optional output path. If NULL, uses tempdir().
#' @return Path to the created BigWig file. Caller must clean up.
make_synthetic_bigwig <- function(regions_gr, scores, path = NULL) {
  if (base::is.null(path)) {
    path <- base::tempfile(fileext = ".bw")
  }

  # Attach scores as the score column required by BigWig format
  regions_gr$score <- scores

  # BigWig export requires seqlengths
  si <- GenomeInfoDb::Seqinfo(
    seqnames = base::as.character(
      base::unique(GenomicRanges::seqnames(regions_gr))
    ),
    seqlengths = base::rep(250000000L,
      base::length(base::unique(GenomicRanges::seqnames(regions_gr)))
    )
  )
  GenomeInfoDb::seqinfo(regions_gr) <- si

  rtracklayer::export(regions_gr, path, format = "BigWig")
  base::return(path)
}

# ----------------------------------------------------------------------------
# 5. make_synthetic_bed_csv: Create a CSV manifest for build_database
# ----------------------------------------------------------------------------
#' Create a CSV manifest file suitable for build_database
#'
#' @param bed_entries List of lists, each with: name (character), type (character),
#'   content (GRanges). Writes BED files to tempdir and creates the CSV manifest.
#' @param genome Character genome assembly name (default "hg38")
#' @param path Optional output CSV path. If NULL, uses tempdir().
#' @return Path to the created CSV manifest file. Caller must clean up.
make_synthetic_bed_csv <- function(bed_entries, genome = "hg38", path = NULL) {
  if (base::is.null(path)) {
    path <- base::tempfile(fileext = ".csv")
  }

  manifest_rows <- base::list()
  temp_dir <- base::tempdir()

  for (i in base::seq_along(bed_entries)) {
    entry <- bed_entries[[i]]
    bed_path <- base::file.path(temp_dir,
      base::paste0(entry$name, ".bed"))
    rtracklayer::export(entry$content, bed_path, format = "BED")

    manifest_rows[[i]] <- base::data.frame(
      path = bed_path,
      genome = genome,
      name = entry$name,
      format = "bed",
      type = entry$type,
      stringsAsFactors = FALSE
    )
  }

  manifest <- base::do.call(base::rbind, manifest_rows)
  utils::write.csv(manifest, path, row.names = FALSE)
  base::return(path)
}

# ----------------------------------------------------------------------------
# 6. make_synthetic_putative_enhancers: Mock putative enhancers data.frame
# ----------------------------------------------------------------------------
#' Create a data.frame matching the output format of find_putative_enhancers
#'
#' @param n Integer number of putative enhancer rows to generate (default 5)
#' @return A data.frame with columns matching find_putative_enhancers output
make_synthetic_putative_enhancers <- function(n = 5L) {
  starts <- base::seq(1000L, by = 5000L, length.out = n)
  ends <- starts + 1000L

  base::data.frame(
    putative_id = base::paste0("PE_", base::seq_len(n)),
    chr = base::rep("chr1", n),
    start = starts,
    end = ends,
    width = ends - starts,
    source = base::rep("co_marks", n),
    chromatin_state = base::rep("active", n),
    chromatin_state_detail = base::rep("H3K4me1+H3K27ac", n),
    histone_marks = base::rep("h3k4me1,h3k27ac", n),
    n_histone_marks = base::rep(2L, n),
    h2az = base::rep(FALSE, n),
    tf_names = base::rep("", n),
    n_tfs = base::rep(0L, n),
    SYMBOL = base::paste0("GENE", base::seq_len(n)),
    stringsAsFactors = FALSE
  )
}

# ----------------------------------------------------------------------------
# 7. save_regression_baseline: Save .rds baseline
# ----------------------------------------------------------------------------
#' Save an R object as an .rds regression baseline
#'
#' @param result The R object to save
#' @param test_name Character name for the baseline (used as filename stem)
#' @param subdir Character subdirectory under tests/testthat/ (default "testdata")
save_regression_baseline <- function(result, test_name, subdir = "testdata") {
  dir_path <- testthat::test_path(subdir)
  if (!base::dir.exists(dir_path)) {
    base::dir.create(dir_path, recursive = TRUE)
  }
  baseline_path <- base::file.path(dir_path, base::paste0(test_name, ".rds"))
  base::saveRDS(result, baseline_path)
  base::invisible(baseline_path)
}

# ----------------------------------------------------------------------------
# 8. load_regression_baseline: Load .rds baseline
# ----------------------------------------------------------------------------
#' Load a regression baseline .rds file
#'
#' @param test_name Character name of the baseline
#' @param subdir Character subdirectory under tests/testthat/ (default "testdata")
#' @return The loaded R object, or NULL if the baseline file does not exist
load_regression_baseline <- function(test_name, subdir = "testdata") {
  baseline_path <- base::file.path(
    testthat::test_path(subdir),
    base::paste0(test_name, ".rds")
  )
  if (!base::file.exists(baseline_path)) {
    base::return(NULL)
  }
  base::readRDS(baseline_path)
}

# ----------------------------------------------------------------------------
# 9. expect_regression_match: Compare current output against saved baseline
# ----------------------------------------------------------------------------
#' Compare current function output against a saved regression baseline
#'
#' On first run (no baseline exists), saves current output as the baseline and
#' passes with an informational message. On subsequent runs, compares current
#' output against the saved baseline using testthat::expect_equal().
#'
#' @param current The current function output to compare
#' @param test_name Character name for the baseline
#' @param subdir Character subdirectory under tests/testthat/ (default "testdata")
#' @param tolerance Optional numeric tolerance for floating-point comparison
#' @return Invisibly returns the current output
expect_regression_match <- function(current, test_name, subdir = "testdata",
                                    tolerance = NULL) {
  baseline_path <- base::file.path(
    testthat::test_path(subdir),
    base::paste0(test_name, ".rds")
  )

  if (!base::file.exists(baseline_path)) {
    # First run: save baseline
    dir_path <- base::dirname(baseline_path)
    if (!base::dir.exists(dir_path)) {
      base::dir.create(dir_path, recursive = TRUE)
    }
    base::saveRDS(current, baseline_path)
    testthat::succeed(base::paste("Baseline saved:", baseline_path))
    base::return(base::invisible(current))
  }

  baseline <- base::readRDS(baseline_path)
  if (!base::is.null(tolerance)) {
    testthat::expect_equal(current, baseline, tolerance = tolerance)
  } else {
    testthat::expect_equal(current, baseline)
  }
  base::invisible(current)
}
