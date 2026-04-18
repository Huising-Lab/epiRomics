library(testthat)
library(epiRomics)
suppressPackageStartupMessages(library(GenomicRanges))

# Build database using helper
database <- build_test_dB()

# Generate enhanceosome using real data (guard against NULL when extdata missing)
enhancers <- NULL
enhanceosome <- NULL
enhanceosome_chr11 <- NULL
if (!is.null(database)) {
  histone_marks <- database@meta$name[database@meta$type == "histone"]
  enhancers <- epiRomics::find_enhancers_by_comarks(database,
    histone_marks[1], histone_marks[2])
  enhanceosome <- epiRomics::find_enhanceosomes(enhancers, database)
  # Restrict to the exact toy BigWig coverage window (hg38 chr11 INS locus,
  # ~2.10–2.22 Mb). Annotations outside this span produce harmless "Failed to
  # summarize range" warnings from rtracklayer.
  toy_window <- GRanges("chr11", IRanges(2099779, 2221221))
  enhanceosome_chr11 <- IRanges::subsetByOverlaps(
    enhanceosome@annotations, toy_window)
}

# --- call_accessible_regions ---
# (plot_super_enhancers, plot_signal_enrichment, plot_tf_chromatin_state,
#  plot_tf_chromatin_celltype are not in stable — see experimental fork)

test_that("call_accessible_regions returns logical vector", {
  require_extdata()
  bw_file <- system.file("extdata", "toy", "BigWigs",
    "Beta_ATAC.toy.bw", package = "epiRomics")
  expect_true(nzchar(bw_file) && file.exists(bw_file),
    info = "Toy Beta_ATAC.toy.bw must ship with the package")

  test_regions <- enhanceosome_chr11[seq_len(min(50, length(enhanceosome_chr11)))]
  expect_gt(length(test_regions), 0L,
    label = "toy enhanceosome must contain chr11 annotations")
  result <- call_accessible_regions(bw_file, test_regions, z_threshold = 2.0)

  expect_true(is.logical(result))
  expect_equal(length(result), length(test_regions))
})

test_that("call_accessible_regions with strict threshold returns fewer TRUE", {
  require_extdata()
  bw_file <- system.file("extdata", "toy", "BigWigs",
    "Beta_ATAC.toy.bw", package = "epiRomics")
  expect_true(nzchar(bw_file) && file.exists(bw_file),
    info = "Toy Beta_ATAC.toy.bw must ship with the package")

  test_regions <- enhanceosome_chr11[seq_len(min(50, length(enhanceosome_chr11)))]
  expect_gt(length(test_regions), 0L,
    label = "toy enhanceosome must contain chr11 annotations")
  result_loose <- call_accessible_regions(bw_file, test_regions, z_threshold = 1.0)
  result_strict <- call_accessible_regions(bw_file, test_regions, z_threshold = 3.0)

  expect_true(sum(result_strict) <= sum(result_loose))
})

test_that("call_accessible_regions errors on invalid inputs", {
  require_extdata()
  expect_gte(length(enhanceosome_chr11), 5L,
    label = "toy enhanceosome chr11 subset must contain at least 5 annotations")
  expect_error(call_accessible_regions("nonexistent.bw",
    enhanceosome_chr11[1:5]))
  bw_file <- system.file("extdata", "toy", "BigWigs",
    "Beta_ATAC.toy.bw", package = "epiRomics")
  expect_true(nzchar(bw_file) && file.exists(bw_file),
    info = "Toy Beta_ATAC.toy.bw must ship with the package")
  expect_error(call_accessible_regions(bw_file, "not_a_granges"))
})
