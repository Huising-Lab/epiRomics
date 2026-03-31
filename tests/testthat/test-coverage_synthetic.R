## Coverage tests for exported functions that lack dedicated test files
## All tests are synthetic (no real data, no database rebuild)

# ---- Helper: minimal epiRomicsS4 for validation tests ----
.make_empty_db <- function(genome = "hg38") {
  epiRomicsS4(
    annotations = GenomicRanges::GRanges(),
    meta = base::data.frame(
      name = base::character(),
      type = base::character(),
      file = base::character(),
      stringsAsFactors = FALSE
    ),
    genome = genome
  )
}

# Helper: DB with a single dummy annotation (passes validate_epiRomics_dB)
.make_minimal_db <- function(genome = "hg38") {
  gr <- GenomicRanges::GRanges("chr1:1-100")
  gr$type <- paste0(genome, "_custom_dummy")
  epiRomicsS4(
    annotations = gr,
    meta = base::data.frame(
      name = "dummy",
      type = "histone",
      file = "dummy.bed",
      stringsAsFactors = FALSE
    ),
    genome = genome
  )
}

# Helper: DB with chip data in meta (passes "No ChIP" check)
.make_chip_db <- function(genome = "hg38") {
  gr <- GenomicRanges::GRanges(c("chr1:1-100", "chr1:200-300"))
  gr$type <- c(paste0(genome, "_custom_TF1"), paste0(genome, "_custom_TF2"))
  epiRomicsS4(
    annotations = gr,
    meta = base::data.frame(
      name = c("TF1", "TF2"),
      type = c("chip", "chip"),
      file = c("tf1.bed", "tf2.bed"),
      stringsAsFactors = FALSE
    ),
    genome = genome
  )
}

# ---- 1. epiRomics_track_layer_gene: validation ----
test_that("epiRomics_track_layer_gene rejects non-character gene_symbol", {
  expect_error(
    epiRomics_track_layer_gene(gene_symbol = 123, epiRomics_dB = .make_empty_db(),
      epiRomics_track_connection = data.frame()),
    "gene_symbol"
  )
})

test_that("epiRomics_track_layer_gene rejects vector gene_symbol", {
  expect_error(
    epiRomics_track_layer_gene(gene_symbol = c("INS", "GCG"),
      epiRomics_dB = .make_empty_db(),
      epiRomics_track_connection = data.frame()),
    "single"
  )
})

# ---- 4. epiRomics_enhancer_predictor_to_ref: validation ----
test_that("epiRomics_enhancer_predictor_to_ref validates dB", {
  expect_error(
    epiRomics_enhancer_predictor_to_ref(list()),
    "epiRomics_dB|epiRomicsS4"
  )
})

test_that("epiRomics_enhancer_predictor_to_ref validates dB annotations before params", {
  db <- .make_empty_db()
  # Empty annotations error fires before curated_database validation
  expect_error(
    epiRomics_enhancer_predictor_to_ref(db, epiRomics_curated_database = 123),
    "empty|NULL|annotations"
  )
})

# ---- 5. epiRomics_enhancers_co_marks: validation ----
test_that("epiRomics_enhancers_co_marks validates dB", {
  expect_error(
    epiRomics_enhancers_co_marks(list()),
    "epiRomics_dB|epiRomicsS4"
  )
})

test_that("epiRomics_enhancers_co_marks validates dB annotations before params", {
  db <- .make_empty_db()
  # Empty annotations error fires before histone param validation
  expect_error(
    epiRomics_enhancers_co_marks(db, epiRomics_histone_mark_1 = 123),
    "empty|NULL|annotations"
  )
})

# ---- 6. maxCovFilesCached: validation ----
test_that("maxCovFilesCached validates file paths", {
  gr <- GenomicRanges::GRanges("chr1:1-100")
  expect_error(
    maxCovFilesCached(c("/nonexistent/a.bw"), gr),
    "not found|exist|valid"
  )
})

test_that("maxCovFilesCached validates file existence before GRanges", {
  # File path validation fires before GRanges validation
  expect_error(
    maxCovFilesCached("dummy.bw", "not_a_granges"),
    "not exist|not found"
  )
})

# ---- 7. classify_celltype_accessibility: validation ----
test_that("classify_celltype_accessibility requires named bw_paths", {
  gr <- GenomicRanges::GRanges("chr1:1-100")
  expect_error(
    classify_celltype_accessibility(c("a.bw", "b.bw"), gr),
    "named"
  )
})

test_that("classify_celltype_accessibility rejects non-character input", {
  gr <- GenomicRanges::GRanges("chr1:1-100")
  expect_error(
    classify_celltype_accessibility(123, gr),
    "character|named"
  )
})

# ---- 8. epiRomics_annotate_putative: validation ----
test_that("epiRomics_annotate_putative requires non-empty data.frame", {
  db <- .make_empty_db()
  expect_error(
    epiRomics_annotate_putative(data.frame(), db),
    "non-empty"
  )
})

test_that("epiRomics_annotate_putative rejects non-data.frame", {
  db <- .make_empty_db()
  expect_error(
    epiRomics_annotate_putative("not_a_df", db),
    "data.frame"
  )
})

# ---- 9. epiRomics_tf_overlap: validation ----
test_that("epiRomics_tf_overlap validates dB class", {
  db <- .make_empty_db()
  expect_error(
    epiRomics_tf_overlap(list(), db),
    "epiRomics_dB|epiRomicsS4"
  )
  expect_error(
    epiRomics_tf_overlap(db, list()),
    "epiRomics_dB|epiRomicsS4"
  )
})

test_that("epiRomics_tf_overlap errors when no ChIP data in meta", {
  # Use DB that has annotations (passes validate_epiRomics_dB) but no chip meta
  db <- .make_minimal_db()
  expect_error(
    epiRomics_tf_overlap(db, db),
    "No ChIP"
  )
})

test_that("epiRomics_tf_overlap rejects non-GRanges regions", {
  # Need DB with chip meta to get past "No ChIP" check
  db <- .make_chip_db()
  expect_error(
    epiRomics_tf_overlap(db, db, regions = "not_granges"),
    "GRanges"
  )
})

test_that("epiRomics_tf_overlap works with synthetic chip data", {
  db <- .make_chip_db()
  result <- epiRomics_tf_overlap(db, db)
  expect_type(result, "list")
  expect_true("overlap_matrix" %in% names(result))
  expect_true("jaccard_matrix" %in% names(result))
  expect_true("summary" %in% names(result))
  expect_equal(length(result$tf_names), 2)
  # Overlap matrix should be 2x2
  expect_equal(dim(result$overlap_matrix), c(2, 2))
})

# ---- 10. epiRomics_chromatin_states_categories: validation ----
test_that("epiRomics_chromatin_states_categories validates empty annotations", {
  db <- .make_empty_db()
  expect_error(
    epiRomics_chromatin_states_categories(db),
    "empty|NULL|annotations"
  )
})

test_that("epiRomics_chromatin_states_categories errors with no histone marks", {
  # DB with annotations but no histone marks in meta -> error
  db <- .make_chip_db()  # only chip in meta, no histone
  expect_error(
    epiRomics_chromatin_states_categories(db),
    "No histone|histone_marks"
  )
})

# ---- 11. plot_signal_histogram: validation ----
test_that("plot_signal_histogram rejects empty bw_paths", {
  gr <- GenomicRanges::GRanges("chr1:1-100")
  expect_error(
    plot_signal_histogram(character(0), gr),
    "character vector"
  )
})

test_that("plot_signal_histogram rejects non-character bw_paths", {
  gr <- GenomicRanges::GRanges("chr1:1-100")
  expect_error(
    plot_signal_histogram(123, gr),
    "character"
  )
})
