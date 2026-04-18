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

# Helper: DB with a single dummy annotation (passes validate_database)
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

# ---- 1. plot_gene_tracks: validation ----
test_that("plot_gene_tracks rejects non-character gene_symbol", {
  expect_error(
    plot_gene_tracks(gene_symbol = 123, database = .make_empty_db(),
      track_connection = data.frame()),
    "gene_symbol"
  )
})

test_that("plot_gene_tracks rejects vector gene_symbol", {
  expect_error(
    plot_gene_tracks(gene_symbol = c("INS", "GCG"),
      database = .make_empty_db(),
      track_connection = data.frame()),
    "single"
  )
})

# ---- 4. benchmark_enhancer_predictor: validation ----
test_that("benchmark_enhancer_predictor validates dB", {
  expect_error(
    benchmark_enhancer_predictor(list()),
    "database|epiRomicsS4"
  )
})

test_that("benchmark_enhancer_predictor validates dB annotations before params", {
  db <- .make_empty_db()
  # Empty annotations error fires before curated_database validation
  expect_error(
    benchmark_enhancer_predictor(db, curated_database = 123),
    "empty|NULL|annotations"
  )
})

# ---- 5. find_enhancers_by_comarks: validation ----
test_that("find_enhancers_by_comarks validates dB", {
  expect_error(
    find_enhancers_by_comarks(list()),
    "database|epiRomicsS4"
  )
})

test_that("find_enhancers_by_comarks validates dB annotations before params", {
  db <- .make_empty_db()
  # Empty annotations error fires before histone param validation
  expect_error(
    find_enhancers_by_comarks(db, histone_mark_1 = 123),
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

# ---- 8. annotate_enhancers: validation ----
test_that("annotate_enhancers requires non-empty data.frame", {
  db <- .make_empty_db()
  expect_error(
    annotate_enhancers(data.frame(), db),
    "non-empty"
  )
})

test_that("annotate_enhancers rejects non-data.frame", {
  db <- .make_empty_db()
  expect_error(
    annotate_enhancers("not_a_df", db),
    "data.frame"
  )
})

# ---- 9. analyze_tf_overlap: validation ----
test_that("analyze_tf_overlap validates dB class", {
  db <- .make_empty_db()
  expect_error(
    analyze_tf_overlap(list(), db),
    "database|epiRomicsS4"
  )
  expect_error(
    analyze_tf_overlap(db, list()),
    "database|epiRomicsS4"
  )
})

test_that("analyze_tf_overlap errors when no ChIP data in meta", {
  # Use DB that has annotations (passes validate_database) but no chip meta
  db <- .make_minimal_db()
  expect_error(
    analyze_tf_overlap(db, db),
    "No ChIP"
  )
})

test_that("analyze_tf_overlap rejects non-GRanges regions", {
  # Need DB with chip meta to get past "No ChIP" check
  db <- .make_chip_db()
  expect_error(
    analyze_tf_overlap(db, db, regions = "not_granges"),
    "GRanges"
  )
})

test_that("analyze_tf_overlap works with synthetic chip data", {
  db <- .make_chip_db()
  result <- analyze_tf_overlap(db, db)
  expect_type(result, "list")
  expect_true("overlap_matrix" %in% names(result))
  expect_true("jaccard_matrix" %in% names(result))
  expect_true("summary" %in% names(result))
  expect_equal(length(result$tf_names), 2)
  # Overlap matrix should be 2x2
  expect_equal(dim(result$overlap_matrix), c(2, 2))
})

# ---- 10. chromatin_state_categories: validation ----
test_that("chromatin_state_categories validates empty annotations", {
  db <- .make_empty_db()
  expect_error(
    chromatin_state_categories(db),
    "empty|NULL|annotations"
  )
})

test_that("chromatin_state_categories errors with no histone marks", {
  # DB with annotations but no histone marks in meta -> error
  db <- .make_chip_db()  # only chip in meta, no histone
  expect_error(
    chromatin_state_categories(db),
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
