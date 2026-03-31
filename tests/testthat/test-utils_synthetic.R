# ============================================================================
# SYNTHETIC DATA TESTS for utils.R internal functions
# Tests: validation functions, chromatin state palette, signal aggregation,
#        resolve_txdb, aggregate_weighted_signal
# ============================================================================

library(testthat)
library(epiRomics)
library(GenomicRanges)
library(IRanges)

# ============================================================================
# Helper to create valid mock epiRomicsS4 objects
# ============================================================================
make_valid_dB <- function(genome = "hg38") {
  gr <- GRanges("chr1", IRanges(1000, 2000))
  gr$type <- paste0(genome, "_custom_h3k27ac")
  meta <- data.frame(name = "h3k27ac", type = "histone",
    stringsAsFactors = FALSE)
  methods::new("epiRomicsS4", annotations = gr, meta = meta,
    genome = genome,
    txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene",
    organism = "org.Hs.eg.db")
}

# ============================================================================
# validate_epiRomics_dB tests
# ============================================================================

test_that("validate_epiRomics_dB accepts valid objects", {
  db <- make_valid_dB()
  expect_true(epiRomics:::validate_epiRomics_dB(db))
})

test_that("validate_epiRomics_dB rejects non-epiRomicsS4 objects", {
  expect_error(epiRomics:::validate_epiRomics_dB("not_a_db"),
    "must be an epiRomicsS4 object")
  expect_error(epiRomics:::validate_epiRomics_dB(42),
    "must be an epiRomicsS4 object")
  expect_error(epiRomics:::validate_epiRomics_dB(NULL),
    "must be an epiRomicsS4 object")
  expect_error(epiRomics:::validate_epiRomics_dB(list()),
    "must be an epiRomicsS4 object")
})

test_that("validate_epiRomics_dB rejects empty annotations", {
  db <- make_valid_dB()
  db@annotations <- GRanges()
  expect_error(epiRomics:::validate_epiRomics_dB(db),
    "annotations is empty or NULL")
})

test_that("validate_epiRomics_dB rejects empty meta", {
  db <- make_valid_dB()
  db@meta <- data.frame()
  expect_error(epiRomics:::validate_epiRomics_dB(db),
    "meta is empty or NULL")
})

test_that("validate_epiRomics_dB rejects empty genome", {
  db <- make_valid_dB()
  db@genome <- character(0)
  expect_error(epiRomics:::validate_epiRomics_dB(db),
    "genome is empty or NULL")
})

test_that("validate_epiRomics_dB includes function_name in error", {
  expect_error(
    epiRomics:::validate_epiRomics_dB("bad", function_name = "my_func"),
    "my_func")
})

# ============================================================================
# validate_genomic_ranges tests
# ============================================================================

test_that("validate_genomic_ranges accepts valid GRanges", {
  gr <- GRanges("chr1", IRanges(1, 100))
  expect_true(epiRomics:::validate_genomic_ranges(gr))
})

test_that("validate_genomic_ranges rejects non-GRanges", {
  expect_error(epiRomics:::validate_genomic_ranges("chr1:1-100"),
    "must be a GRanges object")
  expect_error(epiRomics:::validate_genomic_ranges(42),
    "must be a GRanges object")
  expect_error(epiRomics:::validate_genomic_ranges(data.frame(chr = "chr1")),
    "must be a GRanges object")
})

test_that("validate_genomic_ranges rejects empty GRanges", {
  expect_error(epiRomics:::validate_genomic_ranges(GRanges()),
    "gr is empty")
})

test_that("validate_genomic_ranges includes function_name in error", {
  expect_error(
    epiRomics:::validate_genomic_ranges("bad", function_name = "test_fn"),
    "test_fn")
})

# ============================================================================
# validate_file_paths tests
# ============================================================================

test_that("validate_file_paths accepts existing files", {
  tf <- tempfile()
  writeLines("test", tf)
  on.exit(unlink(tf))
  expect_true(epiRomics:::validate_file_paths(tf))
})

test_that("validate_file_paths rejects non-character", {
  expect_error(epiRomics:::validate_file_paths(42),
    "must be a character vector")
  expect_error(epiRomics:::validate_file_paths(NULL),
    "must be a character vector")
})

test_that("validate_file_paths rejects empty vector", {
  expect_error(epiRomics:::validate_file_paths(character(0)),
    "is empty")
})

test_that("validate_file_paths rejects missing files", {
  expect_error(epiRomics:::validate_file_paths("/nonexistent/file.txt"),
    "do not exist")
})

test_that("validate_file_paths skips existence check when told to", {
  expect_true(epiRomics:::validate_file_paths("/any/path.txt",
    check_exists = FALSE))
})

test_that("validate_file_paths includes function_name in error", {
  expect_error(
    epiRomics:::validate_file_paths(42, function_name = "my_fn"),
    "my_fn")
})

# ============================================================================
# validate_character_param tests
# ============================================================================

test_that("validate_character_param accepts valid strings", {
  expect_true(epiRomics:::validate_character_param("hello", "param1"))
})

test_that("validate_character_param rejects non-character", {
  expect_error(epiRomics:::validate_character_param(42, "param1"),
    "must be a character string")
})

test_that("validate_character_param rejects vector of length > 1", {
  expect_error(epiRomics:::validate_character_param(c("a", "b"), "param1"),
    "must be a single character string")
})

test_that("validate_character_param rejects empty string by default", {
  expect_error(epiRomics:::validate_character_param("", "param1"),
    "cannot be empty or NA")
})

test_that("validate_character_param rejects NA by default", {
  expect_error(epiRomics:::validate_character_param(NA_character_, "param1"),
    "cannot be empty or NA")
})

test_that("validate_character_param allows empty when told to", {
  expect_true(epiRomics:::validate_character_param("", "param1",
    allow_empty = TRUE))
})

# ============================================================================
# validate_numeric_param tests
# ============================================================================

test_that("validate_numeric_param accepts valid numbers", {
  expect_true(epiRomics:::validate_numeric_param(5, "param1"))
  expect_true(epiRomics:::validate_numeric_param(0.5, "param1"))
  expect_true(epiRomics:::validate_numeric_param(-10, "param1"))
})

test_that("validate_numeric_param rejects non-numeric", {
  expect_error(epiRomics:::validate_numeric_param("5", "param1"),
    "must be numeric")
  expect_error(epiRomics:::validate_numeric_param(TRUE, "param1"),
    "must be numeric")
})

test_that("validate_numeric_param rejects vector of length > 1", {
  expect_error(epiRomics:::validate_numeric_param(c(1, 2), "param1"),
    "must be a single numeric value")
})

test_that("validate_numeric_param rejects NA", {
  expect_error(epiRomics:::validate_numeric_param(NA_real_, "param1"),
    "cannot be NA")
})

test_that("validate_numeric_param enforces min_val", {
  expect_error(epiRomics:::validate_numeric_param(0, "param1",
    min_val = 1), "must be >= 1")
  expect_true(epiRomics:::validate_numeric_param(1, "param1", min_val = 1))
})

test_that("validate_numeric_param enforces max_val", {
  expect_error(epiRomics:::validate_numeric_param(10, "param1",
    max_val = 5), "must be <= 5")
  expect_true(epiRomics:::validate_numeric_param(5, "param1", max_val = 5))
})

test_that("validate_numeric_param enforces both min and max", {
  expect_true(epiRomics:::validate_numeric_param(5, "param1",
    min_val = 1, max_val = 10))
  expect_error(epiRomics:::validate_numeric_param(0, "param1",
    min_val = 1, max_val = 10), "must be >= 1")
  expect_error(epiRomics:::validate_numeric_param(11, "param1",
    min_val = 1, max_val = 10), "must be <= 10")
})

# ============================================================================
# validate_logical_param tests
# ============================================================================

test_that("validate_logical_param accepts TRUE and FALSE", {
  expect_true(epiRomics:::validate_logical_param(TRUE, "param1"))
  expect_true(epiRomics:::validate_logical_param(FALSE, "param1"))
})

test_that("validate_logical_param rejects non-logical", {
  expect_error(epiRomics:::validate_logical_param(1, "param1"),
    "must be logical")
  expect_error(epiRomics:::validate_logical_param("TRUE", "param1"),
    "must be logical")
})

test_that("validate_logical_param rejects vector of length > 1", {
  expect_error(epiRomics:::validate_logical_param(c(TRUE, FALSE), "param1"),
    "must be a single logical value")
})

test_that("validate_logical_param rejects NA", {
  expect_error(epiRomics:::validate_logical_param(NA, "param1"),
    "cannot be NA")
})

# ============================================================================
# .chromatin_state_palette tests
# ============================================================================

test_that("chromatin_state_palette returns named character vector", {
  pal <- epiRomics:::.chromatin_state_palette()
  expect_type(pal, "character")
  expect_true(length(pal) > 0)
  expect_true(!is.null(names(pal)))
})

test_that("chromatin_state_palette contains all 6 required states", {
  pal <- epiRomics:::.chromatin_state_palette()
  required <- c("active", "bivalent", "poised", "primed", "repressed", "unmarked")
  expect_true(all(required %in% names(pal)))
})

test_that("chromatin_state_palette colors are valid hex codes", {
  pal <- epiRomics:::.chromatin_state_palette()
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", pal)))
})

test_that("chromatin_state_palette has no duplicate colors for distinct states", {
  pal <- epiRomics:::.chromatin_state_palette()
  # Each of the 6 states should have a different color
  expect_true(pal["active"] != pal["repressed"])
  expect_true(pal["active"] != pal["bivalent"])
  expect_true(pal["poised"] != pal["primed"])
})

# ============================================================================
# .aggregate_signal_over_regions tests
# ============================================================================

test_that("aggregate_signal returns correct mean per region", {
  # BigWig-like data: 3 bins on chr1 with known scores
  bw <- GRanges("chr1", IRanges(c(100, 200, 300), c(199, 299, 399)),
    score = c(10, 20, 30))

  # One region overlapping bins 1 and 2
  regions <- GRanges("chr1", IRanges(100, 299))

  result <- epiRomics:::.aggregate_signal_over_regions(bw, regions)
  expect_equal(result, mean(c(10, 20)))
})

test_that("aggregate_signal returns 0 for non-overlapping regions", {
  bw <- GRanges("chr1", IRanges(100, 200), score = 10)
  regions <- GRanges("chr2", IRanges(100, 200))  # different chromosome

  result <- epiRomics:::.aggregate_signal_over_regions(bw, regions)
  expect_equal(result, 0)
})

test_that("aggregate_signal handles empty bw_data", {
  bw <- GRanges()
  regions <- GRanges("chr1", IRanges(100, 200))

  result <- epiRomics:::.aggregate_signal_over_regions(bw, regions)
  expect_equal(result, 0)
})

test_that("aggregate_signal handles empty regions", {
  bw <- GRanges("chr1", IRanges(100, 200), score = 10)
  regions <- GRanges()

  result <- epiRomics:::.aggregate_signal_over_regions(bw, regions)
  expect_equal(length(result), 0)
})

test_that("aggregate_signal handles multiple regions correctly", {
  bw <- GRanges("chr1",
    IRanges(c(100, 200, 500, 600), c(199, 299, 599, 699)),
    score = c(10, 20, 50, 60))

  regions <- GRanges("chr1", IRanges(c(100, 500), c(299, 699)))

  result <- epiRomics:::.aggregate_signal_over_regions(bw, regions)
  expect_equal(length(result), 2)
  expect_equal(result[1], mean(c(10, 20)))  # region 1 overlaps bins 1,2
  expect_equal(result[2], mean(c(50, 60)))  # region 2 overlaps bins 3,4
})

test_that("aggregate_signal respects custom agg_fun", {
  bw <- GRanges("chr1",
    IRanges(c(100, 200, 300), c(199, 299, 399)),
    score = c(10, 20, 30))

  regions <- GRanges("chr1", IRanges(100, 399))

  result_max <- epiRomics:::.aggregate_signal_over_regions(bw, regions,
    agg_fun = max)
  expect_equal(result_max, 30)

  result_sum <- epiRomics:::.aggregate_signal_over_regions(bw, regions,
    agg_fun = sum)
  expect_equal(result_sum, 60)
})

# ============================================================================
# .aggregate_weighted_signal tests
# ============================================================================

test_that("aggregate_weighted_signal computes score*width sums", {
  # 2 bins: score=10,width=100; score=20,width=100
  bw <- GRanges("chr1",
    IRanges(c(100, 200), c(199, 299)),
    score = c(10, 20))

  regions <- GRanges("chr1", IRanges(100, 299))

  result <- epiRomics:::.aggregate_weighted_signal(bw, regions)
  # 10*100 + 20*100 = 3000
  expect_equal(result, 10 * 100 + 20 * 100)
})

test_that("aggregate_weighted_signal returns 0 for non-overlapping", {
  bw <- GRanges("chr1", IRanges(100, 200), score = 10)
  regions <- GRanges("chr2", IRanges(100, 200))

  result <- epiRomics:::.aggregate_weighted_signal(bw, regions)
  expect_equal(result, 0)
})

test_that("aggregate_weighted_signal handles empty inputs", {
  expect_equal(
    epiRomics:::.aggregate_weighted_signal(GRanges(), GRanges("chr1", IRanges(1, 100))),
    0)
  expect_equal(
    length(epiRomics:::.aggregate_weighted_signal(GRanges("chr1", IRanges(1, 100), score = 1), GRanges())),
    0)
})

test_that("aggregate_weighted_signal handles multiple regions", {
  bw <- GRanges("chr1",
    IRanges(c(100, 500), c(199, 599)),
    score = c(5, 10))

  regions <- GRanges("chr1", IRanges(c(100, 500), c(199, 599)))

  result <- epiRomics:::.aggregate_weighted_signal(bw, regions)
  expect_equal(result[1], 5 * 100)    # score=5, width=100
  expect_equal(result[2], 10 * 100)   # score=10, width=100
})

# ============================================================================
# resolve_txdb tests
# ============================================================================

test_that("resolve_txdb errors on invalid format", {
  expect_error(epiRomics:::resolve_txdb("no_double_colon"),
    "Invalid TxDb format")
  expect_error(epiRomics:::resolve_txdb(""),
    "Invalid TxDb format")
})

test_that("resolve_txdb correctly parses package::object format", {
  # We can't easily test the actual resolution without the package loaded,

  # but we can verify it doesn't error on a valid format
  # TxDb.Hsapiens.UCSC.hg38.knownGene should be available in test env
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  txdb <- epiRomics:::resolve_txdb(
    "TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene")
  expect_s4_class(txdb, "TxDb")
})

# ============================================================================
# .fast_bw_signal / .fast_bw_replicate_mean edge case tests
# ============================================================================

test_that(".fast_bw_signal returns zeros for empty regions", {
  empty_gr <- GRanges()
  result <- epiRomics:::.fast_bw_signal("nonexistent.bw", empty_gr, type = "mean")
  expect_equal(length(result), 0L)
  expect_is(result, "numeric")
})

test_that(".fast_bw_replicate_mean handles empty inputs", {
  empty_gr <- GRanges()
  result <- epiRomics:::.fast_bw_replicate_mean(
    character(0), empty_gr, type = "mean")
  expect_equal(length(result), 0L)
  result2 <- epiRomics:::.fast_bw_replicate_mean(
    "nonexistent.bw", empty_gr, type = "mean")
  expect_equal(length(result2), 0L)
})

test_that(".fast_bw_replicate_mean with single file matches .fast_bw_signal", {
  # Can only test with actual BigWig file
  bw_file <- system.file("extdata", "BigWigs",
    "Kaestner_Alpha.atac.signal.subset.multi_chr.bigwig", package = "epiRomics")
  skip_if(bw_file == "", "No BigWig test file available")
  regions <- GRanges("chr1", IRanges(c(1000000, 2000000), c(1010000, 2010000)))
  single <- epiRomics:::.fast_bw_signal(bw_file, regions, type = "mean")
  replicate <- epiRomics:::.fast_bw_replicate_mean(bw_file, regions, type = "mean")
  expect_equal(single, replicate, tolerance = 1e-6)
})

test_that(".fast_bw_weighted_signal returns zeros for empty regions", {
  empty_gr <- GRanges()
  result <- epiRomics:::.fast_bw_weighted_signal("nonexistent.bw", empty_gr)
  expect_equal(length(result), 0L)
})

test_that(".aggregate_signal_over_regions handles empty GRanges", {
  empty_bw <- GRanges()
  regions <- GRanges("chr1", IRanges(1000, 2000))
  result <- epiRomics:::.aggregate_signal_over_regions(empty_bw, regions)
  expect_equal(result, 0)

  # Both empty
  result2 <- epiRomics:::.aggregate_signal_over_regions(empty_bw, GRanges())
  expect_equal(length(result2), 0L)
})

# ============================================================
# .detect_cores tests
# ============================================================

test_that(".detect_cores returns positive integer", {
  n <- epiRomics:::.detect_cores()
  expect_true(is.integer(n))
  expect_true(n >= 1L)
})

test_that(".detect_cores respects max_cores cap", {
  n <- epiRomics:::.detect_cores(max_cores = 2L)
  expect_true(n >= 1L)
  expect_true(n <= 2L)
})

test_that(".detect_cores never returns NA", {
  n <- epiRomics:::.detect_cores()
  expect_false(is.na(n))
})

test_that(".detect_cores with max_cores=1 returns 1", {
  n <- epiRomics:::.detect_cores(max_cores = 1L)
  expect_equal(n, 1L)
})
