# Synthetic tests for get_regions_of_interest() expanded modes
# No database rebuild — uses in-memory GRanges + mock epiRomicsS4

# --- Helpers ---
make_mock_dB <- function(n = 10, with_symbol = TRUE) {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = base::seq(1000L, by = 5000L, length.out = n),
      width = 500L
    )
  )
  if (with_symbol) {
    gr$SYMBOL <- base::paste0("GENE", base::seq_len(n))
  }
  methods::new("epiRomicsS4",
    annotations = gr,
    meta = base::data.frame(dummy = "test", stringsAsFactors = FALSE),
    txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    organism = "Homo sapiens",
    genome = "hg38"
  )
}

make_test_gr <- function(start = 1000L, width = 500L, n = 1) {
  GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
      start = base::seq(start, by = 5000L, length.out = n),
      width = width
    )
  )
}

make_temp_bed <- function(regions_gr) {
  tf <- base::tempfile(fileext = ".bed")
  rtracklayer::export(regions_gr, tf, format = "BED")
  tf
}

# ============================================================
# Validation tests
# ============================================================
test_that("rejects non-epiRomicsS4 input", {
  expect_error(
    get_regions_of_interest("not_a_dB", make_test_gr()),
    "epiRomicsS4"
  )
})

test_that("rejects empty annotations", {
  dB <- make_mock_dB(10)
  dB@annotations <- GenomicRanges::GRanges()
  expect_error(
    get_regions_of_interest(dB, make_test_gr()),
    "empty|no annotations"
  )
})

test_that("rejects non-GRanges test_regions", {
  dB <- make_mock_dB(10)
  expect_error(
    get_regions_of_interest(dB, "not_granges"),
    "GRanges"
  )
})

test_that("granges mode requires test_regions", {
  dB <- make_mock_dB(10)
  expect_error(
    get_regions_of_interest(dB, input_type = "granges"),
    "test_regions.*required"
  )
})

test_that("bed mode requires bed_path", {
  dB <- make_mock_dB(10)
  expect_error(
    get_regions_of_interest(dB, input_type = "bed"),
    "bed_path required"
  )
})

test_that("genelist mode requires gene_list", {
  dB <- make_mock_dB(10)
  expect_error(
    get_regions_of_interest(dB, input_type = "genelist"),
    "gene_list required"
  )
})

test_that("genelist mode rejects empty gene_list", {
  dB <- make_mock_dB(10)
  expect_error(
    get_regions_of_interest(dB, input_type = "genelist",
      gene_list = character(0)),
    "gene_list required"
  )
})

test_that("invalid input_type rejected", {
  dB <- make_mock_dB(10)
  expect_error(
    get_regions_of_interest(dB, input_type = "invalid"),
    "arg"
  )
})

# ============================================================
# GRanges mode (backward-compatible)
# ============================================================
test_that("GRanges overlap returns correct subset", {
  dB <- make_mock_dB(10)
  # First region starts at 1000, second at 6000, third at 11000...
  test_gr <- make_test_gr(start = 1000L, width = 600L, n = 1)
  result <- get_regions_of_interest(dB, test_gr)
  expect_s4_class(result, "epiRomicsS4")
  expect_equal(length(result@annotations), 1)
  expect_equal(GenomicRanges::start(result@annotations), 1000L)
})

test_that("GRanges no overlap warns", {
  dB <- make_mock_dB(10)
  test_gr <- make_test_gr(start = 999000L, width = 100L, n = 1)
  expect_warning(
    get_regions_of_interest(dB, test_gr),
    "No overlapping"
  )
})

test_that("GRanges multiple overlaps work", {
  dB <- make_mock_dB(10)
  # Cover first two regions: 1000-1500 and 6000-6500
  test_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(900L, 5900L), width = 700L)
  )
  result <- get_regions_of_interest(dB, test_gr)
  expect_equal(length(result@annotations), 2)
})

test_that("GRanges mode preserves metadata columns", {
  dB <- make_mock_dB(10)
  test_gr <- make_test_gr(start = 1000L, width = 600L, n = 1)
  result <- get_regions_of_interest(dB, test_gr)
  expect_true("SYMBOL" %in% names(S4Vectors::mcols(result@annotations)))
  expect_equal(result@annotations$SYMBOL, "GENE1")
})

test_that("GRanges mode ignores input_type when test_regions provided", {
  dB <- make_mock_dB(10)
  test_gr <- make_test_gr(start = 1000L, width = 600L, n = 1)
  # Even with input_type = "bed", should use test_regions directly
  result <- get_regions_of_interest(dB, test_gr, input_type = "bed")
  expect_equal(length(result@annotations), 1)
})

# ============================================================
# BED mode
# ============================================================
test_that("BED mode imports and filters correctly", {
  dB <- make_mock_dB(10)
  bed_gr <- make_test_gr(start = 1000L, width = 600L, n = 1)
  bed_file <- make_temp_bed(bed_gr)
  on.exit(base::unlink(bed_file))

  result <- get_regions_of_interest(dB, input_type = "bed",
    bed_path = bed_file)
  expect_equal(length(result@annotations), 1)
})

test_that("BED mode missing file errors", {
  dB <- make_mock_dB(10)
  expect_error(
    get_regions_of_interest(dB, input_type = "bed",
      bed_path = "/nonexistent/file.bed"),
    "not found"
  )
})

test_that("BED mode no overlap warns", {
  dB <- make_mock_dB(10)
  bed_gr <- make_test_gr(start = 999000L, width = 100L, n = 1)
  bed_file <- make_temp_bed(bed_gr)
  on.exit(base::unlink(bed_file))

  expect_warning(
    get_regions_of_interest(dB, input_type = "bed",
      bed_path = bed_file),
    "No overlapping"
  )
})

# ============================================================
# Genelist mode
# ============================================================
test_that("genelist filters by SYMBOL", {
  dB <- make_mock_dB(10)
  result <- get_regions_of_interest(dB, input_type = "genelist",
    gene_list = c("GENE1", "GENE3", "GENE5"))
  expect_equal(length(result@annotations), 3)
  expect_true(all(result@annotations$SYMBOL %in% c("GENE1", "GENE3", "GENE5")))
})

test_that("genelist no matches warns", {
  dB <- make_mock_dB(10)
  expect_warning(
    get_regions_of_interest(dB, input_type = "genelist",
      gene_list = c("NONEXISTENT1", "NONEXISTENT2")),
    "No annotations matched"
  )
})

test_that("genelist all match returns all", {
  dB <- make_mock_dB(5)
  genes <- paste0("GENE", 1:5)
  result <- get_regions_of_interest(dB, input_type = "genelist",
    gene_list = genes)
  expect_equal(length(result@annotations), 5)
})

test_that("genelist requires SYMBOL column", {
  dB <- make_mock_dB(10, with_symbol = FALSE)
  expect_error(
    get_regions_of_interest(dB, input_type = "genelist",
      gene_list = c("GENE1")),
    "SYMBOL column"
  )
})

# ============================================================
# Combined mode
# ============================================================
test_that("combined mode unions BED + genelist evidence", {
  dB <- make_mock_dB(10)
  # BED covers region 1 (start=1000)
  bed_gr <- make_test_gr(start = 1000L, width = 600L, n = 1)
  bed_file <- make_temp_bed(bed_gr)
  on.exit(base::unlink(bed_file))

  # Genelist covers GENE5 (different region)
  result <- get_regions_of_interest(dB, input_type = "combined",
    bed_path = bed_file, gene_list = c("GENE5"))
  # Should retain GENE1 (BED overlap) + GENE5 (genelist match) = 2
  expect_equal(length(result@annotations), 2)
  expect_true(all(c("GENE1", "GENE5") %in% result@annotations$SYMBOL))
})

test_that("combined mode with only BED works", {
  dB <- make_mock_dB(10)
  bed_gr <- make_test_gr(start = 1000L, width = 600L, n = 1)
  bed_file <- make_temp_bed(bed_gr)
  on.exit(base::unlink(bed_file))

  result <- get_regions_of_interest(dB, input_type = "combined",
    bed_path = bed_file)
  expect_equal(length(result@annotations), 1)
})

test_that("combined mode with only genelist works", {
  dB <- make_mock_dB(10)
  result <- get_regions_of_interest(dB, input_type = "combined",
    gene_list = c("GENE2", "GENE4"))
  expect_equal(length(result@annotations), 2)
})

test_that("combined mode no evidence warns", {
  dB <- make_mock_dB(10)
  expect_warning(
    get_regions_of_interest(dB, input_type = "combined"),
    "No overlapping"
  )
})

test_that("combined mode deduplicates overlapping evidence", {
  dB <- make_mock_dB(10)
  # BED covers region 1, genelist also includes GENE1
  bed_gr <- make_test_gr(start = 1000L, width = 600L, n = 1)
  bed_file <- make_temp_bed(bed_gr)
  on.exit(base::unlink(bed_file))

  result <- get_regions_of_interest(dB, input_type = "combined",
    bed_path = bed_file, gene_list = c("GENE1"))
  # GENE1 hit by both sources, should still be 1 region not 2
  expect_equal(length(result@annotations), 1)
})

# ============================================================
# Return structure
# ============================================================
test_that("all modes return epiRomicsS4 with preserved slots", {
  dB <- make_mock_dB(10)
  test_gr <- make_test_gr(start = 1000L, width = 600L, n = 1)
  result <- get_regions_of_interest(dB, test_gr)
  expect_s4_class(result, "epiRomicsS4")
  expect_equal(result@txdb, "TxDb.Hsapiens.UCSC.hg38.knownGene")
  expect_equal(result@organism, "Homo sapiens")
  expect_equal(result@genome, "hg38")
})
