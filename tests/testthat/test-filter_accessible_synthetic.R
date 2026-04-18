## Synthetic tests for filter_accessible_regions() — 4 modes
## Tests validation, genelist mode, bed mode overlap, scope, and combined logic.
## Signal mode requires real BigWig files and is tested in integration tests.
##
## NOTE on suppressMessages() usage throughout this file:
## filter_accessible_regions() emits info-level progress messages
## ("Matched N genes", "Retained N regions", "Loading BED"). These are
## intentionally noisy for interactive use but are NOT under test here —
## every test_that() block asserts a structural property of the returned
## data.frame (e.g. result$atac_accessible logical vector). The suppression
## keeps test output clean without hiding any behavior the tests pin.
## suppressWarnings() — where used — is individually justified inline.

# --- Helper: build a minimal putative enhancers data.frame ---
.make_pe <- function(n = 5, with_symbol = TRUE, with_tss = FALSE) {
  pe <- data.frame(
    chr = rep("chr1", n),
    start = seq(1000, by = 10000, length.out = n),
    end = seq(2000, by = 10000, length.out = n),
    stringsAsFactors = FALSE
  )
  if (with_symbol) {
    pe$SYMBOL <- c("INS", "GCG", "SST", "PPY", "MAFA")[seq_len(n)]
  }
  if (with_tss) {
    pe$distanceToTSS <- c(100, 500, 5000, 50000, 100000)[seq_len(n)]
  }
  pe
}

# ============================================================================
# VALIDATION TESTS
# ============================================================================

test_that("filter_accessible rejects empty data.frame", {
  expect_error(
    filter_accessible_regions(data.frame(), mode = "genelist",
      gene_list = "INS"),
    "non-empty data.frame"
  )
})

test_that("filter_accessible rejects invalid mode", {
  pe <- .make_pe()
  expect_error(
    filter_accessible_regions(pe, mode = "invalid"),
    "'arg' should be one of"
  )
})

test_that("filter_accessible rejects invalid scope", {
  pe <- .make_pe()
  expect_error(
    filter_accessible_regions(pe, mode = "genelist", scope = "bogus",
      gene_list = "INS"),
    "'arg' should be one of"
  )
})

test_that("filter_accessible rejects missing chr/start/end", {
  pe <- data.frame(gene = "INS", start = 1, end = 100,
    stringsAsFactors = FALSE)
  expect_error(
    filter_accessible_regions(pe, mode = "genelist", gene_list = "INS"),
    "missing required columns.*chr"
  )
})

test_that("signal mode requires track_connection", {
  pe <- .make_pe()
  expect_error(
    filter_accessible_regions(pe, mode = "signal"),
    "track_connection required"
  )
})

test_that("bed mode requires bed_path", {
  pe <- .make_pe()
  expect_error(
    filter_accessible_regions(pe, mode = "bed"),
    "bed_path required"
  )
})

test_that("genelist mode requires gene_list", {
  pe <- .make_pe()
  expect_error(
    filter_accessible_regions(pe, mode = "genelist"),
    "gene_list required"
  )
})

test_that("genelist mode requires SYMBOL column", {
  pe <- .make_pe(with_symbol = FALSE)
  expect_error(
    filter_accessible_regions(pe, mode = "genelist",
      gene_list = c("INS")),
    "SYMBOL"
  )
})

# ============================================================================
# GENELIST MODE TESTS
# ============================================================================

test_that("genelist mode flags matching genes", {
  pe <- .make_pe()
  result <- suppressMessages(
    filter_accessible_regions(pe, mode = "genelist",
      gene_list = c("INS", "GCG"))
  )
  expect_true("atac_accessible" %in% names(result))
  expect_equal(result$atac_accessible, c(TRUE, TRUE, FALSE, FALSE, FALSE))
})

test_that("genelist mode with no matches flags all FALSE", {
  pe <- .make_pe()
  result <- suppressMessages(
    filter_accessible_regions(pe, mode = "genelist",
      gene_list = c("BRCA1", "TP53"))
  )
  expect_true(all(!result$atac_accessible))
})

test_that("genelist mode with all matches flags all TRUE", {
  pe <- .make_pe()
  result <- suppressMessages(
    filter_accessible_regions(pe, mode = "genelist",
      gene_list = c("INS", "GCG", "SST", "PPY", "MAFA"))
  )
  expect_true(all(result$atac_accessible))
})

# ============================================================================
# SCOPE TESTS (filter_distal vs filter_all)
# ============================================================================

test_that("filter_distal retains promoter-proximal regions", {
  pe <- .make_pe(with_tss = TRUE)
  # Regions 1,2 are within 2000bp of TSS (100, 500)
  # With genelist only matching "SST" (region 3, distal)
  result <- suppressMessages(
    filter_accessible_regions(pe, mode = "genelist",
      scope = "filter_distal", gene_list = c("SST"))
  )
  # Regions 1,2 retained (promoter), region 3 retained (gene match),
  # regions 4,5 filtered out
  expect_equal(result$atac_accessible, c(TRUE, TRUE, TRUE, FALSE, FALSE))
})

test_that("filter_all does not retain promoter-proximal", {
  pe <- .make_pe(with_tss = TRUE)
  result <- suppressMessages(
    filter_accessible_regions(pe, mode = "genelist",
      scope = "filter_all", gene_list = c("SST"))
  )
  # Only region 3 matches
  expect_equal(result$atac_accessible, c(FALSE, FALSE, TRUE, FALSE, FALSE))
})

test_that("filter_distal without distanceToTSS column marks none as promoter", {
  pe <- .make_pe(with_tss = FALSE)
  result <- suppressMessages(
    filter_accessible_regions(pe, mode = "genelist",
      scope = "filter_distal", gene_list = c("INS"))
  )
  # No distanceToTSS -> no promoter protection -> only INS passes
  expect_equal(result$atac_accessible, c(TRUE, FALSE, FALSE, FALSE, FALSE))
})

# ============================================================================
# BED MODE TESTS (using temporary BED files)
# ============================================================================

test_that("bed mode rejects missing file", {
  pe <- .make_pe()
  expect_error(
    filter_accessible_regions(pe, mode = "bed",
      bed_path = "/nonexistent/file.bed"),
    "BED file not found"
  )
})

test_that("bed mode finds overlapping regions", {
  pe <- .make_pe()
  # Create a temp BED file overlapping regions 1 and 3
  bed_file <- tempfile(fileext = ".bed")
  writeLines(c(
    "chr1\t1000\t2000\tpeak1",
    "chr1\t21000\t22000\tpeak2"
  ), bed_file)
  on.exit(unlink(bed_file))

  result <- suppressMessages(
    filter_accessible_regions(pe, mode = "bed", bed_path = bed_file)
  )
  expect_true(result$atac_accessible[1])   # overlaps peak1
  expect_true(result$atac_accessible[3])   # overlaps peak2
  expect_false(result$atac_accessible[5])  # no overlap
})

test_that("bed mode with non-overlapping BED gives all FALSE", {
  pe <- .make_pe()
  bed_file <- tempfile(fileext = ".bed")
  writeLines("chr2\t1000\t2000\tpeak1", bed_file)
  on.exit(unlink(bed_file))

  # Expected 'no sequence levels in common' is the behavior under test
  result <- suppressWarnings(suppressMessages(
    filter_accessible_regions(pe, mode = "bed", bed_path = bed_file)
  ))
  expect_true(all(!result$atac_accessible))
})

# ============================================================================
# COMBINED MODE TESTS
# ============================================================================

test_that("combined mode unions genelist and bed evidence", {
  pe <- .make_pe()
  bed_file <- tempfile(fileext = ".bed")
  writeLines("chr1\t1000\t2000\tpeak1", bed_file)  # overlaps region 1
  on.exit(unlink(bed_file))

  result <- suppressMessages(
    filter_accessible_regions(pe, mode = "combined",
      bed_path = bed_file, gene_list = c("SST"))  # region 3
  )
  # Region 1 from bed, region 3 from genelist
  expect_equal(result$atac_accessible, c(TRUE, FALSE, TRUE, FALSE, FALSE))
})

test_that("combined mode with no evidence gives all FALSE", {
  pe <- .make_pe()
  bed_file <- tempfile(fileext = ".bed")
  writeLines("chr2\t1000\t2000\tpeak1", bed_file)
  on.exit(unlink(bed_file))

  # Expected 'no sequence levels in common' is the behavior under test
  result <- suppressWarnings(suppressMessages(
    filter_accessible_regions(pe, mode = "combined",
      bed_path = bed_file, gene_list = c("BRCA1"))
  ))
  expect_true(all(!result$atac_accessible))
})

test_that("combined mode with scope=filter_distal retains promoters", {
  pe <- .make_pe(with_tss = TRUE)
  result <- suppressMessages(
    filter_accessible_regions(pe, mode = "combined",
      scope = "filter_distal", gene_list = c("MAFA"))
  )
  # Regions 1,2 from promoter proximity, region 5 from genelist
  expect_equal(result$atac_accessible, c(TRUE, TRUE, FALSE, FALSE, TRUE))
})

# ============================================================================
# RETURN STRUCTURE
# ============================================================================

test_that("result preserves original columns plus atac_accessible", {
  pe <- .make_pe()
  pe$custom_col <- 1:5
  result <- suppressMessages(
    filter_accessible_regions(pe, mode = "genelist",
      gene_list = c("INS"))
  )
  expect_true("custom_col" %in% names(result))
  expect_true("atac_accessible" %in% names(result))
  expect_equal(nrow(result), 5)
})

# ============================================================================
# SIGNAL CONDITION TESTS
# ============================================================================

test_that("signal mode emits proper warning on BigWig import failure", {
  pe <- .make_pe()
  tc <- data.frame(
    path = "/tmp/nonexistent_test.bw",
    name = "test_atac",
    color = "#000000",
    type = "atac",
    stringsAsFactors = FALSE
  )
  # Capture all warnings to check that our code emits a proper warning()
  # (rtracklayer may also emit its own warnings about the missing file)
  warnings_caught <- character(0)
  # suppressMessages: progress messages (see file-level note). Warnings are
  # NOT suppressed — they are captured via withCallingHandlers below because
  # the "failed to import" warning IS under test.
  result <- suppressMessages(
    withCallingHandlers(
      filter_accessible_regions(pe, mode = "signal",
        track_connection = tc),
      warning = function(w) {
        warnings_caught <<- c(warnings_caught, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  )
  # Our code must emit a warning matching "failed to import"
  import_warnings <- grep("failed to import", warnings_caught, value = TRUE)
  expect_true(length(import_warnings) >= 1,
    info = paste("Expected 'failed to import' warning, got:",
      paste(warnings_caught, collapse = "; ")))
  # The warning must include the file path
  expect_true(any(grepl("/tmp/nonexistent_test\\.bw", import_warnings)),
    info = "Warning should include the BigWig file path")
  # Function returns a data.frame (graceful degradation)
  expect_true(is.data.frame(result))
})

test_that("signal mode returns zero signal for failed track", {
  pe <- .make_pe()
  tc <- data.frame(
    path = "/tmp/nonexistent_test.bw",
    name = "test_atac",
    color = "#000000",
    type = "atac",
    stringsAsFactors = FALSE
  )
  # Expected "failed to import" warning (tested elsewhere — see the
  # preceding test_that block) plus info-level progress messages. Suppressing
  # both here because this test asserts the graceful-degradation contract:
  # the function returns a data.frame with the new column containing zeros.
  result <- suppressMessages(suppressWarnings(
    filter_accessible_regions(pe, mode = "signal",
      track_connection = tc)
  ))
  expect_true(is.data.frame(result))
  expect_true("atac_test_atac" %in% names(result))
  expect_equal(result$atac_test_atac, rep(0, 5))
})
