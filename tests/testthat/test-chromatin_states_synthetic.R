# ============================================================================
# SYNTHETIC DATA TESTS for classify_chromatin_states
# These tests use small, hand-crafted GRanges to verify classification rules.
# No database rebuild, no file I/O — pure logic verification in seconds.
#
# 6-state scheme: active, poised, repressed, bivalent, primed, unmarked
# Genomic context (promoter/gene_body/intergenic) is a separate column.
# ============================================================================

library(testthat)
library(epiRomics)
library(GenomicRanges)

# Helper: create a minimal mock epiRomicsS4 with synthetic histone marks
make_mock_dB <- function(marks_list, genome = "hg38") {

  # marks_list: named list of GRanges, e.g. list(h3k27ac = GRanges(...), ...)
  all_gr <- GRanges()
  meta_rows <- list()

  for (mark_name in names(marks_list)) {
    gr <- marks_list[[mark_name]]
    gr$type <- paste0(genome, "_custom_", mark_name)
    all_gr <- c(all_gr, gr)
    meta_rows[[mark_name]] <- data.frame(
      name = mark_name, type = "histone",
      stringsAsFactors = FALSE)
  }

  meta <- do.call(rbind, meta_rows)

  methods::new("epiRomicsS4",
    annotations = all_gr,
    meta = meta,
    genome = genome,
    txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene",
    organism = "org.Hs.eg.db")
}

# ============================================================================
# Test 1: active = H3K4me1 + H3K27ac (no H3K27me3)
# ============================================================================
test_that("active classified from H3K4me1 + H3K27ac", {
  marks <- list(
    h3k4me1 = GRanges("chr1", IRanges(1000, 2000)),
    h3k27ac = GRanges("chr1", IRanges(1000, 2000))
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = FALSE)
  expect_true(nrow(result) > 0)
  expect_true("active" %in% result$chromatin_state)
})

# ============================================================================
# Test 2: active = H3K4me3 + H3K27ac (no H3K27me3)
# ============================================================================
test_that("active classified from H3K4me3 + H3K27ac", {
  marks <- list(
    h3k4me3 = GRanges("chr1", IRanges(5000, 6000)),
    h3k27ac = GRanges("chr1", IRanges(5000, 6000))
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = FALSE)
  expect_true(nrow(result) > 0)
  expect_true("active" %in% result$chromatin_state)
})

# ============================================================================
# Test 3: bivalent = H3K4me3 + H3K27me3
# ============================================================================
test_that("bivalent classified from H3K4me3 + H3K27me3", {
  marks <- list(
    h3k4me3 = GRanges("chr1", IRanges(5000, 6000)),
    h3k27me3 = GRanges("chr1", IRanges(5000, 6000))
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = FALSE)
  expect_true(nrow(result) > 0)
  expect_true("bivalent" %in% result$chromatin_state)
})

# ============================================================================
# Test 4: poised = H3K4me1 + H3K27me3
# ============================================================================
test_that("poised classified from H3K4me1 + H3K27me3", {
  marks <- list(
    h3k4me1 = GRanges("chr1", IRanges(10000, 11000)),
    h3k27me3 = GRanges("chr1", IRanges(10000, 11000))
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = FALSE)
  expect_true("poised" %in% result$chromatin_state)
})

# ============================================================================
# Test 5: repressed = H3K27me3 + H3K9me3
# ============================================================================
test_that("repressed classified from H3K27me3 + H3K9me3", {
  marks <- list(
    h3k27me3 = GRanges("chr1", IRanges(20000, 21000)),
    h3k9me3 = GRanges("chr1", IRanges(20000, 21000))
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = FALSE)
  expect_true("repressed" %in% result$chromatin_state)
})

# ============================================================================
# Test 6: repressed = H3K27me3 alone (Polycomb)
# ============================================================================
test_that("repressed classified from H3K27me3 alone", {
  marks <- list(
    h3k27me3 = GRanges("chr1", IRanges(30000, 31000))
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = FALSE)
  expect_true("repressed" %in% result$chromatin_state)
})

# ============================================================================
# Test 7: repressed = H3K9me3 alone (heterochromatin collapsed into repressed)
# ============================================================================
test_that("repressed classified from H3K9me3 alone (heterochromatin)", {
  marks <- list(
    h3k9me3 = GRanges("chr1", IRanges(40000, 41000))
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = FALSE)
  expect_true("repressed" %in% result$chromatin_state)
})

# ============================================================================
# Test 8: unmarked = H3K36me3 alone (transcribed gene body, not active regulatory)
# ============================================================================
test_that("unmarked classified from H3K36me3 alone (gene body)", {
  marks <- list(
    h3k36me3 = GRanges("chr1", IRanges(50000, 51000))
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = FALSE)
  expect_true("unmarked" %in% result$chromatin_state)
})

# ============================================================================
# Test 9: primed = H3K4me1 alone (no H3K27ac, no H3K27me3)
# ============================================================================
test_that("primed classified from H3K4me1 alone", {
  marks <- list(
    h3k4me1 = GRanges("chr1", IRanges(60000, 61000))
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = FALSE)
  expect_true("primed" %in% result$chromatin_state)
})

# ============================================================================
# Test 10: active = H3K4me1 + H3K27ac + H3K36me3 (strong enhancer → active)
# ============================================================================
test_that("active classified from H3K4me1 + H3K27ac + H3K36me3", {
  marks <- list(
    h3k4me1 = GRanges("chr1", IRanges(70000, 71000)),
    h3k27ac = GRanges("chr1", IRanges(70000, 71000)),
    h3k36me3 = GRanges("chr1", IRanges(70000, 71000))
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = FALSE)
  expect_true("active" %in% result$chromatin_state)
})

# ============================================================================
# Test 11: Priority — H3K27me3 + H3K9me3 = repressed
# ============================================================================
test_that("H3K27me3 + H3K9me3 is repressed (highest repressive priority)", {
  marks <- list(
    h3k27me3 = GRanges("chr1", IRanges(80000, 81000)),
    h3k9me3 = GRanges("chr1", IRanges(80000, 81000))
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = FALSE)
  region_state <- result$chromatin_state[1]
  expect_equal(region_state, "repressed")
})

# ============================================================================
# Test 12: Priority — H3K4me3 + H3K27me3 + H3K27ac = bivalent (H3K27me3 wins)
# ============================================================================
test_that("H3K4me3 + H3K27me3 + H3K27ac = bivalent (repressive overrides active)", {
  marks <- list(
    h3k4me3 = GRanges("chr1", IRanges(90000, 91000)),
    h3k27me3 = GRanges("chr1", IRanges(90000, 91000)),
    h3k27ac = GRanges("chr1", IRanges(90000, 91000))
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = FALSE)
  expect_true("bivalent" %in% result$chromatin_state)
})

# ============================================================================
# Test 13: genomic_context is populated when refine_by_tss = TRUE
# ============================================================================
test_that("genomic_context populated with TSS refinement enabled", {
  marks <- list(
    h3k4me3 = GRanges("chr1", IRanges(50000000, 50001000)),
    h3k27ac = GRanges("chr1", IRanges(50000000, 50001000))
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = TRUE)
  # State should always be 'active' regardless of TSS proximity
  expect_true("active" %in% result$chromatin_state)
  # genomic_context should be one of the valid values
  expect_true(all(result$genomic_context %in% c("promoter", "gene_body", "intergenic")))
})

# ============================================================================
# Test 14: All 6 valid states — no old state names ever appear
# ============================================================================
test_that("only 6 valid state labels are ever produced", {
  marks <- list(
    h3k4me3 = GRanges("chr1", IRanges(5000, 6000)),
    h3k27me3 = GRanges("chr1", IRanges(c(5000, 50000000), c(6000, 50001000)))
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = TRUE)
  valid_states <- c("active", "bivalent", "poised", "primed", "repressed", "unmarked")
  expect_true(all(result$chromatin_state %in% valid_states))
  # Explicitly check old names NEVER appear
  old_names <- c("active_promoter", "active_enhancer", "strong_enhancer",
    "active_gene_body", "poised_enhancer", "primed_enhancer",
    "bivalent_promoter", "bivalent_enhancer", "heterochromatin",
    "quiescent", "unclassified")
  for (old in old_names) {
    expect_false(old %in% result$chromatin_state,
      info = paste0("Old state '", old, "' should never appear"))
  }
})

# ============================================================================
# Test 15: Multiple non-overlapping regions classified independently
# ============================================================================
test_that("multiple non-overlapping regions get independent states", {
  marks <- list(
    h3k4me1 = GRanges("chr1", IRanges(c(1000, 50000), c(2000, 51000))),
    h3k27ac = GRanges("chr1", IRanges(1000, 2000)),           # only first
    h3k27me3 = GRanges("chr1", IRanges(50000, 51000))         # only second
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = FALSE)

  expect_true(nrow(result) >= 2)
  # First region: H3K4me1 + H3K27ac = active
  expect_true("active" %in% result$chromatin_state)
  # Second region: H3K4me1 + H3K27me3 = poised
  expect_true("poised" %in% result$chromatin_state)
})

# ============================================================================
# Test 16: Output structure is always correct
# ============================================================================
test_that("output has required columns regardless of input", {
  marks <- list(
    h3k4me1 = GRanges("chr1", IRanges(1000, 2000))
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = FALSE)

  required_cols <- c("seqnames", "start", "end", "width", "chromatin_state",
    "genomic_context", "marks_present", "n_marks", "is_hotspot")
  expect_true(all(required_cols %in% names(result)))
  expect_true(is.data.frame(result))
})

# ============================================================================
# Test 17: n_marks is consistent with marks_present
# ============================================================================
test_that("n_marks matches comma count in marks_present", {
  marks <- list(
    h3k4me1 = GRanges("chr1", IRanges(1000, 2000)),
    h3k27ac = GRanges("chr1", IRanges(1000, 2000)),
    h3k36me3 = GRanges("chr1", IRanges(1000, 2000))
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = FALSE)

  for (i in seq_len(nrow(result))) {
    if (result$marks_present[i] == "") {
      expect_equal(result$n_marks[i], 0L)
    } else {
      expected_n <- length(strsplit(result$marks_present[i], ",")[[1]])
      expect_equal(result$n_marks[i], expected_n,
        info = paste0("Row ", i, ": n_marks=", result$n_marks[i],
          " but marks_present='", result$marks_present[i], "'"))
    }
  }
})

# ============================================================================
# Test 18: Empty regions input returns empty data.frame
# ============================================================================
test_that("empty regions returns empty data.frame with warning", {
  marks <- list(
    h3k4me1 = GRanges("chr1", IRanges(1000, 2000))
  )
  db <- make_mock_dB(marks)
  empty_gr <- GRanges()
  expect_warning(
    result <- classify_chromatin_states(db, regions = empty_gr),
    "No regions to classify"
  )
  expect_equal(nrow(result), 0)
})

# ============================================================================
# Test 19: genomic_context values are always valid
# ============================================================================
test_that("genomic_context is always one of promoter/gene_body/intergenic", {
  marks <- list(
    h3k4me1 = GRanges("chr1", IRanges(c(1000, 50000, 100000), c(2000, 51000, 101000))),
    h3k27ac = GRanges("chr1", IRanges(1000, 2000))
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = TRUE)
  expect_true(all(result$genomic_context %in% c("promoter", "gene_body", "intergenic")))
})

# ============================================================================
# Test 20: is_hotspot = TRUE when 3+ marks
# ============================================================================
test_that("is_hotspot TRUE for regions with 3+ overlapping marks", {
  marks <- list(
    h3k4me1 = GRanges("chr1", IRanges(1000, 2000)),
    h3k27ac = GRanges("chr1", IRanges(1000, 2000)),
    h3k36me3 = GRanges("chr1", IRanges(1000, 2000))
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = FALSE)
  expect_true(any(result$is_hotspot))
  hotspots <- result[result$is_hotspot, ]
  expect_true(all(hotspots$n_marks >= 3))
})

# ============================================================================
# Test 21: is_hotspot = FALSE when < 3 marks
# ============================================================================
test_that("is_hotspot FALSE for regions with < 3 marks", {
  marks <- list(
    h3k4me1 = GRanges("chr1", IRanges(1000, 2000)),
    h3k27ac = GRanges("chr1", IRanges(1000, 2000))
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = FALSE)
  expect_false(any(result$is_hotspot))
})

# ============================================================================
# Test 22: H2A.Z alone = unmarked
# ============================================================================
test_that("H2A.Z alone classified as unmarked", {
  marks <- list(
    H2A.Z = GRanges("chr1", IRanges(1000, 2000))
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = FALSE)
  expect_true("unmarked" %in% result$chromatin_state)
})

# ============================================================================
# Test 23: H2A.Z + H3K27ac (no H3K4me3) = active
# ============================================================================
test_that("H2A.Z + H3K27ac classified as active", {
  marks <- list(
    H2A.Z = GRanges("chr1", IRanges(1000, 2000)),
    h3k27ac = GRanges("chr1", IRanges(1000, 2000))
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = FALSE)
  expect_true("active" %in% result$chromatin_state)
})

# ============================================================================
# Test 24: H2A.Z + H3K27me3 = poised
# ============================================================================
test_that("H2A.Z + H3K27me3 classified as poised", {
  marks <- list(
    H2A.Z = GRanges("chr1", IRanges(1000, 2000)),
    h3k27me3 = GRanges("chr1", IRanges(1000, 2000))
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = FALSE)
  expect_true("poised" %in% result$chromatin_state)
})

# ============================================================================
# Test 25: H3K4me1 + H3K27ac + H3K27me3 = poised (H3K27me3 takes precedence)
# ============================================================================
test_that("H3K4me1 + H3K27ac + H3K27me3 = poised (repressive wins)", {
  marks <- list(
    h3k4me1 = GRanges("chr1", IRanges(1000, 2000)),
    h3k27ac = GRanges("chr1", IRanges(1000, 2000)),
    h3k27me3 = GRanges("chr1", IRanges(1000, 2000))
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = FALSE)
  # H3K4me1 + H3K27me3 = poised (checked before H3K4me1 + H3K27ac)
  expect_true("poised" %in% result$chromatin_state)
})

# ============================================================================
# Test 26: No histone marks found → error
# ============================================================================
test_that("error when no histone marks in database", {
  db <- methods::new("epiRomicsS4",
    annotations = GRanges(),
    meta = data.frame(name = character(0), type = character(0),
      stringsAsFactors = FALSE),
    genome = "hg38",
    txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene",
    organism = "org.Hs.eg.db")
  expect_error(classify_chromatin_states(db), "empty|No histone")
})

# ============================================================================
# Test 27: Custom regions parameter works
# ============================================================================
test_that("custom regions parameter classifies specified regions", {
  marks <- list(
    h3k4me1 = GRanges("chr1", IRanges(1000, 2000)),
    h3k27ac = GRanges("chr1", IRanges(1000, 2000))
  )
  db <- make_mock_dB(marks)
  # Query a region that overlaps the marks
  custom_gr <- GRanges("chr1", IRanges(500, 2500))
  result <- classify_chromatin_states(db, regions = custom_gr, refine_by_tss = FALSE)
  expect_equal(nrow(result), 1)
  expect_equal(result$chromatin_state[1], "active")
})

# ============================================================================
# Test 28: Custom regions with no overlap = unmarked
# ============================================================================
test_that("custom regions with no mark overlap = unmarked", {
  marks <- list(
    h3k4me1 = GRanges("chr1", IRanges(1000, 2000))
  )
  db <- make_mock_dB(marks)
  # Query a region that does NOT overlap the marks
  custom_gr <- GRanges("chr1", IRanges(50000, 51000))
  result <- classify_chromatin_states(db, regions = custom_gr, refine_by_tss = FALSE)
  expect_equal(nrow(result), 1)
  expect_equal(result$chromatin_state[1], "unmarked")
})

# ============================================================================
# Test 29: h2az_overlap column is present in output
# ============================================================================
test_that("h2az_overlap column present in output", {
  marks <- list(
    h3k4me1 = GRanges("chr1", IRanges(1000, 2000))
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = FALSE)
  expect_true("h2az_overlap" %in% names(result))
})

# ============================================================================
# Test 30: All mark combinations produce valid states
# ============================================================================
test_that("comprehensive mark combinations all produce valid 6-state labels", {
  valid_states <- c("active", "bivalent", "poised", "primed", "repressed", "unmarked")

  # Test each biologically meaningful combination
  combos <- list(
    list(marks = list(h3k4me1 = GRanges("chr1", IRanges(1000, 2000)),
                      h3k27ac = GRanges("chr1", IRanges(1000, 2000))),
         expected = "active"),
    list(marks = list(h3k4me3 = GRanges("chr1", IRanges(1000, 2000)),
                      h3k27ac = GRanges("chr1", IRanges(1000, 2000))),
         expected = "active"),
    list(marks = list(h3k36me3 = GRanges("chr1", IRanges(1000, 2000))),
         expected = "unmarked"),  # gene body transcription, not active regulatory
    list(marks = list(h3k4me3 = GRanges("chr1", IRanges(1000, 2000)),
                      h3k27me3 = GRanges("chr1", IRanges(1000, 2000))),
         expected = "bivalent"),
    list(marks = list(h3k4me1 = GRanges("chr1", IRanges(1000, 2000)),
                      h3k27me3 = GRanges("chr1", IRanges(1000, 2000))),
         expected = "poised"),
    list(marks = list(h3k4me1 = GRanges("chr1", IRanges(1000, 2000))),
         expected = "primed"),
    list(marks = list(h3k27me3 = GRanges("chr1", IRanges(1000, 2000))),
         expected = "repressed"),
    list(marks = list(h3k9me3 = GRanges("chr1", IRanges(1000, 2000))),
         expected = "repressed"),
    list(marks = list(h3k27me3 = GRanges("chr1", IRanges(1000, 2000)),
                      h3k9me3 = GRanges("chr1", IRanges(1000, 2000))),
         expected = "repressed"),
    # H3K27ac alone = active (Creyghton et al. 2010: H3K27ac is sufficient)
    list(marks = list(h3k27ac = GRanges("chr1", IRanges(1000, 2000))),
         expected = "active"),
    # H3K4me3 alone (no activation or repression mark) = primed
    list(marks = list(h3k4me3 = GRanges("chr1", IRanges(1000, 2000))),
         expected = "primed"),
    # EC-1.1: H3K36me3 + H3K27ac = active (gene body WITH active mark)
    list(marks = list(h3k36me3 = GRanges("chr1", IRanges(1000, 2000)),
                      h3k27ac = GRanges("chr1", IRanges(1000, 2000))),
         expected = "active"),
    # H3K36me3 + H3K27me3 = repressed (gene body with repressive mark)
    list(marks = list(h3k36me3 = GRanges("chr1", IRanges(1000, 2000)),
                      h3k27me3 = GRanges("chr1", IRanges(1000, 2000))),
         expected = "repressed"),
    # H3K36me3 + H3K4me1 + H3K27ac = active (strong enhancer/gene body)
    list(marks = list(h3k36me3 = GRanges("chr1", IRanges(1000, 2000)),
                      h3k4me1 = GRanges("chr1", IRanges(1000, 2000)),
                      h3k27ac = GRanges("chr1", IRanges(1000, 2000))),
         expected = "active"),
    # H3K36me3 + H3K4me3 + H3K27ac = active (promoter in gene body)
    list(marks = list(h3k36me3 = GRanges("chr1", IRanges(1000, 2000)),
                      h3k4me3 = GRanges("chr1", IRanges(1000, 2000)),
                      h3k27ac = GRanges("chr1", IRanges(1000, 2000))),
         expected = "active"),
    # H3K36me3 + H3K9me3 = repressed (heterochromatic gene body)
    list(marks = list(h3k36me3 = GRanges("chr1", IRanges(1000, 2000)),
                      h3k9me3 = GRanges("chr1", IRanges(1000, 2000))),
         expected = "repressed")
  )

  for (i in seq_along(combos)) {
    db <- make_mock_dB(combos[[i]]$marks)
    result <- classify_chromatin_states(db, refine_by_tss = FALSE)
    expect_true(all(result$chromatin_state %in% valid_states),
      info = paste0("Combo ", i, ": got ", paste(unique(result$chromatin_state), collapse = ", ")))
    expect_true(combos[[i]]$expected %in% result$chromatin_state,
      info = paste0("Combo ", i, ": expected '", combos[[i]]$expected,
        "' but got '", paste(unique(result$chromatin_state), collapse = ", "), "'"))
  }
})

# ============================================================================
# Test 31: H3K36me3 + H3K27ac = active (gene body with activation)
# ============================================================================
test_that("H3K36me3 + H3K27ac is active (active gene body)", {
  marks <- list(
    h3k36me3 = GRanges("chr1", IRanges(1000, 2000)),
    h3k27ac = GRanges("chr1", IRanges(1000, 2000))
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = FALSE)
  expect_equal(result$chromatin_state[1], "active")
  expect_true(grepl("h3k36me3", result$marks_present[1]))
  expect_true(grepl("h3k27ac", result$marks_present[1]))
})

# ============================================================================
# Test 32: H3K36me3 alone marks_present shows h3k36me3
# ============================================================================
test_that("H3K36me3 alone shows in marks_present even when unmarked", {
  marks <- list(
    h3k36me3 = GRanges("chr1", IRanges(1000, 2000))
  )
  db <- make_mock_dB(marks)
  result <- classify_chromatin_states(db, refine_by_tss = FALSE)
  expect_equal(result$chromatin_state[1], "unmarked")
  expect_true(grepl("h3k36me3", result$marks_present[1]),
    info = "H3K36me3 gene body should appear in marks_present for provenance")
})
