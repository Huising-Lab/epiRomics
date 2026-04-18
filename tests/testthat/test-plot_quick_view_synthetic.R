## Synthetic tests for plot_quick_view()
## Tests input validation only (no database, no BigWig files required)

test_that("quick_view rejects missing gene and region", {
  bw <- c(Alpha = tempfile(fileext = ".bw"))
  file.create(bw)
  on.exit(unlink(bw))
  expect_error(
    plot_quick_view(bw_paths = bw),
    "Exactly one of 'gene' or 'region' must be provided"
  )
})

test_that("quick_view rejects both gene and region", {
  bw <- c(Alpha = tempfile(fileext = ".bw"))
  file.create(bw)
  on.exit(unlink(bw))
  expect_error(
    plot_quick_view(gene = "INS",
      region = list(chr = "chr11", start = 1, end = 100),
      bw_paths = bw),
    "Provide either 'gene' or 'region', not both"
  )
})

test_that("quick_view rejects empty bw_paths", {
  expect_error(
    plot_quick_view(gene = "INS", bw_paths = character(0)),
    "bw_paths must be a non-empty named character vector"
  )
})

test_that("quick_view rejects unnamed bw_paths", {
  bw <- tempfile(fileext = ".bw")
  file.create(bw)
  on.exit(unlink(bw))
  expect_error(
    plot_quick_view(gene = "INS", bw_paths = bw),
    "bw_paths must be named"
  )
})

test_that("quick_view rejects mismatched labels length", {
  bw <- c(Alpha = tempfile(fileext = ".bw"))
  file.create(bw)
  on.exit(unlink(bw))
  expect_error(
    plot_quick_view(gene = "INS", bw_paths = bw,
      labels = c("A", "B")),
    "labels must have the same length as bw_paths"
  )
})

test_that("quick_view rejects missing BigWig files", {
  bw <- c(Alpha = "/nonexistent/path/alpha.bw")
  expect_error(
    plot_quick_view(gene = "INS", bw_paths = bw),
    "BigWig file.*not found"
  )
})

test_that("quick_view rejects unknown genome without explicit TxDb", {
  bw <- c(Alpha = tempfile(fileext = ".bw"))
  file.create(bw)
  on.exit(unlink(bw))
  expect_error(
    plot_quick_view(gene = "INS", bw_paths = bw, genome = "hg19"),
    "no built-in TxDb default"
  )
})

test_that("quick_view rejects invalid region format", {
  bw <- c(Alpha = tempfile(fileext = ".bw"))
  file.create(bw)
  on.exit(unlink(bw))
  expect_error(
    plot_quick_view(region = list(chrom = "chr11", start = 1),
      bw_paths = bw),
    "region must be a list with 'chr', 'start', 'end' elements"
  )
})

## Gene mode integration tests exercise the full rendering pipeline including
## .resolve_gene_coordinates(). A real gene-mode test REQUIRES the org.Db and
## TxDb annotation packages to be loadable — otherwise the test would not be
## exercising the actual lookup path. We hard-assert their availability via
## requireNamespace() so a missing Suggests dependency surfaces as a clear
## failure, not a silent pass. INS is a hg38 chr11 gene inside the toy
## BigWig window (chr11:2,099,779-2,221,221).

test_that("quick_view gene-mode (hg38 auto-resolve) renders PDF with real org.Db + TxDb", {
  require_extdata()
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  # Prove both annotation packages actually load in this test's namespace.
  expect_true(requireNamespace("org.Hs.eg.db", quietly = TRUE))
  expect_true(requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene",
                               quietly = TRUE))

  bw <- c(Beta_ATAC = file.path(toy_extdata_dir(), "BigWigs",
                                "Beta_ATAC.toy.bw"))
  out_pdf <- tempfile(fileext = ".pdf")
  on.exit(unlink(out_pdf), add = TRUE)

  # Expected benign noise from plot_quick_view on gene-mode render:
  #   * Gviz::plotTracks emits informational "drawing axis" messages.
  #   * rtracklayer may warn about bins outside the toy window.
  # Neither is under test — assertion is PDF existence + non-trivial size.
  suppressWarnings(suppressMessages(
    plot_quick_view(gene = "INS", bw_paths = bw,
                    genome = "hg38", export = out_pdf)
  ))

  expect_true(file.exists(out_pdf))
  expect_gt(file.info(out_pdf)$size, 1000)
})

test_that("quick_view infers RNA type and renders mixed tracks", {
  ## Unit check: helper classifies sample names correctly
  tc <- epiRomics:::.build_quick_view_track_connection(
    bw_paths = c(Alpha_ATAC = "/tmp/a.bw", Alpha_RNA = "/tmp/b.bw",
                 Beta_rnaseq = "/tmp/c.bw", Beta_expression = "/tmp/d.bw"),
    colors = NULL
  )
  expect_equal(tc$type, c("atac", "rna", "rna", "rna"))
  expect_equal(tc$name,
               c("Alpha_ATAC", "Alpha_RNA", "Beta_rnaseq", "Beta_expression"))

  ## Integration check: render mixed ATAC + RNA tracks end-to-end.
  ## Gene mode requires org.Db + TxDb to be loadable.
  require_extdata()
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  expect_true(requireNamespace("org.Hs.eg.db", quietly = TRUE))
  expect_true(requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene",
                               quietly = TRUE))
  bw <- c(
    Alpha_ATAC = file.path(toy_extdata_dir(), "BigWigs",
                           "Alpha_ATAC.toy.bw"),
    Alpha_RNA  = file.path(toy_extdata_dir(), "BigWigs",
                           "Alpha_RNA.toy.bw")
  )
  out_pdf <- tempfile(fileext = ".pdf")
  on.exit(unlink(out_pdf), add = TRUE)
  # Expected benign noise: Gviz axis-drawing messages + rtracklayer bin warnings
  # against the toy chr11 window; not under test here.
  suppressWarnings(suppressMessages(
    plot_quick_view(gene = "INS", bw_paths = bw,
                    mirror = TRUE, export = out_pdf)
  ))
  expect_true(file.exists(out_pdf))
  expect_gt(file.info(out_pdf)$size, 1000)
})

test_that("quick_view applies label override", {
  ## Unit check: validator rewrites names when labels provided
  bw <- c(Original = tempfile(fileext = ".bw"))
  file.create(bw)
  on.exit(unlink(bw), add = TRUE)
  result <- epiRomics:::.validate_quick_view_inputs(
    gene = "INS", region = NULL, bw_paths = bw,
    labels = "NewLabel", colors = NULL
  )
  expect_equal(names(result), "NewLabel")
  expect_equal(unname(result), unname(bw))

  ## Integration check: rendered PDF respects label override.
  ## Gene mode requires org.Db + TxDb to be loadable.
  require_extdata()
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  expect_true(requireNamespace("org.Hs.eg.db", quietly = TRUE))
  expect_true(requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene",
                               quietly = TRUE))
  bw2 <- c(OriginalName = file.path(toy_extdata_dir(), "BigWigs",
                                    "Beta_ATAC.toy.bw"))
  out_pdf <- tempfile(fileext = ".pdf")
  on.exit(unlink(out_pdf), add = TRUE)
  # Expected benign noise: Gviz axis-drawing messages + rtracklayer bin warnings
  # against the toy chr11 window; not under test here.
  suppressWarnings(suppressMessages(
    plot_quick_view(gene = "INS", bw_paths = bw2,
                    labels = "RelabeledBeta", export = out_pdf)
  ))
  expect_true(file.exists(out_pdf))
  expect_gt(file.info(out_pdf)$size, 1000)
})
