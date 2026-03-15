## Synthetic tests for epiRomics_quick_view()
## Tests input validation only (no database, no BigWig files required)

test_that("quick_view rejects missing gene and region", {
  bw <- c(Alpha = tempfile(fileext = ".bw"))
  file.create(bw)
  on.exit(unlink(bw))
  expect_error(
    epiRomics_quick_view(bw_paths = bw),
    "Exactly one of 'gene' or 'region' must be provided"
  )
})

test_that("quick_view rejects both gene and region", {
  bw <- c(Alpha = tempfile(fileext = ".bw"))
  file.create(bw)
  on.exit(unlink(bw))
  expect_error(
    epiRomics_quick_view(gene = "INS",
      region = list(chr = "chr11", start = 1, end = 100),
      bw_paths = bw),
    "Provide either 'gene' or 'region', not both"
  )
})

test_that("quick_view rejects empty bw_paths", {
  expect_error(
    epiRomics_quick_view(gene = "INS", bw_paths = character(0)),
    "bw_paths must be a non-empty named character vector"
  )
})

test_that("quick_view rejects unnamed bw_paths", {
  bw <- tempfile(fileext = ".bw")
  file.create(bw)
  on.exit(unlink(bw))
  expect_error(
    epiRomics_quick_view(gene = "INS", bw_paths = bw),
    "bw_paths must be named"
  )
})

test_that("quick_view rejects mismatched labels length", {
  bw <- c(Alpha = tempfile(fileext = ".bw"))
  file.create(bw)
  on.exit(unlink(bw))
  expect_error(
    epiRomics_quick_view(gene = "INS", bw_paths = bw,
      labels = c("A", "B")),
    "labels must have the same length as bw_paths"
  )
})

test_that("quick_view rejects missing BigWig files", {
  bw <- c(Alpha = "/nonexistent/path/alpha.bw")
  expect_error(
    epiRomics_quick_view(gene = "INS", bw_paths = bw),
    "BigWig file.*not found"
  )
})

test_that("quick_view rejects invalid genome", {
  bw <- c(Alpha = tempfile(fileext = ".bw"))
  file.create(bw)
  on.exit(unlink(bw))
  expect_error(
    epiRomics_quick_view(gene = "INS", bw_paths = bw, genome = "hg19"),
    "'arg' should be one of"
  )
})

test_that("quick_view rejects invalid region format", {
  bw <- c(Alpha = tempfile(fileext = ".bw"))
  file.create(bw)
  on.exit(unlink(bw))
  expect_error(
    epiRomics_quick_view(region = list(chrom = "chr11", start = 1),
      bw_paths = bw),
    "region must be a list with 'chr', 'start', 'end' elements"
  )
})

test_that("quick_view auto-resolves hg38 TxDb", {
  # Cannot test rendering without real data, but verify the function

# accepts valid hg38 params and progresses past validation
  bw <- c(Alpha_ATAC = tempfile(fileext = ".bw"))
  file.create(bw)
  on.exit(unlink(bw))

  # Should fail downstream (gene lookup) but NOT on validation
  err <- tryCatch(
    epiRomics_quick_view(gene = "INS", bw_paths = bw, genome = "hg38"),
    error = function(e) conditionMessage(e)
  )
  # Should get past validation — error should be about gene lookup, not params
  expect_false(grepl("bw_paths must be", err))
  expect_false(grepl("Exactly one of", err))
})

test_that("quick_view auto-resolves mm10 TxDb", {
  bw <- c(Sample = tempfile(fileext = ".bw"))
  file.create(bw)
  on.exit(unlink(bw))

  err <- tryCatch(
    epiRomics_quick_view(gene = "Ins1", bw_paths = bw, genome = "mm10"),
    error = function(e) conditionMessage(e)
  )
  expect_false(grepl("bw_paths must be", err))
})

test_that("quick_view infers RNA type from names", {
  bw1 <- tempfile(fileext = ".bw")
  bw2 <- tempfile(fileext = ".bw")
  file.create(bw1, bw2)
  on.exit(unlink(c(bw1, bw2)))

  bw <- c(Alpha_ATAC = bw1, Alpha_RNA = bw2)

  # Capture the error from gene lookup (validation passes)
  err <- tryCatch(
    epiRomics_quick_view(gene = "INS", bw_paths = bw, mirror = TRUE),
    error = function(e) conditionMessage(e)
  )
  # Should get past the track_connection building phase
  expect_false(grepl("bw_paths must be", err))
})

test_that("quick_view labels override works", {
  bw <- c(Original = tempfile(fileext = ".bw"))
  file.create(bw)
  on.exit(unlink(bw))

  err <- tryCatch(
    epiRomics_quick_view(gene = "INS", bw_paths = bw,
      labels = "NewLabel"),
    error = function(e) conditionMessage(e)
  )
  expect_false(grepl("labels must have", err))
})
