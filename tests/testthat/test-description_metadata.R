# Test DESCRIPTION metadata requirements (META-01, META-02)
# Guards R version dependency and funder acknowledgments
# against regression

test_that("META-01: R version dependency is R (>= 4.5.0)", {
  desc_path <- file.path(
    testthat::test_path("..", ".."), "DESCRIPTION"
  )
  if (!file.exists(desc_path)) {
    desc_path <- system.file(
      "DESCRIPTION", package = "epiRomics"
    )
  }
  expect_true(
    file.exists(desc_path),
    info = "DESCRIPTION file must exist"
  )

  dcf <- base::read.dcf(desc_path, fields = "Depends")
  depends_field <- dcf[1, "Depends"]

  expect_true(
    base::grepl("R \\(>= 4\\.5\\.0\\)", depends_field),
    info = base::paste0(
      "Depends must contain 'R (>= 4.5.0)', got: ",
      depends_field
    )
  )
})

test_that("META-02: exactly 6 funder entries with role fnd", {
  desc_path <- file.path(
    testthat::test_path("..", ".."), "DESCRIPTION"
  )
  if (!file.exists(desc_path)) {
    desc_path <- system.file(
      "DESCRIPTION", package = "epiRomics"
    )
  }

  desc_lines <- base::readLines(desc_path)
  fnd_count <- base::length(
    base::grep('role = "fnd"', desc_lines)
  )
  expect_equal(fnd_count, 6L)
})

test_that("META-02: all grant identifiers in funder entries", {
  desc_path <- file.path(
    testthat::test_path("..", ".."), "DESCRIPTION"
  )
  if (!file.exists(desc_path)) {
    desc_path <- system.file(
      "DESCRIPTION", package = "epiRomics"
    )
  }

  desc_text <- base::paste(
    base::readLines(desc_path), collapse = "\n"
  )

  # 8 grant identifiers from REQUIREMENTS.md
  grant_ids <- c(
    "NIDDK110276",
    "CDA-2-2013-54",
    "2-SRA-2021-1054-M-N",
    "Sims",
    "Machine Learning",
    "CA093373",
    "OD018223"
  )
  for (gid in grant_ids) {
    expect_true(
      base::grepl(gid, desc_text, fixed = TRUE),
      info = base::paste0(
        "Grant '", gid,
        "' must appear in DESCRIPTION"
      )
    )
  }

  # Case-insensitive check for start-up funds
  expect_true(
    base::grepl("start-up", desc_text, ignore.case = TRUE),
    info = "Start-up funds must appear in DESCRIPTION"
  )
})
