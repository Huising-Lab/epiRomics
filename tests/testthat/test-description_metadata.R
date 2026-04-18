# Test DESCRIPTION metadata requirements (META-01)
# Guards R version dependency against regression.

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
