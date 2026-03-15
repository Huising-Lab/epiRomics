library(testthat)
library(epiRomics)

# Test: epiRomicsS4 class definition and structure
testthat::test_that("epiRomicsS4 S4 class is defined correctly", {
  # Check class existence
  testthat::expect_true(isClass("epiRomicsS4"))

  # Get class slots
  slots <- getSlots("epiRomicsS4")

  # Check slot names
  testthat::expect_true(all(c("annotations", "meta", "txdb", "organism", "genome") %in% names(slots)))

  # Check slot types
  testthat::expect_equal(slots[["annotations"]], "GRanges")
  testthat::expect_equal(slots[["meta"]], "data.frame")
  testthat::expect_equal(slots[["txdb"]], "character")
  testthat::expect_equal(slots[["organism"]], "character")
  testthat::expect_equal(slots[["genome"]], "character")

  # Test class instantiation
  obj <- methods::new("epiRomicsS4")
  testthat::expect_s4_class(obj, "epiRomicsS4")
  testthat::expect_true(inherits(obj@annotations, "GRanges"))
  testthat::expect_true(inherits(obj@meta, "data.frame"))
  testthat::expect_true(is.character(obj@txdb))
  testthat::expect_true(is.character(obj@organism))
  testthat::expect_true(is.character(obj@genome))
})
