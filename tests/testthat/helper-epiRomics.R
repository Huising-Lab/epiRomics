# Helper to create a working db_sheet with full paths from the relative-path CSV
# testthat sources helper-*.R files before running tests

# Package-local cache environment (avoids global env anti-pattern per TEST-02)
.test_cache <- base::new.env(parent = base::emptyenv())

#' Locate the shipped toy data directory (inst/extdata/toy/)
#' @return Absolute path (character) to the toy data directory; "" if missing.
toy_extdata_dir <- function() {
  system.file("extdata", "toy", package = "epiRomics")
}

#' Check if example toy data is available
#' @return TRUE if toy extdata is shipped, FALSE otherwise
has_extdata <- function() {
  toy_dir <- toy_extdata_dir()
  nzchar(toy_dir) && dir.exists(toy_dir) &&
    file.exists(file.path(toy_dir, "example_epiRomics_Db_sheet.csv"))
}

#' Require toy example data to be available (hard assertion).
#' @details Called inside test_that() blocks. Toy data ships with the
#'   package tarball (inst/extdata/toy/) so its absence is a packaging
#'   regression that MUST fail loudly rather than silently skip — any
#'   Bioconductor CI run that ships a tarball without toy/ needs to
#'   break the build, not emit a green skip.
require_extdata <- function() {
  testthat::expect_true(has_extdata(),
    info = "inst/extdata/toy must ship with the package tarball")
  base::stopifnot(has_extdata())
}

#' Create a temp CSV with resolved full paths for build_database
#' @return Path to temp CSV with absolute paths, or NULL if toy extdata missing
make_db_sheet <- function() {
  if (!has_extdata()) return(NULL)
  toy_dir <- toy_extdata_dir()
  db_csv <- file.path(toy_dir, "example_epiRomics_Db_sheet.csv")
  df <- utils::read.csv(db_csv, stringsAsFactors = FALSE)
  df$path <- file.path(toy_dir, df$path)
  tmp <- tempfile(fileext = ".csv")
  utils::write.csv(df, tmp, row.names = FALSE)
  tmp
}

#' Build the standard test database (cached per session)
#' @return epiRomicsS4 object, or NULL if extdata missing
build_test_dB <- function() {
  if (base::exists(".test_epiRomics_dB", envir = .test_cache)) {
    return(base::get(".test_epiRomics_dB", envir = .test_cache))
  }
  db_sheet <- make_db_sheet()
  if (is.null(db_sheet)) return(NULL)
  # Expected benign warnings from build_database() during fixture construction:
  #   * annotatr builtin-annotation fetching emits "No seqlevels in common"
  #     against the toy chr11-only subset.
  #   * ChIPseeker emits "This peak has no overlap with gene" for toy peaks
  #     outside known gene bodies.
  # Neither is under test here; the caller is producing a fixture dB. Warnings
  # from package code under test are instead captured by tests that assert
  # specific behavior (e.g. test-filter_accessible_synthetic.R signal mode).
  suppressWarnings({
    dB <- epiRomics::build_database(
      db_file = db_sheet,
      txdb_organism = "TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene",
      genome = "hg38",
      organism = "org.Hs.eg.db"
    )
  })
  base::assign(".test_epiRomics_dB", dB, envir = .test_cache)
  dB
}

#' Load track connection data frame with full paths (toy BigWigs)
#' @return data.frame with path, name, color, type columns, or NULL if missing
make_track_connection <- function() {
  if (!has_extdata()) return(NULL)
  toy_dir <- toy_extdata_dir()
  bw_csv <- file.path(toy_dir, "example_epiRomics_BW_sheet.csv")
  tc <- utils::read.csv(bw_csv, stringsAsFactors = FALSE)
  tc$path <- file.path(toy_dir, tc$path)
  tc
}
