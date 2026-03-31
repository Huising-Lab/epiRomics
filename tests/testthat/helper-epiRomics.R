# Helper to create a working db_sheet with full paths from the relative-path CSV
# testthat sources helper-*.R files before running tests

# Package-local cache environment (avoids global env anti-pattern per TEST-02)
.test_cache <- base::new.env(parent = base::emptyenv())

#' Check if example data is available (inst/extdata present)
#' @return TRUE if extdata is available, FALSE otherwise
has_extdata <- function() {
  extdata_dir <- system.file("extdata", package = "epiRomics")
  nzchar(extdata_dir) && dir.exists(extdata_dir) &&
    file.exists(file.path(extdata_dir, "example_epiRomics_Db_sheet.csv"))
}

#' Skip test if example data is not available
#' @details Called inside test_that() blocks to skip integration tests
#'   when extdata is not bundled (e.g., during R CMD check on built package).
skip_if_no_extdata <- function() {
  testthat::skip_if(!has_extdata(),
                    "Example data not available (inst/extdata excluded from build)")
}

#' Create a temp CSV with resolved full paths for epiRomics_build_dB
#' @return Path to temp CSV with absolute paths, or NULL if extdata missing
make_db_sheet <- function() {
  if (!has_extdata()) return(NULL)
  extdata_dir <- system.file("extdata", package = "epiRomics")
  db_csv <- file.path(extdata_dir, "example_epiRomics_Db_sheet.csv")
  df <- utils::read.csv(db_csv, stringsAsFactors = FALSE)
  df$path <- file.path(extdata_dir, df$path)
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
  suppressWarnings({
    dB <- epiRomics::epiRomics_build_dB(
      epiRomics_db_file = db_sheet,
      txdb_organism = "TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene",
      epiRomics_genome = "hg38",
      epiRomics_organism = "org.Hs.eg.db"
    )
  })
  base::assign(".test_epiRomics_dB", dB, envir = .test_cache)
  dB
}

#' Load track connection data frame with full paths
#' @return data.frame with path, name, color, type columns, or NULL if missing
make_track_connection <- function() {
  if (!has_extdata()) return(NULL)
  extdata_dir <- system.file("extdata", package = "epiRomics")
  bw_csv <- file.path(extdata_dir, "example_epiRomics_BW_sheet.csv")
  tc <- utils::read.csv(bw_csv, stringsAsFactors = FALSE)
  tc$path <- file.path(extdata_dir, tc$path)
  tc
}
