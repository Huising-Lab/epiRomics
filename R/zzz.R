# Suppress R CMD check NOTEs for non-standard evaluation column references
utils::globalVariables(c("broad", "broad_score", "detail", "score", "ws"))

#' Package startup message
#'
#' Displays version, cache status, and citation info on
#' package attach.
#'
#' @param libname character. Library path.
#' @param pkgname character. Package name.
#' @return Called for side effect (startup message).
#' @noRd
.onAttach <- function(libname, pkgname) {
  ver <- utils::packageVersion("epiRomics")

  if (epiRomics_has_cache()) {
    data_msg <- "Example data: cached \u2713"
  } else {
    data_msg <- paste0(
      "Run epiRomics_cache_data() to download ",
      "example datasets (~1.3 GB).")
  }

  base::packageStartupMessage(
    "epiRomics v", ver, " \u2014 Epigenomic Analysis Package\n",
    data_msg, "\n",
    "Citation: citation(\"epiRomics\")"
  )
}
