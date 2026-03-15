#' @useDynLib epiRomics, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"

# Import annotation packages to satisfy CRAN requirements
#' @importFrom TxDb.Hsapiens.UCSC.hg38.knownGene TxDb.Hsapiens.UCSC.hg38.knownGene
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom TxDb.Mmusculus.UCSC.mm10.knownGene TxDb.Mmusculus.UCSC.mm10.knownGene
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @importFrom utils download.file untar head

# The following block is used by usethis to automatically manage roxygen namespace tags. Modify with
# care!  usethis namespace: start usethis namespace: end
NULL
