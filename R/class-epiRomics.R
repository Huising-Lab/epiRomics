#' An S4 class to manage epiRomics databases and downstream results
#'
#' @slot annotations GRanges
#' @slot meta data.frame
#' @slot txdb txdb string name
#' @slot organism org.db string name
#' @slot genome genome name, e.g. 'mm10' or 'hg38'
#' @export
epiRomicsS4 <- methods::setClass(Class = "epiRomicsS4", slots = c(
  annotations = "GRanges", meta = "data.frame",
  txdb = "character", organism = "character", genome = "character"
), prototype = list(
  annotations = GenomicRanges::GRanges(),
  meta = base::data.frame(), txdb = base::character(), organism = base::character(), genome = base::character()
))
