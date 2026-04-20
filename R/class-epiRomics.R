#' An S4 class to manage epiRomics databases and downstream results
#'
#' The \code{epiRomicsS4} class holds a genomic annotation set along
#' with the metadata needed to interpret it (data-source catalog,
#' TxDb identifier, organism annotation package, genome assembly).
#' Every slot is exposed through a public getter and setter — users
#' should access slot contents through those accessors rather than
#' \code{obj@slot} or \code{methods::slot(obj, "slot")}.
#'
#' @slot annotations GRanges object containing genomic annotations.
#'   Access via \code{\link{annotations}()}.
#' @slot meta data.frame containing metadata about loaded data sources.
#'   Access via \code{\link{meta}()}.
#' @slot txdb character string of TxDb package::object name.
#'   Access via \code{\link{txdb}()}.
#' @slot organism character string of org.db package name.
#'   Access via \code{\link{organism}()}.
#' @slot genome character string of genome assembly name
#'   (e.g., 'mm10', 'hg38'). Access via \code{\link{genome}()}.
#' @return An S4 class definition for epiRomicsS4
#' @seealso
#'   Public accessors for this class:
#'   \code{\link{annotations}},
#'   \code{\link[=annotations<-]{annotations<-}},
#'   \code{\link{meta}},
#'   \code{\link[=meta<-]{meta<-}},
#'   \code{\link{txdb}},
#'   \code{\link[=txdb<-]{txdb<-}},
#'   \code{\link{organism}},
#'   \code{\link[=organism<-]{organism<-}},
#'   \code{\link{genome}},
#'   \code{\link[=genome<-]{genome<-}}.
#'   Overview: \code{\link{epiRomicsS4-accessors}}.
#' @examples
#' showClass("epiRomicsS4")
#' db <- make_example_database()
#' genome(db)
#' organism(db)
#' length(annotations(db))
#' @export
#' @importFrom methods setClass setValidity setMethod setGeneric
epiRomicsS4 <- methods::setClass(Class = "epiRomicsS4", slots = c(
  annotations = "GRanges", meta = "data.frame",
  txdb = "character", organism = "character", genome = "character"
), prototype = list(
  annotations = GenomicRanges::GRanges(),
  meta = base::data.frame(),
  txdb = base::character(),
  organism = base::character(),
  genome = base::character()
))

methods::setValidity("epiRomicsS4", function(object) {
  errors <- base::character()
  if (base::length(object@genome) > 0 &&
      !base::all(base::nchar(object@genome) > 0)) {
    errors <- base::c(errors, "genome must be a non-empty character string")
  }
  if (base::length(object@txdb) > 0 &&
      !base::all(base::nchar(object@txdb) > 0)) {
    errors <- base::c(errors, "txdb must be a non-empty character string")
  }
  if (base::length(errors) == 0) TRUE else errors
})

methods::setMethod("show", "epiRomicsS4", function(object) {
  base::cat("epiRomicsS4 object\n")
  base::cat(
    "  Genome:",
    if (base::length(object@genome) > 0) {
      object@genome
    } else {
      "(not set)"
    }, "\n")
  base::cat("  Annotations:", base::length(object@annotations), "ranges\n")
  base::cat("  Meta:", base::nrow(object@meta), "rows\n")
  base::cat(
    "  Organism:",
    if (base::length(object@organism) > 0) {
      object@organism
    } else {
      "(not set)"
    }, "\n")
  base::cat(
    "  TxDb:",
    if (base::length(object@txdb) > 0) {
      object@txdb
    } else {
      "(not set)"
    }, "\n")
})
