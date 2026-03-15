#' Visualizes data from epiRomics results (deprecated)
#'
#' This function is deprecated. Use \code{\link{epiRomics_track_layer}} instead,
#' which automatically detects the genome from the database object and works
#' for both human and mouse data.
#'
#' @inheritParams epiRomics_track_layer
#' @return GViz plot
#' @export
#' @examples
#' db <- epiRomicsS4(annotations = GenomicRanges::GRanges(),
#'   meta = data.frame(name = character(), type = character(),
#'     file = character(), stringsAsFactors = FALSE),
#'   genome = "hg38")
#' tryCatch(epiRomics_track_layer_human(db, 1, db,
#'   data.frame(file = character(), name = character(),
#'     color = character(), type = character())),
#'   error = function(e) message(e$message))
#' \donttest{
#' epiRomics_track_layer_human("INS", epiRomics_dB, epiRomics_enhancers,
#'   epiRomics_track_connection)
#' }
epiRomics_track_layer_human <- function(epiRomics_putative_enhanceosome, epiRomics_index, epiRomics_dB, epiRomics_track_connection, epiRomics_keep_epitracks = TRUE) {
  .Deprecated("epiRomics_track_layer",
    msg = "epiRomics_track_layer_human() is deprecated. Use epiRomics_track_layer() which auto-detects genome."
  )
  epiRomics_track_layer(epiRomics_putative_enhanceosome, epiRomics_index, epiRomics_dB, epiRomics_track_connection, epiRomics_keep_epitracks)
}
