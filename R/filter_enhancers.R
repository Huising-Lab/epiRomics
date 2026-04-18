#' Filters putative enhancers called by
#' epiRomics_enhancers by crossing against curated
#' FANTOM data
#'
#' @param putative_enhancers epiRomics class
#'   database containing putative enhancer calls
#' @param database epiRomics class database
#'   containing all data initially loaded
#' @param type epiRomics reference containing
#'   database to validate putative enhancers against.
#'   Default NULL derives from genome:
#'   \code{paste0(genome, "_custom_fantom")}.
#' @return Variable of class epiRomics with filtered candidate enhancer regions
#' @export
#' @examples
#' db <- make_example_database()
#' eso <- make_example_enhanceosome(db)
#' ## filter_enhancers requires fantom annotations in database
#' tryCatch(
#'   filter_enhancers(eso, db),
#'   error = function(e) message(e$message)
#' )
filter_enhancers <- function(
    putative_enhancers,
    database,
    type = NULL) {
  # Parameter validation
  validate_database(
    putative_enhancers,
    "filter_enhancers"
  )
  validate_database(database, "filter_enhancers")
  if (base::is.null(type)) {
    type <- base::paste0(
      database@genome, "_custom_fantom"
    )
  }
  validate_character_param(
    type, "type",
    "filter_enhancers"
  )

  # Check if putative enhancers have annotations
  if (base::length(putative_enhancers@annotations) == 0) {
    base::stop(
      "filter_enhancers: ",
      "putative_enhancers ",
      "has no annotations"
    )
  }

  tryCatch(
    {
      enhanceosome_raw_calls <-
        putative_enhancers@annotations
      epiRomics_functional <- database@annotations

      # Check if the filter type exists in the database
      filter_annotations <- epiRomics_functional[
        epiRomics_functional$type == type, ]
      if (base::length(filter_annotations) == 0) {
        fmt <- "filter_enhancers: No annotations found for filter type '%s'"
        base::stop(base::sprintf(fmt, type))
      }

      # Filtered output object
      putative_enhancers_filtered <- putative_enhancers
      putative_enhancers_filtered@annotations <-
        IRanges::subsetByOverlaps(
          enhanceosome_raw_calls,
          filter_annotations
        )

      # Check if filtering resulted in any data
      n_filtered <- base::length(
        putative_enhancers_filtered@annotations
      )
      if (n_filtered == 0) {
        fmt <- paste0("filter_enhancers: No overlapping regions found ",
                      "with filter type '%s'")
        base::warning(base::sprintf(fmt, type))
      }

      base::return(putative_enhancers_filtered)
    },
    error = function(e) {
      base::stop(base::sprintf("filter_enhancers: %s", e$message))
    }
  )
}
