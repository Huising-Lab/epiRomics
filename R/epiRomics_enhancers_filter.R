#' Filters putative enhancers called by epiRomics_enhancers by crossing against curated FANTOM data
#'
#' @param epiRomics_putative_enhancers epiRomics class database containing putative enhancer calls
#' @param epiRomics_dB epiRomics class database containing all data initially loaded
#' @param epiRomics_type epiRomics reference containing database to validate putative enhancers against
#' @return Variable of class epiRomics with filtered candidate enhancer regions
#' @export
#' @examples
#' db <- epiRomicsS4(annotations = GenomicRanges::GRanges(),
#'   meta = data.frame(name = character(), type = character(),
#'     file = character(), stringsAsFactors = FALSE),
#'   genome = "hg38")
#' tryCatch(epiRomics_enhancers_filter(db, db), error = function(e) message(e$message))
#' \donttest{
#' filtered <- epiRomics_enhancers_filter(epiRomics_putative_enhancers, epiRomics_dB)
#' }
epiRomics_enhancers_filter <- function(epiRomics_putative_enhancers, epiRomics_dB, epiRomics_type = "mm10_custom_fantom") {
  # Parameter validation
  validate_epiRomics_dB(epiRomics_putative_enhancers, "epiRomics_enhancers_filter")
  validate_epiRomics_dB(epiRomics_dB, "epiRomics_enhancers_filter")
  validate_character_param(epiRomics_type, "epiRomics_type", "epiRomics_enhancers_filter")

  # Check if putative enhancers have annotations
  if (base::length(epiRomics_putative_enhancers@annotations) == 0) {
    base::stop("epiRomics_enhancers_filter: epiRomics_putative_enhancers has no annotations")
  }

  tryCatch(
    {
      epiRomics_enhanceosome_raw_calls <- epiRomics_putative_enhancers@annotations
      epiRomics_functional <- epiRomics_dB@annotations

      # Check if the filter type exists in the database
      filter_annotations <- epiRomics_functional[epiRomics_functional$type == epiRomics_type, ]
      if (base::length(filter_annotations) == 0) {
        base::stop(base::sprintf("epiRomics_enhancers_filter: No annotations found for filter type '%s'", epiRomics_type))
      }

      # Filtered output object
      epiRomics_putative_enhancers_filtered <- epiRomics_putative_enhancers
      epiRomics_putative_enhancers_filtered@annotations <- IRanges::subsetByOverlaps(
        epiRomics_enhanceosome_raw_calls,
        filter_annotations
      )

      # Check if filtering resulted in any data
      if (base::length(epiRomics_putative_enhancers_filtered@annotations) == 0) {
        base::warning(base::sprintf("epiRomics_enhancers_filter: No overlapping regions found with filter type '%s'", epiRomics_type))
      }

      base::return(epiRomics_putative_enhancers_filtered)
    },
    error = function(e) {
      base::stop(base::sprintf("epiRomics_enhancers_filter: %s", e$message))
    }
  )
}
