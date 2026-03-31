#' Identifies putative enhancer regions utilizing
#' select histone mark co-occurrence
#'
#' @param epiRomics_dB epiRomics class database
#'   containing all data initially loaded
#' @param epiRomics_histone_mark_1 name of first
#'   histone mark, must match name in
#'   epiROmics_dB@meta, default set to h3k4me1
#' @param epiRomics_histone_mark_2 name of second
#'   histone mark, must match name in
#'   epiROmics_dB@meta default set to h3k27ac
#' @return Variable of class epiRomics further
#'   exploring candidate enhancer regions identified
#'   after histone integration
#' @export
#' @examples
#' db <- epiRomicsS4(annotations = GenomicRanges::GRanges(),
#'   meta = data.frame(name = character(), type = character(),
#'     file = character(), stringsAsFactors = FALSE),
#'   genome = "hg38")
#' tryCatch(
#'   epiRomics_enhancers_co_marks(db),
#'   error = function(e) message(e$message)
#' )
#' \donttest{
#' enhancers <- epiRomics_enhancers_co_marks(epiRomics_dB)
#' }
epiRomics_enhancers_co_marks <- function(
    epiRomics_dB,
    epiRomics_histone_mark_1 = "h3k4me1",
    epiRomics_histone_mark_2 = "h3k27ac") {
  # Parameter validation
  validate_epiRomics_dB(
    epiRomics_dB, "epiRomics_enhancers_co_marks"
  )
  validate_character_param(
    epiRomics_histone_mark_1,
    "epiRomics_histone_mark_1",
    "epiRomics_enhancers_co_marks"
  )
  validate_character_param(
    epiRomics_histone_mark_2,
    "epiRomics_histone_mark_2",
    "epiRomics_enhancers_co_marks"
  )

  # Check if histone marks exist in the database
  # H2A.Z is type "histone" (not separate variant type), so users can call
  # epiRomics_enhancers(dB, "h2az", "h3k27ac") for H2A.Z-based enhancer ID
  available_marks <- epiRomics_dB@meta$name[epiRomics_dB@meta$type %in%
    base::c("histone", "chip")]
  if (!(epiRomics_histone_mark_1 %in% available_marks)) {
    mark_fmt <- paste0(
      "epiRomics_enhancers_co_marks: ",
      "histone mark '%s' not found in ",
      "database. Available marks: %s"
    )
    base::stop(base::sprintf(
      mark_fmt,
      epiRomics_histone_mark_1,
      base::paste(available_marks, collapse = ", ")
    ))
  }
  if (!(epiRomics_histone_mark_2 %in% available_marks)) {
    mark_fmt <- paste0(
      "epiRomics_enhancers_co_marks: ",
      "histone mark '%s' not found in ",
      "database. Available marks: %s"
    )
    base::stop(base::sprintf(
      mark_fmt,
      epiRomics_histone_mark_2,
      base::paste(available_marks, collapse = ", ")
    ))
  }

  tryCatch(
    {
      epiRomics_putative <- epiRomics_dB
      mark1 <- base::paste0(
        epiRomics_dB@genome, "_custom_",
        epiRomics_histone_mark_1
      )
      mark2 <- base::paste0(
        epiRomics_dB@genome, "_custom_",
        epiRomics_histone_mark_2
      )

      # Check if the marks exist in annotations
      mark1_exists <- base::any(epiRomics_dB@annotations$type == mark1)
      mark2_exists <- base::any(epiRomics_dB@annotations$type == mark2)

      if (!mark1_exists) {
        base::stop(base::sprintf(
          paste0(
            "epiRomics_enhancers_co_marks: ",
            "No annotations found for mark '%s'"
          ),
          mark1
        ))
      }
      if (!mark2_exists) {
        base::stop(base::sprintf(
          paste0(
            "epiRomics_enhancers_co_marks: ",
            "No annotations found for mark '%s'"
          ),
          mark2
        ))
      }

      mark1_gr <- GenomicRanges::reduce(
        epiRomics_dB@annotations[
          epiRomics_dB@annotations$type == mark1, ]
      )
      mark2_gr <- GenomicRanges::reduce(
        epiRomics_dB@annotations[
          epiRomics_dB@annotations$type == mark2, ]
      )
      epiRomics_putative@annotations <-
        GenomicRanges::intersect(mark1_gr, mark2_gr)

      if (base::length(epiRomics_putative@annotations) == 0) {
        base::warning(paste0(
          "epiRomics_enhancers_co_marks: ",
          "No overlapping regions found ",
          "between the two histone marks"
        ))
      }

      base::return(epiRomics_putative)
    },
    error = function(e) {
      base::stop(base::sprintf("epiRomics_enhancers_co_marks: %s", e$message))
    }
  )
}
