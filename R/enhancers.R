#' Identifies putative enhancer regions utilizing
#' select histone mark co-occurrence
#'
#' @param database epiRomics class database
#'   containing all data initially loaded
#' @param histone_mark_1 name of first
#'   histone mark, must match name in
#'   epiROmics_dB@meta, default set to h3k4me1
#' @param histone_mark_2 name of second
#'   histone mark, must match name in
#'   epiROmics_dB@meta default set to h3k27ac
#' @return Variable of class epiRomics further
#'   exploring candidate enhancer regions identified
#'   after histone integration
#' @export
#' @examples
#' db <- make_example_database()
#' enhancers <- find_enhancers_by_comarks(db)
#' length(methods::slot(enhancers, "annotations"))
find_enhancers_by_comarks <- function(
    database,
    histone_mark_1 = "h3k4me1",
    histone_mark_2 = "h3k27ac") {
  # Parameter validation
  validate_database(
    database, "find_enhancers_by_comarks"
  )
  validate_character_param(
    histone_mark_1,
    "histone_mark_1",
    "find_enhancers_by_comarks"
  )
  validate_character_param(
    histone_mark_2,
    "histone_mark_2",
    "find_enhancers_by_comarks"
  )

  # Check if histone marks exist in the database
  # H2A.Z is type "histone" (not separate variant type), so users can call
  # epiRomics_enhancers(dB, "h2az", "h3k27ac") for H2A.Z-based enhancer ID
  available_marks <- database@meta$name[database@meta$type %in%
    base::c("histone", "chip")]
  if (!(histone_mark_1 %in% available_marks)) {
    mark_fmt <- paste0(
      "find_enhancers_by_comarks: ",
      "histone mark '%s' not found in ",
      "database. Available marks: %s"
    )
    base::stop(base::sprintf(
      mark_fmt,
      histone_mark_1,
      base::paste(available_marks, collapse = ", ")
    ))
  }
  if (!(histone_mark_2 %in% available_marks)) {
    mark_fmt <- paste0(
      "find_enhancers_by_comarks: ",
      "histone mark '%s' not found in ",
      "database. Available marks: %s"
    )
    base::stop(base::sprintf(
      mark_fmt,
      histone_mark_2,
      base::paste(available_marks, collapse = ", ")
    ))
  }

  tryCatch(
    {
      epiRomics_putative <- database
      mark1 <- base::paste0(
        database@genome, "_custom_",
        histone_mark_1
      )
      mark2 <- base::paste0(
        database@genome, "_custom_",
        histone_mark_2
      )

      # Check if the marks exist in annotations
      mark1_exists <- base::any(database@annotations$type == mark1)
      mark2_exists <- base::any(database@annotations$type == mark2)

      fmt <- "find_enhancers_by_comarks: No annotations found for mark '%s'"
      if (!mark1_exists) {
        base::stop(base::sprintf(fmt, mark1))
      }
      if (!mark2_exists) {
        base::stop(base::sprintf(fmt, mark2))
      }

      mark1_gr <- GenomicRanges::reduce(
        database@annotations[
          database@annotations$type == mark1, ]
      )
      mark2_gr <- GenomicRanges::reduce(
        database@annotations[
          database@annotations$type == mark2, ]
      )
      epiRomics_putative@annotations <-
        GenomicRanges::intersect(mark1_gr, mark2_gr)

      if (base::length(epiRomics_putative@annotations) == 0) {
        base::warning(
          "find_enhancers_by_comarks: ",
          "No overlapping regions found ",
          "between the two histone marks"
        )
      }

      base::return(epiRomics_putative)
    },
    error = function(e) {
      base::stop(base::sprintf("find_enhancers_by_comarks: %s", e$message))
    }
  )
}
