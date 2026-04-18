#' Compute overlap fractions between histone marks and a curated reference
#'
#' Splits annotations by type, prepares the reduced curated reference,
#' then for each histone mark computes the overlap fraction against the
#' reference. Parallelises when marks > 2 and parallel is available.
#' Also warns for marks with zero annotations.
#'
#' @param database epiRomics class database.
#' @param histone Character vector of histone mark names.
#' @param curated_database Character string of curated database name.
#' @return data.frame with columns Histone_Mark and Fraction_of_Overlap,
#'   sorted by fraction descending.
#' @noRd
.compute_overlap_fractions <- function(
    database,
    histone,
    curated_database) {
  # Pre-split annotations by type for O(1) lookup
  annot_by_type <- base::split(
    database@annotations,
    database@annotations$type
  )

  # Calculate curated database count once
  mark2 <- base::paste0(
    database@genome, "_custom_",
    curated_database
  )
  mark2_annotations <- annot_by_type[[mark2]]
  if (base::is.null(mark2_annotations)) {
    mark2_annotations <- database@annotations[0]
  }
  mark2_count <- base::length(mark2_annotations)

  if (mark2_count == 0) {
    msg <- paste0(
      "benchmark_enhancer_predictor: ",
      "No annotations found for curated database '%s'"
    )
    base::stop(base::sprintf(msg, mark2))
  }

  # Pre-compute reduce(mark2) ONCE — this is constant across all marks.
  reduced_mark2 <- GenomicRanges::reduce(mark2_annotations)

  # Pre-build all mark1 keys
  genome_prefix <- base::paste0(database@genome, "_custom_")
  mark1_keys <- base::paste0(genome_prefix, histone)

  # Overlap computation per mark
  .compute_overlap <- function(idx) {
    mark1_annots <- annot_by_type[[mark1_keys[idx]]]
    if (base::is.null(mark1_annots) || base::length(mark1_annots) == 0L) {
      return(0)
    }
    overlap_count <- base::length(
      GenomicRanges::intersect(
        GenomicRanges::reduce(mark1_annots),
        reduced_mark2
      )
    )
    overlap_count / mark2_count
  }

  n_marks <- base::length(histone)
  mark_indices <- base::seq_len(n_marks)

  if (n_marks > 2L &&
      .Platform$OS.type == "unix" &&
      base::requireNamespace(
        "parallel", quietly = TRUE
      )) {
    n_cores <- .detect_cores(
      max_cores = n_marks
    )
    fractions <- base::unlist(
      parallel::mclapply(
        mark_indices, .compute_overlap,
        mc.cores = n_cores
      )
    )
  } else {
    fractions <- base::vapply(mark_indices, .compute_overlap, base::numeric(1L))
  }

  # Warn for marks with zero annotations
  zero_marks <- histone[fractions == 0]
  for (zm in zero_marks) {
    zm_annots <- annot_by_type[[base::paste0(genome_prefix, zm)]]
    if (base::is.null(zm_annots) || base::length(zm_annots) == 0L) {
      fmt <- paste0("benchmark_enhancer_predictor: No annotations found ",
                    "for histone mark '%s'")
      base::warning(base::sprintf(fmt, zm))
    }
  }

  # Build and return sorted result data.frame
  histone_test <- base::data.frame(
    Histone_Mark = histone,
    Fraction_of_Overlap = base::as.character(fractions),
    stringsAsFactors = FALSE
  )
  histone_test <- histone_test[
    base::order(histone_test[, 2],
      decreasing = TRUE), ]

  base::return(histone_test)
}

#' Evaluate histone marks against any curated reference database
#'
#' Returns overlap fractions per histone mark against a reference database,
#' allowing comparison of which marks best predict regions in the reference.
#' Supports any database in database, not just FANTOM.
#'
#' @param database epiRomics class database
#'   containing all data initially loaded
#' @param histone name or vector of histone
#'   mark(s), must match name in epiROmics_dB@meta,
#'   default set to h3k4me1
#' @param curated_database database to test
#'   histone marks against, must match name in
#'   epiROmics_dB@meta default set to fantom
#' @return Variable of class dataframe further exploring top histone marks that
#'   may determine enhancer regions
#' @export
#' @examples
#' db <- make_example_database()
#' ## Use h3k27ac as the curated reference; compare h3k4me1 against it
#' test <- benchmark_enhancer_predictor(
#'   db,
#'   histone = "h3k4me1",
#'   curated_database = "h3k27ac"
#' )
#' test
benchmark_enhancer_predictor <- function(
    database,
    histone = "h3k4me1",
    curated_database = "fantom") {
  # Parameter validation
  validate_database(database, "benchmark_enhancer_predictor")
  validate_character_param(
    curated_database,
    "curated_database",
    "benchmark_enhancer_predictor"
  )

  # Validate histone parameter (can be single string or vector)
  if (!base::is.character(histone)) {
    base::stop(
      "benchmark_enhancer_predictor: ",
      "histone must be a ",
      "character vector"
    )
  }
  if (base::length(histone) == 0) {
    base::stop(
      "benchmark_enhancer_predictor: ",
      "histone cannot be empty"
    )
  }
  if (base::any(
    base::is.na(histone) |
      histone == ""
  )) {
    base::stop(
      "benchmark_enhancer_predictor: ",
      "histone cannot contain ",
      "NA or empty values"
    )
  }

  tryCatch(
    {
      # Check if histone marks exist in the database
      available_marks <- database@meta$name[
        database@meta$type %in%
          base::c("histone", "chip")
      ]
      missing_marks <- histone[
        !(histone %in% available_marks)
      ]
      if (base::length(missing_marks) > 0) {
        msg <- paste0(
          "benchmark_enhancer_predictor",
          ": Histone marks not found in ",
          "database: %s. Available marks: %s"
        )
        base::stop(base::sprintf(
          msg,
          base::paste(missing_marks, collapse = ", "),
          base::paste(available_marks, collapse = ", ")
        ))
      }

      # Check if curated database exists
      if (!(curated_database %in% available_marks)) {
        msg <- paste0(
          "benchmark_enhancer_predictor",
          ": Curated database '%s' not ",
          "found in database. ",
          "Available marks: %s"
        )
        base::stop(base::sprintf(
          msg,
          curated_database,
          base::paste(available_marks, collapse = ", ")
        ))
      }

      # Compute overlap fractions and build result via extracted helper
      histone_test <- .compute_overlap_fractions(
        database,
        histone,
        curated_database
      )

      base::return(histone_test)
    },
    error = function(e) {
      base::stop(base::sprintf(
        "benchmark_enhancer_predictor: %s",
        e$message
      ))
    }
  )
}
