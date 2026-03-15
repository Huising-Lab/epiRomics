#' Evaluate histone marks against any curated reference database
#'
#' Returns overlap fractions per histone mark against a reference database,
#' allowing comparison of which marks best predict regions in the reference.
#' Supports any database in epiRomics_dB, not just FANTOM.
#'
#' @param epiRomics_dB epiRomics class database containing all data initially loaded
#' @param epiRomics_histone name or vector of histone mark(s), must match name in
#'   epiROmics_dB@meta, default set to h3k4me1
#' @param epiRomics_curated_database database to test histone marks against, must
#'   match name in epiROmics_dB@meta default set to fantom
#' @return Variable of class dataframe further exploring top histone marks that
#'   may determine enhancer regions
#' @export
#' @examples
#' db <- epiRomicsS4(annotations = GenomicRanges::GRanges(),
#'   meta = data.frame(name = character(), type = character(),
#'     file = character(), stringsAsFactors = FALSE),
#'   genome = "hg38")
#' tryCatch(epiRomics_enhancer_predictor_to_ref(db), error = function(e) message(e$message))
#' \donttest{
#' test <- epiRomics_enhancer_predictor_to_ref(epiRomics_dB)
#' }
epiRomics_enhancer_predictor_to_ref <- function(epiRomics_dB, epiRomics_histone = "h3k4me1", epiRomics_curated_database = "fantom") {
  # Parameter validation
  validate_epiRomics_dB(epiRomics_dB, "epiRomics_enhancer_predictor_to_ref")
  validate_character_param(epiRomics_curated_database, "epiRomics_curated_database", "epiRomics_enhancer_predictor_to_ref")

  # Validate histone parameter (can be single string or vector)
  if (!base::is.character(epiRomics_histone)) {
    base::stop("epiRomics_enhancer_predictor_to_ref: epiRomics_histone must be a character vector")
  }
  if (base::length(epiRomics_histone) == 0) {
    base::stop("epiRomics_enhancer_predictor_to_ref: epiRomics_histone cannot be empty")
  }
  if (base::any(base::is.na(epiRomics_histone) | epiRomics_histone == "")) {
    base::stop("epiRomics_enhancer_predictor_to_ref: epiRomics_histone cannot contain NA or empty values")
  }

  tryCatch(
    {
      # Check if histone marks exist in the database
      available_marks <- epiRomics_dB@meta$name[epiRomics_dB@meta$type %in% base::c("histone", "chip")]
      missing_marks <- epiRomics_histone[!(epiRomics_histone %in% available_marks)]
      if (base::length(missing_marks) > 0) {
        base::stop(base::sprintf(
          "epiRomics_enhancer_predictor_to_ref: Histone marks not found in database: %s. Available marks: %s",
          base::paste(missing_marks, collapse = ", "),
          base::paste(available_marks, collapse = ", ")
        ))
      }

      # Check if curated database exists
      if (!(epiRomics_curated_database %in% available_marks)) {
        base::stop(base::sprintf(
          "epiRomics_enhancer_predictor_to_ref: Curated database '%s' not found in database. Available marks: %s",
          epiRomics_curated_database,
          base::paste(available_marks, collapse = ", ")
        ))
      }

      # Pre-allocate matrix with proper dimensions
      histone_test <- base::matrix(data = NA, nrow = base::length(epiRomics_histone), ncol = 2)
      base::colnames(histone_test) <- base::c("Histone_Mark", "Fraction_of_Overlap")
      histone_test[, 1] <- epiRomics_histone

      # Pre-split annotations by type for O(1) lookup
      annot_by_type <- base::split(
        epiRomics_dB@annotations,
        epiRomics_dB@annotations$type
      )

      # Calculate curated database count once
      mark2 <- base::paste0(epiRomics_dB@genome, "_custom_", epiRomics_curated_database)
      mark2_annotations <- annot_by_type[[mark2]]
      if (base::is.null(mark2_annotations)) mark2_annotations <- epiRomics_dB@annotations[0]
      mark2_count <- base::length(mark2_annotations)

      if (mark2_count == 0) {
        base::stop(base::sprintf("epiRomics_enhancer_predictor_to_ref: No annotations found for curated database '%s'", mark2))
      }

      # Pre-compute reduce(mark2) ONCE — this is constant across all marks.
      # For huge datasets (100k+ regions), reduce() is expensive (sort + merge).
      # Moving it out of the loop avoids redundant O(n log n) per iteration.
      reduced_mark2 <- GenomicRanges::reduce(mark2_annotations)

      # Pre-build all mark1 keys and fetch annotations in bulk
      genome_prefix <- base::paste0(epiRomics_dB@genome, "_custom_")
      mark1_keys <- base::paste0(genome_prefix, epiRomics_histone)

      # Compute overlap fractions — parallelise when marks > 2 and parallel available
      .compute_overlap <- function(idx) {
        mark1_annots <- annot_by_type[[mark1_keys[idx]]]
        if (base::is.null(mark1_annots) || base::length(mark1_annots) == 0L) {
          return(0)
        }
        # reduce() + intersect() are C-level GenomicRanges ops
        overlap_count <- base::length(
          GenomicRanges::intersect(
            GenomicRanges::reduce(mark1_annots),
            reduced_mark2
          )
        )
        overlap_count / mark2_count
      }

      n_marks <- base::length(epiRomics_histone)
      mark_indices <- base::seq_len(n_marks)

      if (n_marks > 2L && base::requireNamespace("parallel", quietly = TRUE)) {
        n_cores <- base::min(n_marks, base::max(1L, parallel::detectCores() - 1L))
        fractions <- base::unlist(parallel::mclapply(
          mark_indices, .compute_overlap, mc.cores = n_cores
        ))
      } else {
        fractions <- base::vapply(mark_indices, .compute_overlap, base::numeric(1L))
      }

      # Warn for marks with zero annotations
      zero_marks <- epiRomics_histone[fractions == 0]
      for (zm in zero_marks) {
        zm_annots <- annot_by_type[[base::paste0(genome_prefix, zm)]]
        if (base::is.null(zm_annots) || base::length(zm_annots) == 0L) {
          base::warning(base::sprintf(
            "epiRomics_enhancer_predictor_to_ref: No annotations found for histone mark '%s'", zm
          ))
        }
      }

      histone_test[, 2] <- fractions

      # Convert matrix to data frame and sort (use base::, not BiocGenerics:: — this is a matrix)
      histone_test <- base::as.data.frame(histone_test, stringsAsFactors = FALSE)
      histone_test <- histone_test[base::order(histone_test[, 2], decreasing = TRUE), ]

      base::return(histone_test)
    },
    error = function(e) {
      base::stop(base::sprintf("epiRomics_enhancer_predictor_to_ref: %s", e$message))
    }
  )
}

#' Evaluate histone marks against a curated database (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `epiRomics_enhancer_predictor_test()` was renamed to
#' [epiRomics_enhancer_predictor_to_ref()] for clarity.
#'
#' @inheritParams epiRomics_enhancer_predictor_to_ref
#' @return data.frame with overlap fractions
#' @export
#' @examples
#' db <- epiRomicsS4(annotations = GenomicRanges::GRanges(),
#'   meta = data.frame(name = character(), type = character(),
#'     file = character(), stringsAsFactors = FALSE),
#'   genome = "hg38")
#' tryCatch(epiRomics_enhancer_predictor_test(db),
#'   error = function(e) message(e$message))
#' \donttest{
#' # Use epiRomics_enhancer_predictor_to_ref() instead
#' result <- epiRomics_enhancer_predictor_to_ref(epiRomics_dB)
#' }
epiRomics_enhancer_predictor_test <- function(epiRomics_dB, epiRomics_histone = "h3k4me1", epiRomics_curated_database = "fantom") {
  .Deprecated("epiRomics_enhancer_predictor_to_ref",
    msg = "epiRomics_enhancer_predictor_test() is deprecated. Use epiRomics_enhancer_predictor_to_ref() instead."
  )
  epiRomics_enhancer_predictor_to_ref(epiRomics_dB, epiRomics_histone, epiRomics_curated_database)
}
