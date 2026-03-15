#' Build epiRomics database
#'
#' Reads epigenomic annotation files from a CSV manifest and builds a unified
#' epiRomics database for downstream analysis. Supports optional extra columns
#' for ChIP/histone peak files (signal, pval, qval, peak).
#'
#' @param epiRomics_db_file character string of path to properly formatted csv file
#'   containing epigenetic data. [See vignette for more details]
#' @param txdb_organism a character string containing the TxDB associated with your data.
#' @param epiRomics_genome a character string containing the genome associated with your data. e.g. "mm10" or "hg38".
#' @param epiRomics_organism a character string containing the org.db associated with your data.
#' @param extraCols named character vector of extra columns to read from chip/histone
#'   BED files. Default is NULL (no extra columns). Set to
#'   \code{c(signal = "numeric", pval = "numeric", qval = "numeric", peak = "numeric")}
#'   to read narrowPeak columns.
#' @param data_dir optional character string specifying the root directory for
#'   resolving relative file paths in the CSV manifest. When provided, any
#'   relative path in the \code{path} column is prefixed with \code{data_dir}.
#'   This is especially useful with cached data from
#'   \code{\link{epiRomics_cache_data}} where the CSV uses relative paths.
#'   Default is NULL (paths used as-is).
#' @return Variable of class epiRomics for further downstream analysis
#' @export
#' @importFrom methods new
#' @importFrom data.table %like%
#' @examples
#' tryCatch(epiRomics_build_dB("nonexistent.csv"), error = function(e) message(e$message))
#' \donttest{
#' epiRomics_dB <- epiRomics_build_dB(epiRomics_data_path, epiRomics_meta, "mm10")
#' }
epiRomics_build_dB <-
  function(epiRomics_db_file,
           txdb_organism,
           epiRomics_genome,
           epiRomics_organism,
           extraCols = NULL,
           data_dir = NULL) {
    validate_file_paths(epiRomics_db_file, "epiRomics_build_dB")
    validate_character_param(txdb_organism, "txdb_organism", "epiRomics_build_dB")
    validate_character_param(epiRomics_genome, "epiRomics_genome", "epiRomics_build_dB")
    validate_character_param(epiRomics_organism, "epiRomics_organism", "epiRomics_build_dB")

    tryCatch(
      {
        # --- Fast CSV read with data.table::fread() ---
        # fread() is 5-10x faster than utils::read.csv() for larger manifests
        # and handles column types more efficiently. Falls back to read.csv()
        # if data.table is unavailable (it's already in Imports via %like%).
        if (base::requireNamespace("data.table", quietly = TRUE)) {
          epiRomics_db_data_table <- base::as.data.frame(
            data.table::fread(epiRomics_db_file, stringsAsFactors = FALSE),
            stringsAsFactors = FALSE
          )
        } else {
          epiRomics_db_data_table <- utils::read.csv(epiRomics_db_file,
                                                      stringsAsFactors = FALSE)
        }

        required_cols <- base::c("path", "genome", "name", "format")
        missing_cols <- required_cols[!(required_cols %in% base::colnames(epiRomics_db_data_table))]
        if (base::length(missing_cols) > 0) {
          base::stop(base::sprintf(
            "Missing required columns in CSV file: %s",
            base::paste(missing_cols, collapse = ", ")
          ))
        }

        n_rows <- base::nrow(epiRomics_db_data_table)
        if (n_rows == 0L) {
          base::stop("CSV file is empty")
        }

        # --- Resolve relative paths against data_dir if provided ---
        if (!base::is.null(data_dir)) {
          validate_character_param(data_dir, "data_dir", "epiRomics_build_dB")
          if (!base::dir.exists(data_dir)) {
            base::stop(base::sprintf(
              "data_dir does not exist: %s", data_dir
            ))
          }
          # Only prefix paths that are not already absolute
          is_relative <- !base::grepl("^(/|[A-Za-z]:)", epiRomics_db_data_table$path)
          epiRomics_db_data_table$path[is_relative] <- base::file.path(
            data_dir, epiRomics_db_data_table$path[is_relative]
          )
        }

        validate_file_paths(epiRomics_db_data_table$path, "epiRomics_build_dB")

        if (!base::all(epiRomics_db_data_table$genome == epiRomics_genome)) {
          base::stop(base::sprintf(
            "Genome mismatch. Expected '%s' but found: %s",
            epiRomics_genome, base::paste(base::unique(epiRomics_db_data_table$genome), collapse = ", ")
          ))
        }

        # --- Vectorized validation (replaces per-row loop) ---
        name_vals <- epiRomics_db_data_table$name
        format_vals <- epiRomics_db_data_table$format
        bad_name <- base::which(base::is.na(name_vals) | name_vals == "")
        bad_format <- base::which(base::is.na(format_vals) | format_vals == "")
        if (base::length(bad_name) > 0L) {
          base::stop(base::sprintf("Empty or NA name in row %d", bad_name[1L]))
        }
        if (base::length(bad_format) > 0L) {
          base::stop(base::sprintf("Empty or NA format in row %d", bad_format[1L]))
        }

        # --- Read annotations (sequential — annotatr_cache is shared state) ---
        # Pre-compute which rows need extraCols to avoid per-row conditional
        needs_extra <- !base::is.null(extraCols) &
          format_vals %in% base::c("chip", "histone")

        for (i in base::seq_len(n_rows)) {
          row_extra <- if (needs_extra[i]) extraCols else NULL
          annotatr::read_annotations(
            con = epiRomics_db_data_table$path[i],
            genome = epiRomics_db_data_table$genome[i],
            name = name_vals[i],
            format = format_vals[i],
            extraCols = row_extra
          )
        }

        epiRomics_db_annot_list <- base::c(
          annotatr::builtin_annotations()[annotatr::builtin_annotations() %like% epiRomics_genome],
          (annotatr::annotatr_cache$list_env())
        )

        if (base::length(epiRomics_db_annot_list) == 0L) {
          base::stop(base::sprintf("No annotations found for genome '%s'", epiRomics_genome))
        }

        epiRomics_dB <- methods::new("epiRomicsS4")

        # --- Cache build_annotations result to RDS for fast repeated loads ---
        # Key: hash of genome + annotation list + input file modification times
        cache_dir <- base::file.path(base::tempdir(), "epiRomics_cache")
        if (!base::dir.exists(cache_dir)) base::dir.create(cache_dir, recursive = TRUE)
        file_mtimes <- base::paste(
          base::file.info(epiRomics_db_data_table$path)$mtime,
          collapse = "|"
        )
        cache_key <- digest::digest(base::paste(
          epiRomics_genome,
          base::paste(epiRomics_db_annot_list, collapse = ","),
          file_mtimes,
          sep = "|"
        ))
        cache_file <- base::file.path(cache_dir, base::paste0("annotations_", cache_key, ".rds"))

        if (base::file.exists(cache_file)) {
          base::message("Loading cached annotations...")
          epiRomics_dB@annotations <- base::readRDS(cache_file)
        } else {
          epiRomics_dB@annotations <-
            annotatr::build_annotations(genome = epiRomics_genome, annotations = epiRomics_db_annot_list)
          base::saveRDS(epiRomics_dB@annotations, cache_file)
        }

        if (base::length(epiRomics_dB@annotations) == 0L) {
          base::stop("Failed to build annotations - no data returned")
        }

        epiRomics_dB@meta <- epiRomics_db_data_table
        epiRomics_dB@txdb <- txdb_organism
        epiRomics_dB@organism <- epiRomics_organism
        epiRomics_dB@genome <- epiRomics_genome
        base::return(epiRomics_dB)
      },
      error = function(e) {
        base::stop(base::sprintf("epiRomics_build_dB: %s", e$message))
      }
    )
  }

#' Build epiRomics database (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' Use \code{\link{epiRomics_build_dB}} with the \code{extraCols} parameter instead.
#'
#' @inheritParams epiRomics_build_dB
#' @return Variable of class epiRomics for further downstream analysis
#' @export
#' @examples
#' tryCatch(epiRomics_build_dB_2("nonexistent.csv"), error = function(e) message(e$message))
#' \donttest{
#' epiRomics_dB <- epiRomics_build_dB_2(epiRomics_data_path, epiRomics_meta, "mm10")
#' }
epiRomics_build_dB_2 <-
  function(epiRomics_db_file,
           txdb_organism,
           epiRomics_genome,
           epiRomics_organism,
           data_dir = NULL) {
    .Deprecated("epiRomics_build_dB",
      msg = "epiRomics_build_dB_2() is deprecated. Use epiRomics_build_dB() with extraCols parameter instead."
    )
    epiRomics_build_dB(
      epiRomics_db_file, txdb_organism, epiRomics_genome, epiRomics_organism,
      extraCols = base::c(signal = "numeric", pval = "numeric", qval = "numeric", peak = "numeric"),
      data_dir = data_dir
    )
  }
