#' Annotate enhanceosome peaks via ChIPseeker with RDS caching
#'
#' Sets up a cache directory, computes a digest key from the GRanges
#' coordinates and annotation parameters, and either loads a cached
#' annotatePeak result or runs ChIPseeker::annotatePeak fresh. Returns
#' the annotated GRanges with sequential names.
#'
#' @param epiRomics_enhanceosome_data_table GRanges
#'   of sorted enhanceosome peaks.
#' @param database epiRomics class database (for txdb, organism).
#' @param txdb_obj TxDb object for annotation.
#' @return Annotated GRanges with sequential names reflecting ChIP_Hits rank.
#' @noRd
.annotate_enhanceosome_peaks <- function(
    epiRomics_enhanceosome_data_table,
    database,
    txdb_obj) {
  # Cache annotatePeak result for repeated calls with same input
  cache_dir <- base::file.path(base::tempdir(), "epiRomics_cache")
  if (!base::dir.exists(cache_dir)) {
    base::dir.create(cache_dir, recursive = TRUE)
  }
  annot_key <- digest::digest(base::paste(
    base::paste(
      base::as.character(GenomeInfoDb::seqnames(
        epiRomics_enhanceosome_data_table
      )),
      collapse = ","
    ),
    base::paste(
      BiocGenerics::start(
        epiRomics_enhanceosome_data_table
      ),
      collapse = ","
    ),
    database@txdb, database@organism, sep = "|"
  ))
  annot_cache <- base::file.path(
    cache_dir,
    base::paste0("annotatePeak_", annot_key, ".rds")
  )

  if (base::file.exists(annot_cache)) {
    base::message("Loading cached annotatePeak results...")
    annotated <- base::readRDS(annot_cache)
  } else {
    annotated <- ChIPseeker::annotatePeak(
      epiRomics_enhanceosome_data_table,
      tssRegion = c(-3000, 3000),
      TxDb = txdb_obj,
      annoDb = database@organism
    )
    annotated <- annotated@anno
    base::saveRDS(annotated, annot_cache)
  }

  # Set sequential names (1, 2, ..., N) reflecting ChIP_Hits rank.
  # These names persist through downstream filtering so that filtered
  # subsets retain their original enhanceosome rank.
  base::names(annotated) <- base::as.character(
    base::seq_len(base::length(annotated))
  )

  base::return(annotated)
}

#' Identifies putative enhanceosome regions by
#' cross-referencing candidate enhancer regions
#' against co-TF enrichment
#'
#' @param putative_enhancers epiRomics class
#'   database containing putative enhancer calls
#' @param database epiRomics class database
#'   containing all data initially loaded
#' @return Variable of class epiRomics with enhanceosome annotations
#' @export
#' @examples
#' db <- make_example_database()
#' eso <- make_example_enhanceosome(db)
#' length(methods::slot(eso, "annotations"))
find_enhanceosomes <- function(putative_enhancers, database) {
  epiRomics_chips <- database@meta[database@meta$type == "chip", "name"]
  epiRomics_chips_db_access <- base::paste0(
    database@genome, "_custom_",
    epiRomics_chips
  )
  # Pre-split annotations by type to avoid repeated linear scans
  annot_by_type <- base::split(
    database@annotations,
    database@annotations$type
  )

  query_gr <- putative_enhancers@annotations
  n_chips <- base::length(epiRomics_chips)

  # Build subject list for lapply (avoids per-iteration hash lookup)
  subject_list <- base::lapply(
    epiRomics_chips_db_access,
    function(key) annot_by_type[[key]]
  )

  # Parallel countOverlaps when >2 ChIP targets and parallel available
  use_par <- n_chips > 2L &&
    base::requireNamespace("parallel", quietly = TRUE) &&
    .Platform$OS.type == "unix"

  apply_fn <- if (use_par) {
    n_cores <- .detect_cores(
      max_cores = n_chips
    )
    function(x, f) parallel::mclapply(x, f, mc.cores = n_cores)
  } else {
    base::lapply
  }

  count_cols <- apply_fn(subject_list, function(s_gr) {
    if (base::is.null(s_gr) || base::length(s_gr) == 0L) {
      return(base::rep(0L, base::length(query_gr)))
    }
    GenomicRanges::countOverlaps(query_gr, s_gr)
  })

  # Add counts directly as mcols
  # (avoids slow GRanges->df->GRanges round-trip)
  for (j in base::seq_along(epiRomics_chips)) {
    GenomicRanges::mcols(query_gr)[[epiRomics_chips[j]]] <- count_cols[[j]]
  }
  count_matrix <- base::do.call(base::cbind, count_cols)
  GenomicRanges::mcols(query_gr)$ChIP_Hits <- base::rowSums(count_matrix)

  # Sort by ChIP_Hits descending
  epiRomics_enhanceosome_data_table <- query_gr[
    base::order(GenomicRanges::mcols(query_gr)$ChIP_Hits, decreasing = TRUE)]
  txdb_obj <- resolve_txdb(database@txdb)

  # Annotate peaks via extracted helper
  annotated <- .annotate_enhanceosome_peaks(
    epiRomics_enhanceosome_data_table,
    database, txdb_obj
  )

  enhanceosome <- putative_enhancers
  enhanceosome@annotations <- annotated
  base::return(enhanceosome)
}
