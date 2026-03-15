#' Identifies putative enhanceosome regions by cross-referencing candidate enhancer regions against co-TF enrichment
#'
#' @param epiRomics_putative_enhancers epiRomics class database containing putative enhancer calls
#' @param epiRomics_dB epiRomics class database containing all data initially loaded
#' @return Variable of class epiRomics with enhanceosome annotations
#' @export
#' @examples
#' db <- epiRomicsS4(annotations = GenomicRanges::GRanges(),
#'   meta = data.frame(name = character(), type = character(),
#'     file = character(), stringsAsFactors = FALSE),
#'   genome = "hg38")
#' tryCatch(epiRomics_enhanceosome(db, db), error = function(e) message(e$message))
#' \donttest{
#' enhanceosome <- epiRomics_enhanceosome(epiRomics_putative_enhancers, epiRomics_dB)
#' }
epiRomics_enhanceosome <- function(epiRomics_putative_enhancers, epiRomics_dB) {
  epiRomics_chips <- epiRomics_dB@meta[epiRomics_dB@meta$type == "chip", "name"]
  epiRomics_chips_db_access <- base::paste0(epiRomics_dB@genome, "_custom_", epiRomics_chips)
  # Pre-split annotations by type to avoid repeated linear scans
  annot_by_type <- base::split(
    epiRomics_dB@annotations,
    epiRomics_dB@annotations$type
  )

  query_gr <- epiRomics_putative_enhancers@annotations
  n_chips <- base::length(epiRomics_chips)

  # Build subject list for lapply (avoids per-iteration hash lookup)
  subject_list <- base::lapply(epiRomics_chips_db_access, function(key) annot_by_type[[key]])

  # Parallel countOverlaps when >2 ChIP targets and parallel available
  use_par <- n_chips > 2L &&
    base::requireNamespace("parallel", quietly = TRUE) &&
    .Platform$OS.type == "unix"

  apply_fn <- if (use_par) {
    n_cores <- base::min(n_chips, base::max(1L, parallel::detectCores() - 1L))
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

  # Add counts directly as mcols (avoids slow GRanges->data.frame->GRanges round-trip)
  for (j in base::seq_along(epiRomics_chips)) {
    GenomicRanges::mcols(query_gr)[[epiRomics_chips[j]]] <- count_cols[[j]]
  }
  count_matrix <- base::do.call(base::cbind, count_cols)
  GenomicRanges::mcols(query_gr)$ChIP_Hits <- base::rowSums(count_matrix)

  # Sort by ChIP_Hits descending
  epiRomics_enhanceosome_data_table <- query_gr[
    base::order(GenomicRanges::mcols(query_gr)$ChIP_Hits, decreasing = TRUE)]
  txdb_obj <- resolve_txdb(epiRomics_dB@txdb)
  # Cache annotatePeak result for repeated calls with same input
  cache_dir <- base::file.path(base::tempdir(), "epiRomics_cache")
  if (!base::dir.exists(cache_dir)) base::dir.create(cache_dir, recursive = TRUE)
  annot_key <- digest::digest(base::paste(
    base::paste(base::as.character(GenomeInfoDb::seqnames(epiRomics_enhanceosome_data_table)), collapse = ","),
    base::paste(BiocGenerics::start(epiRomics_enhanceosome_data_table), collapse = ","),
    epiRomics_dB@txdb, epiRomics_dB@organism, sep = "|"
  ))
  annot_cache <- base::file.path(cache_dir, base::paste0("annotatePeak_", annot_key, ".rds"))

  if (base::file.exists(annot_cache)) {
    base::message("Loading cached annotatePeak results...")
    epiRomics_enhanceosome_data_table_annotated <- base::readRDS(annot_cache)
  } else {
    epiRomics_enhanceosome_data_table_annotated <- ChIPseeker::annotatePeak(epiRomics_enhanceosome_data_table,
      tssRegion = c(-3000, 3000), TxDb = txdb_obj, annoDb = epiRomics_dB@organism
    )
    epiRomics_enhanceosome_data_table_annotated <- epiRomics_enhanceosome_data_table_annotated@anno
    base::saveRDS(epiRomics_enhanceosome_data_table_annotated, annot_cache)
  }
  # Set sequential names (1, 2, ..., N) reflecting ChIP_Hits rank.

  # These names persist through downstream filtering (e.g.,
  # epiRomics_regions_of_interest) so that filtered subsets retain
  # their original enhanceosome rank as GRanges names / row names.
  base::names(epiRomics_enhanceosome_data_table_annotated) <- base::as.character(
    base::seq_len(base::length(epiRomics_enhanceosome_data_table_annotated))
  )
  epiRomics_enhanceosome <- epiRomics_putative_enhancers
  epiRomics_enhanceosome@annotations <- epiRomics_enhanceosome_data_table_annotated
  base::return(epiRomics_enhanceosome)
}
