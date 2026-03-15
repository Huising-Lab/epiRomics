#' Gviz-based custom region visualization (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This Gviz-based function has been superseded by
#' \code{\link{epiRomics_track_layer}}, which uses base R graphics
#' for ~20x faster rendering (~2 seconds vs ~40+ seconds).
#' For gene-centered views, use \code{\link{epiRomics_track_layer_gene}}
#' or \code{\link{epiRomics_quick_view}}.
#'
#' Generates a multi-track genomic visualization for a custom genomic region,
#' showing gene models, ATAC-seq and RNA-seq signal tracks. Genome is
#' automatically detected from the database object.
#'
#' @param epiRomics_region GRanges of region to visualize
#' @param epiRomics_dB epiRomics class database containing all data initially loaded
#' @param epiRomics_track_connection data frame containing bigwig track locations and their names
#' @return Gviz plot
#' @export
#' @examples
#' db <- epiRomicsS4(annotations = GenomicRanges::GRanges(),
#'   meta = data.frame(name = character(), type = character(),
#'     file = character(), stringsAsFactors = FALSE),
#'   genome = "hg38")
#' tryCatch(epiRomics_region_tracks(GenomicRanges::GRanges(), db,
#'   data.frame(file = character(), name = character(),
#'     color = character(), type = character())),
#'   error = function(e) message(e$message))
#' \donttest{
#' # Use epiRomics_track_layer() or epiRomics_quick_view() instead:
#' epiRomics_quick_view(region = list(chr = "chr11", start = 2159779,
#'   end = 2161209), bw_paths = c(Signal = "track.bw"))
#' }
epiRomics_region_tracks <-
  function(epiRomics_region,
           epiRomics_dB,
           epiRomics_track_connection) {
    if (!base::requireNamespace("Gviz", quietly = TRUE)) {
      base::stop("Package 'Gviz' is required for this deprecated function. ",
                 "Install it with: BiocManager::install('Gviz')\n",
                 "Or use epiRomics_track_layer() or epiRomics_quick_view() instead (no Gviz dependency).",
                 call. = FALSE)
    }
    .Deprecated("epiRomics_track_layer",
      msg = "epiRomics_region_tracks() is deprecated (slow Gviz rendering). Use epiRomics_track_layer() or epiRomics_quick_view() for ~20x faster base R graphics."
    )
    txdb_obj <- resolve_txdb(epiRomics_dB@txdb)
    genome <- epiRomics_dB@genome

    chr <-
      base::as.character(BiocGenerics::unique(GenomeInfoDb::seqnames(epiRomics_region)))
    gtrack <- Gviz::GenomeAxisTrack()
    txTr <-
      Gviz::GeneRegionTrack(
        txdb_obj,
        start = BiocGenerics::start(epiRomics_region),
        end = BiocGenerics::end(epiRomics_region),
        chromosome = chr,
        collapseTranscripts = TRUE,
        symbol = epiRomics_region$SYMBOL,
        name = epiRomics_region$SYMBOL,
        transcriptAnnotation = "symbol",
        geneSymbol = TRUE,
        howId = TRUE,
        min.height = 10,
        size = 10
      )
    itrack <- get_ideogram_track(genome, chr)
    # Pre-import ATAC BigWig data and build tracks (single-pass I/O)
    epiRomics_track_connection_chromatin <-
      epiRomics_track_connection[epiRomics_track_connection$type == "atac", ]
    atac_result <- .build_data_tracks(
      bw_paths = epiRomics_track_connection_chromatin[, 1],
      region = epiRomics_region,
      sample_names = epiRomics_track_connection_chromatin[, 2],
      colors = epiRomics_track_connection_chromatin[, 3],
      genome = genome,
      chr = chr
    )
    chromatin_tracks <- atac_result$tracks

    # Pre-import RNA BigWig data and build tracks (single-pass I/O)
    rna_tracks <- base::list()
    if (base::any(epiRomics_track_connection$type == "rna")) {
      epiRomics_track_connection_rna <-
        epiRomics_track_connection[epiRomics_track_connection$type == "rna", ]
      rna_result <- .build_data_tracks(
        bw_paths = epiRomics_track_connection_rna[, 1],
        region = epiRomics_region,
        sample_names = epiRomics_track_connection_rna[, 2],
        colors = epiRomics_track_connection_rna[, 3],
        genome = genome,
        chr = chr
      )
      rna_tracks <- rna_result$tracks
    }
    # Pre-import histone BigWig data and build tracks (single-pass I/O)
    histone_signal_tracks <- base::list()
    if (base::any(epiRomics_track_connection$type == "histone")) {
      tc_histone <- epiRomics_track_connection[
        epiRomics_track_connection$type == "histone", , drop = FALSE]
      histone_result <- .build_data_tracks(
        bw_paths = tc_histone[, 1],
        region = epiRomics_region,
        sample_names = tc_histone[, 2],
        colors = tc_histone[, 3],
        genome = genome,
        chr = chr
      )
      histone_signal_tracks <- histone_result$tracks
    }
    Gviz::plotTracks(base::c(itrack, gtrack, txTr, chromatin_tracks, rna_tracks,
                              histone_signal_tracks))
  }
