#' Gviz-based gene-centered chromatin visualization (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This Gviz-based function has been superseded by
#' \code{\link{epiRomics_track_layer_gene}}, which uses base R graphics
#' for ~20x faster rendering (~2 seconds vs ~40+ seconds).
#'
#' Generates a multi-track genomic visualization centered on a gene, showing
#' ATAC-seq and RNA-seq signal tracks. Genome is automatically detected from
#' the database object.
#'
#' @param epiRomics_gene_name character of name of gene to visualize
#' @param epiRomics_dB epiRomics class database containing all data initially loaded
#' @param epiRomics_track_connection data frame containing bigwig track locations and their names
#' @return GViz plot
#' @export
#' @examples
#' db <- epiRomicsS4(annotations = GenomicRanges::GRanges(),
#'   meta = data.frame(name = character(), type = character(),
#'     file = character(), stringsAsFactors = FALSE),
#'   genome = "hg38")
#' tryCatch(epiRomics_chromatin_tracks("INS", db,
#'   data.frame(file = character(), name = character(),
#'     color = character(), type = character())),
#'   error = function(e) message(e$message))
#' \donttest{
#' # Use epiRomics_track_layer_gene() instead:
#' epiRomics_track_layer_gene("Ins1", epiRomics_dB, epiRomics_track_connection)
#' }
epiRomics_chromatin_tracks <-
  function(epiRomics_gene_name,
           epiRomics_dB,
           epiRomics_track_connection) {
    if (!base::requireNamespace("Gviz", quietly = TRUE)) {
      base::stop("Package 'Gviz' is required for this deprecated function. ",
                 "Install it with: BiocManager::install('Gviz')\n",
                 "Or use epiRomics_track_layer_gene() instead (no Gviz dependency).",
                 call. = FALSE)
    }
    .Deprecated("epiRomics_track_layer_gene",
      msg = "epiRomics_chromatin_tracks() is deprecated (slow Gviz rendering). Use epiRomics_track_layer_gene() for ~20x faster base R graphics."
    )
    organism_db <- base::getExportedValue(epiRomics_dB@organism, epiRomics_dB@organism)
    txdb_obj <- resolve_txdb(epiRomics_dB@txdb)
    genome <- epiRomics_dB@genome

    epiRomics_entrez_id <-
      base::as.vector(AnnotationDbi::mapIds(
        organism_db,
        epiRomics_gene_name,
        "ENTREZID",
        "SYMBOL"
      ))
    epiRomics_gene_map <- GenomicFeatures::genes(txdb_obj)
    epiRomics_gene_map_track <-
      epiRomics_gene_map[epiRomics_gene_map$gene_id == epiRomics_entrez_id, ]
    # Resize to fix max width issue for certain genes (e.g., MafA)
    epiRomics_gene_map_track <-
      GenomicRanges::resize(
        epiRomics_gene_map_track,
        width = BiocGenerics::width(epiRomics_gene_map_track) + (4000 * 2),
        fix = "center"
      )
    chr <-
      base::as.character(BiocGenerics::unique(GenomeInfoDb::seqnames(epiRomics_gene_map_track)))
    gtrack <- Gviz::GenomeAxisTrack()
    txTr <-
      Gviz::GeneRegionTrack(
        txdb_obj,
        start = BiocGenerics::start(epiRomics_gene_map_track),
        end = BiocGenerics::end(epiRomics_gene_map_track),
        chromosome = chr,
        symbol = epiRomics_gene_name,
        name = epiRomics_gene_name,
        transcriptAnnotation = "symbol",
        geneSymbol = TRUE,
        showId = TRUE
      )
    itrack <- get_ideogram_track(genome, chr)
    # Pre-import ATAC BigWig data and build tracks (single-pass I/O)
    epiRomics_track_connection_chromatin <-
      epiRomics_track_connection[epiRomics_track_connection$type == "atac", ]
    atac_result <- .build_data_tracks(
      bw_paths = epiRomics_track_connection_chromatin[, 1],
      region = epiRomics_gene_map_track,
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
        region = epiRomics_gene_map_track,
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
        region = epiRomics_gene_map_track,
        sample_names = tc_histone[, 2],
        colors = tc_histone[, 3],
        genome = genome,
        chr = chr
      )
      histone_signal_tracks <- histone_result$tracks
    }
    Gviz::plotTracks(
      base::c(itrack, gtrack, txTr, chromatin_tracks, rna_tracks,
              histone_signal_tracks),
      extend.right = 1000,
      extend.left = 1000
    )
  }
