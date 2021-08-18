#' Visualizes chromatin availability
#'
#' @param epiRomics_gene_name character of name of gene to visualize
#' @param epiRomics_dB epiRomics class database containing all data initially loaded
#' @param epiRomics_track_connection data frame containing bigwig track locations and their names
#' @return GViz plot
#' @export
epiRomics_chromatin_tracks <- function(epiRomics_gene_name, epiRomics_dB, epiRomics_track_connection) {
  Gviz::GeneRegionTrack(eval(parse(text = epiRomics_dB@txdb)))
  epiRomics_entrez_id <- base::as.vector(AnnotationDbi::mapIds(base::eval(base::parse(text = epiRomics_dB@organism)), epiRomics_gene_name, "ENTREZID", "SYMBOL"))
  epiRomics_gene_map <- GenomicFeatures::genes((base::eval(base::parse(text = epiRomics_dB@txdb))))
  epiRomics_gene_map_track <- epiRomics_gene_map[epiRomics_gene_map$gene_id == epiRomics_entrez_id, ]
  # Uncommented the below line to fix max width issue for certain genes, such as MafA
  epiRomics_gene_map_track <- resize(epiRomics_gene_map_track, width = width(epiRomics_gene_map_track) + (4000 * 2), fix = "center")
  chr <- base::as.character(BiocGenerics::unique(GenomeInfoDb::seqnames(epiRomics_gene_map_track)))
  gtrack <- Gviz::GenomeAxisTrack()
  txTr <- Gviz::GeneRegionTrack(base::eval(base::parse(text = epiRomics_dB@txdb)),
    start = as.data.frame(epiRomics_gene_map_track)[, "start"], end = as.data.frame(epiRomics_gene_map_track)[, "end"], chromosome = chr,
    symbol = epiRomics_gene_name, name = epiRomics_gene_name, transcriptAnnotation = "symbol",
    geneSymbol = TRUE, showId = TRUE
  )
  itrack <- Gviz::IdeogramTrack(genome = "mm10", chromosome = chr)
  epiRomics_track_connection_chromatin <- epiRomics_track_connection[epiRomics_track_connection$type == "atac", ]
  range_max <- maxCovFiles(epiRomics_track_connection_chromatin[, 1], epiRomics_gene_map_track)
  max_cov <- plyr::round_any((range_max$X + (range_max$X * 0.1)), accuracy = 10, f = ceiling)
  chromatin_tracks <- list()
  for (i in 1:base::dim(epiRomics_track_connection_chromatin)[1]) {
    chromatin_tracks <- c(chromatin_tracks, Gviz::DataTrack(
      size = 10, epiRomics_track_connection_chromatin[i, 1], col = epiRomics_track_connection_chromatin[i, 3], fill = epiRomics_track_connection_chromatin[i, 3], col.histogram = epiRomics_track_connection_chromatin[i, 3], ylim = c(0, max_cov), range = epiRomics_track_connection_chromatin[i, 1], genome = epiRomics_dB@genome, type = "hist",
      chromosome = chr, strand = epiRomics_gene_map_track@strand, name = epiRomics_track_connection_chromatin[i, 2]
    ))
  }
  if (table(epiRomics_track_connection$type == "rna")[[TRUE]] > 0) {
    rna_tracks <- list()
    epiRomics_track_connection_rna <- epiRomics_track_connection[epiRomics_track_connection$type == "rna", ]
    range_max <- maxCovFiles(epiRomics_track_connection_rna[, 1], epiRomics_gene_map_track)
    max_cov <- plyr::round_any((range_max$X + (range_max$X * 0.1)), accuracy = 10, f = ceiling)
    for (i in 1:base::dim(epiRomics_track_connection_rna)[1]) {
      rna_tracks <- c(rna_tracks, Gviz::DataTrack(
        size = 10, epiRomics_track_connection_rna[i, 1], col = epiRomics_track_connection_rna[i, 3], fill = epiRomics_track_connection_rna[i, 3], col.histogram = epiRomics_track_connection_rna[i, 3], ylim = c(0, max_cov), range = epiRomics_track_connection_rna[i, 1], genome = epiRomics_dB@genome, type = "hist",
        chromosome = chr, strand = epiRomics_gene_map_track@strand, name = epiRomics_track_connection_rna[i, 2]
      ))
    }
  }
  Gviz::plotTracks(base::c(itrack, gtrack, txTr, chromatin_tracks, rna_tracks), extend.right = 1000, extend.left = 1000)
}
