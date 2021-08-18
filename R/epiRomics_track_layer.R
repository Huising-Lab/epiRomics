#' Visualizes data from epiRomics results
#'
#' @param epiRomics_putative_enhanceosome epiRomics class database containing putative enhanceosome calls
#' @param epiRomics_index numeric of row value from epiRomics_putative_enhanceosome to visualize
#' @param epiRomics_dB epiRomics class database containing all data initially loaded
#' @param epiRomics_track_connection data frame containing bigwig track locations and their names
#' @param epiRomics_keep_epitracks logical indicating whether to show enhancer and chip tracks, default is TRUE
#' @return GViz plot
#' @export
epiRomics_track_layer <- function(epiRomics_putative_enhanceosome, epiRomics_index, epiRomics_dB, epiRomics_track_connection, epiRomics_keep_epitracks = TRUE) {
  epiRomics_class_save <- epiRomics_putative_enhanceosome
  epiRomics_putative_enhanceosome <- epiRomics_putative_enhanceosome@annotations[epiRomics_index, ]
  Gviz::GeneRegionTrack(eval(parse(text = epiRomics_class_save@txdb)))
  chr <- base::as.character(BiocGenerics::unique(GenomeInfoDb::seqnames(epiRomics_putative_enhanceosome)))
  gtrack <- Gviz::GenomeAxisTrack()
  txTr <- Gviz::GeneRegionTrack(base::eval(base::parse(text = epiRomics_class_save@txdb)),
    start = epiRomics_putative_enhanceosome$geneStart,
    end = epiRomics_putative_enhanceosome$geneEnd, chromosome = epiRomics_putative_enhanceosome$geneChr,
    symbol = epiRomics_putative_enhanceosome$SYMBOL, name = epiRomics_putative_enhanceosome$SYMBOL, transcriptAnnotation = "symbol",
    geneSymbol = TRUE, howId = TRUE, min.height = 10, size = 10
  )
  atrack <- Gviz::AnnotationTrack(epiRomics_putative_enhanceosome, name = base::paste0(
    "Putative_Enhancer_Index_",
    epiRomics_index
  ), feature = "Enhancer", min.height = 15, size = 15)
  itrack <- Gviz::IdeogramTrack(genome = "mm10", chromosome = chr)
  if (base::min(atrack) < base::min(txTr)) {
    min_track <- base::min(atrack)
  } else {
    min_track <- base::min(txTr)
  }
  if (base::max(atrack) < base::max(txTr)) {
    max_track <- base::max(txTr)
  } else {
    max_track <- base::max(atrack)
  }
  min_track <- min_track - 1000
  max_track <- max_track + 1000
  region_track <- GenomicRanges::GRanges(seqnames = chr, ranges = IRanges::IRanges(start = min_track, end = max_track))
  epiRomics_db_data_table <- epiRomics_class_save@meta
  epiRomics_chips <- epiRomics_db_data_table[epiRomics_db_data_table$type == "chip", "name"]
  epiRomics_chips_db_access <- base::paste0("mm10_custom_", epiRomics_chips)
  chip_tracks <- list()
  for (i in 1:base::length(epiRomics_chips_db_access)) {
    chip_tracks <- base::c(chip_tracks, Gviz::AnnotationTrack(IRanges::subsetByOverlaps(epiRomics_dB@annotations[epiRomics_dB@annotations$type ==
      epiRomics_chips_db_access[i], ], region_track),
    name = epiRomics_chips[i],
    feature = "ChIP", min.height = 10, size = 10
    ))
  }
  epiRomics_track_connection_chromatin <- epiRomics_track_connection[epiRomics_track_connection$type == "atac", ]
  range_max <- maxCovFiles(epiRomics_track_connection_chromatin[, 1], region_track)
  max_cov <- plyr::round_any((range_max$X + (range_max$X * 0.1)), accuracy = 10, f = ceiling)
  chromatin_tracks <- list()
  for (i in 1:base::dim(epiRomics_track_connection_chromatin)[1]) {
    chromatin_tracks <- c(chromatin_tracks, Gviz::DataTrack(size = 10, epiRomics_track_connection_chromatin[i, 1], col = epiRomics_track_connection_chromatin[i, 3], fill = epiRomics_track_connection_chromatin[i, 3], col.histogram = epiRomics_track_connection_chromatin[i, 3], ylim = c(0, max_cov), start = min_track, end = max_track, genome = epiRomics_dB@genome, type = "hist", chromosome = chr, name = epiRomics_track_connection_chromatin[i, 2]))
  }

  if (table(epiRomics_track_connection$type == "rna")[[TRUE]] > 0) {
    rna_tracks <- list()
    epiRomics_track_connection_rna <- epiRomics_track_connection[epiRomics_track_connection$type == "rna", ]
    range_max <- maxCovFiles(epiRomics_track_connection_rna[, 1], region_track)
    max_cov <- plyr::round_any((range_max$X + (range_max$X * 0.1)), accuracy = 10, f = ceiling)
    for (i in 1:base::dim(epiRomics_track_connection_rna)[1]) {
      rna_tracks <- c(rna_tracks, Gviz::DataTrack(size = 10, epiRomics_track_connection_rna[i, 1], col = epiRomics_track_connection_rna[i, 3], fill = epiRomics_track_connection_rna[i, 3], col.histogram = epiRomics_track_connection_rna[i, 3], ylim = c(0, max_cov), start = min_track, end = max_track, genome = epiRomics_dB@genome, type = "hist", chromosome = chr, name = epiRomics_track_connection_rna[i, 2]))
    }
  }
  if (epiRomics_keep_epitracks == TRUE) {
    Gviz::plotTracks(base::c(itrack, gtrack, txTr, chromatin_tracks, rna_tracks, atrack, chip_tracks), from = min_track, to = max_track)
  }
  if (epiRomics_keep_epitracks == FALSE) {
    Gviz::plotTracks(base::c(itrack, gtrack, txTr, chromatin_tracks, rna_tracks), from = min_track, to = max_track)
  }
}
