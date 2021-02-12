#' Visualizes chromatin availability
#'
#' @param epiRomics_gene_name character of name of gene to visualize
#' @param epiRomics_dB epiRomics class database containing all data initially loaded
#' @param epiRomics_track_connection data frame containing bigwig track locations and their names
#' @return GViz plot
#' @export
epiRomics_chromatin_tracks <- function(epiRomics_gene_name, epiRomics_dB, epiRomics_track_connection) {
    Gviz::GeneRegionTrack(eval(parse(text = epiRomics_dB@txdb)))
    epiRomics_entrez_id <- base::as.vector(AnnotationDbi::mapIds(base::eval(base::parse(text=epiRomics_dB@organism)), epiRomics_gene_name, 'ENTREZID', 'SYMBOL'))
    epiRomics_gene_map <- GenomicFeatures::genes((base::eval(base::parse(text = epiRomics_dB@txdb))))
    epiRomics_gene_map_track <- epiRomics_gene_map[epiRomics_gene_map$gene_id==epiRomics_entrez_id, ]
    epiRomics_gene_map_track <- resize(epiRomics_gene_map_track, width = width(epiRomics_gene_map_track)+(4000*2), fix = "center")
    chr <- base::as.character(BiocGenerics::unique(GenomeInfoDb::seqnames(epiRomics_gene_map_track)))
    gtrack <- Gviz::GenomeAxisTrack()
    txTr <- Gviz::GeneRegionTrack(base::eval(base::parse(text = epiRomics_dB@txdb)), start = as.data.frame(epiRomics_gene_map_track)[,"start"], end=as.data.frame(epiRomics_gene_map_track)[,"end"], chromosome = chr,
        symbol = epiRomics_gene_name, name = epiRomics_gene_name, transcriptAnnotation = "symbol",
        geneSymbol = TRUE, showId = TRUE)
    itrack <- Gviz::IdeogramTrack(genome = "mm10", chromosome = chr)
    range_max <- maxCovFiles(epiRomics_track_connection[,1], epiRomics_gene_map_track)
    max_cov <- range_max$X+5
    chromatin_tracks <- list()
    for (i in 1:base::dim(epiRomics_track_connection)[1]) {
        chromatin_tracks <- c(chromatin_tracks, Gviz::DataTrack(size = 10, epiRomics_track_connection[i,1], col=epiRomics_track_connection[i,3], fill=epiRomics_track_connection[i,3], col.histogram=epiRomics_track_connection[i,3], ylim=c(0, max_cov),range = epiRomics_track_connection[i,1], genome = epiRomics_dB@genome, type = "hist",
                                                                chromosome = chr, strand=epiRomics_gene_map_track@strand, name = epiRomics_track_connection[i,2]))
    }


    Gviz::plotTracks(base::c(itrack, gtrack, txTr, chromatin_tracks), extend.right = 1000, extend.left = 1000)
}

