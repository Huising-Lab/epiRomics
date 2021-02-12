#' Identifies putative enhanceosome regions by cross-referencing candidate enhancer regions against co-TF enrichment
#'
#' @param epiRomics_putative_enhancers epiRomics class database containing putative enhancer calls
#' @param epiRomics_dB epiRomics class database containing all data initially loaded
#' @return Variable of class epiRomics further exploring candidate enhanceosome regions using co-ChIP hits
#' @export
epiRomics_enhanceosome <- function(epiRomics_putative_enhancers, epiRomics_dB) {
    epiRomics_chips <- epiRomics_dB@meta[epiRomics_dB@meta$type == "chip", "name"]
    epiRomics_chips_db_access <- base::paste0(epiRomics_dB@genome, "_custom_", epiRomics_chips)
    epiRomics_enhanceosome_data_table <- base::matrix(ncol = base::length(epiRomics_chips), nrow = base::length(epiRomics_putative_enhancers@annotations))
    base::colnames(epiRomics_enhanceosome_data_table) <- epiRomics_chips
    for (c in 1:base::length(epiRomics_chips)) {
        epiRomics_enhanceosome_data_table[, c] <- GenomicRanges::countOverlaps(epiRomics_putative_enhancers@annotations, 
            epiRomics_dB@annotations[epiRomics_dB@annotations$type == epiRomics_chips_db_access[c], ])
    }
    epiRomics_enhanceosome_data_table <- base::cbind(epiRomics_enhanceosome_data_table, ChIP_Hits = base::rowSums(epiRomics_enhanceosome_data_table), 
        base::as.data.frame(epiRomics_putative_enhancers@annotations))
    
    epiRomics_enhanceosome_data_table <- epiRomics_enhanceosome_data_table[base::order(epiRomics_enhanceosome_data_table$ChIP_Hits, 
        decreasing = TRUE), ]
    epiRomics_enhanceosome_data_table <- GenomicRanges::GRanges(epiRomics_enhanceosome_data_table)
    epiRomics_enhanceosome_data_table_annotated <- ChIPseeker::annotatePeak(epiRomics_enhanceosome_data_table, 
        tssRegion = c(-3000, 3000), TxDb = base::eval(base::parse(text = epiRomics_dB@txdb)), annoDb = epiRomics_dB@organism)
    epiRomics_enhanceosome_data_table_annotated <- epiRomics_enhanceosome_data_table_annotated@anno
    epiRomics_enhanceosome <- epiRomics_putative_enhancers
    epiRomics_enhanceosome@annotations <- epiRomics_enhanceosome_data_table_annotated
    base::return(epiRomics_enhanceosome)
}
