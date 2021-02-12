#' Filters putative enhancers called by epiRomics_enhancers by crossing against curated FANTOM data
#'
#' @param epiRomics_putative_enhancers epiRomics class database containing putative enhancer calls
#' @param epiRomics_dB epiRomics class database containing all data initially loaded
#' @param epiRomics_type epiRomics reference containing database to validate putative enhancers against
#' @return Variable of class epiRomics with filtered candidate enhancer regions
#' @export


epiRomics_enhancers_filter <- function(epiRomics_putative_enhancers, epiRomics_dB, epiRomics_type="mm10_custom_fantom") {
    epiRomics_enhanceosome_raw_calls <- epiRomics_putative_enhancers@annotations
    epiRomics_functional <- epiRomics_dB@annotations
    # mm10_enhancers_fantom
    epiRomics_putative_enhancers_filtered <- epiRomics_putative_enhancers
    epiRomics_putative_enhancers_filtered@annotations <- IRanges::subsetByOverlaps(epiRomics_enhanceosome_raw_calls,
        epiRomics_functional[epiRomics_functional$type == epiRomics_type,
            ])
    base::return(epiRomics_putative_enhancers_filtered)
}
