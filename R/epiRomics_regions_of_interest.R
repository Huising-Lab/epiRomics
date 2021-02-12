#' Evaluates whether regions of interest derived from external experiments, such as ATAC-Seq, correspond with enhanceosome regions
#'
#' @param epiRomics_putative_enhanceosome epiRomics class database containing putative enhanceosome calls
#' @param epiRomics_test_regions GRanges containing regions of interest
#' @return Variable of class epiRomics with enhanceosome regions overlapping with regions of interest
#' @export


epiRomics_regions_of_interest <- function(epiRomics_putative_enhanceosome, epiRomics_test_regions) {
    epiRomics_identified_regions_of_interest <- epiRomics_putative_enhanceosome
    epiRomics_identified_regions_of_interest@annotations <- IRanges::subsetByOverlaps(epiRomics_putative_enhanceosome@annotations, 
        epiRomics_test_regions)
    
    base::return(epiRomics_identified_regions_of_interest)
}
