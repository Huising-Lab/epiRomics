#' Identifies putative enhancer regions utilizing select histone marks
#'
#' @param epiRomics_dB epiRomics class database containing all data initially loaded
#' @param epiRomics_histone_mark_1 name of first histone mark, must match name in epiROmics_dB@meta, default set to h3k4me1
#' @param epiRomics_histone_mark_2 name of second histone mark, must match name in epiROmics_dB@meta default set to h3k27ac
#' @return Variable of class epiRomics further exploring candidate enhancer regions identified after histone integration
#' @export
epiRomics_enhancers <- function(epiRomics_dB, epiRomics_histone_mark_1 = "h3k4me1", epiRomics_histone_mark_2 = "h3k27ac") {
  epiRomics_putative <- epiRomics_dB
  mark1 <- base::paste0(epiRomics_dB@genome, "_custom_", epiRomics_histone_mark_1)
  mark2 <- base::paste0(epiRomics_dB@genome, "_custom_", epiRomics_histone_mark_2)

  epiRomics_putative@annotations <- GenomicRanges::intersect(GenomicRanges::reduce(epiRomics_dB@annotations[epiRomics_dB@annotations$type ==
    mark1, ]), GenomicRanges::reduce(epiRomics_dB@annotations[epiRomics_dB@annotations$type == mark2, ]))
  base::return(epiRomics_putative)
}
