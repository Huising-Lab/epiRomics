#' Interrogates various histone marks against a curated database to determine which are most informative
#'
#' @param epiRomics_dB epiRomics class database containing all data initially loaded
#' @param epiRomics_histone name or vector of histone mark(s), must match name in epiROmics_dB@meta, default set to h3k4me1
#' @param epiRomics_curated_database database to test histone marks against, must match name in epiROmics_dB@meta default set to fantom
#' @return Variable of class dataframe further exploring top histone marks that may determine enhancer regions
#' @export
epiRomics_enhancer_predictor_test <- function(epiRomics_dB, epiRomics_histone = "h3k4me1", epiRomics_curated_database = "fantom") {
  epiRomics_putative <- epiRomics_dB
  histone_test <- base::matrix(data = NA, nrow = base::length(epiRomics_histone), ncol = 2)
  colnames(histone_test) <- base::c("Histone_Mark", "Fraction_of_Overlap")
  histone_test[, 1] <- epiRomics_histone
  mark2 <- base::paste0(epiRomics_dB@genome, "_custom_", epiRomics_curated_database)
  mark2_count <- base::length(epiRomics_dB@annotations[epiRomics_dB@annotations$type == mark2, ])
  for (i in 1:base::length(epiRomics_histone)) {
    mark1 <- base::paste0(epiRomics_dB@genome, "_custom_", epiRomics_histone[i])
    epiRomics_putative@annotations <-
      GenomicRanges::intersect(
        GenomicRanges::reduce(epiRomics_dB@annotations[epiRomics_dB@annotations$type ==
          mark1, ]),
        GenomicRanges::reduce(epiRomics_dB@annotations[epiRomics_dB@annotations$type == mark2, ])
      )
    mark1_count <- base::length(epiRomics_putative@annotations)

    histone_test[i, 2] <- mark1_count / mark2_count
  }
  histone_test <- BiocGenerics::as.data.frame(histone_test)
  histone_test <- histone_test[base::order(histone_test[, 2], decreasing = TRUE), ]
  base::return(histone_test)
}
