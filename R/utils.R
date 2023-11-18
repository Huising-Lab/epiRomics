#' Utils
#'
#' Utilities to normalize bw for visualization
#'
#' @param bw The bigwig file to examine
#'
#' @param gr The genomic ranges to evaluate for bw
#'
#'  @param bws List of bigwig files to evaluate
#'
#' @noRd
maxCovBw <- function(bw, gr) {
  ovlp <- IRanges::subsetByOverlaps(rtracklayer::import(bw, format="BigWig", selection=BigWigSelection(gr)), gr)
  if (length(ovlp) > 0) {
    print("not empty")
    max_cov <- max(ovlp$score)
  } else {
    print("WARNING: The selected genomic region has no coverage value in the BigWig")
    print("WARNING: Coverage value is arbitrary set to Zero.")
    max_cov <- 0
  }
  print(max_cov)
  return(max_cov)
}

maxCovFiles <- function(bws, gr) {
  # bws <- lapply(bws, rtracklayer:::import)
  max_cov <- c()
  for (i in 1:length(gr)) {
    my_feat <- gr[i, ]
    max_cov[i] <- round(
      max(
        sapply(bws, maxCovBw, gr = my_feat)
      ),
      2
    )
  }
  GenomicRanges::values(gr) <- max_cov
  return(gr)
}
