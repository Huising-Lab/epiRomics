#' Analyze pairwise and multi-way overlap between transcription factor binding sites
#'
#' For N TFs in the enhanceosome, computes pairwise overlap fractions,
#' unique region counts per TF, shared region counts, and UpSet-style
#' intersection data for all TF combinations. Uses Jaccard index for
#' symmetric overlap quantification and overlap coefficient for asymmetric
#' assessment (Church & Hanks, 1990). The presence/absence matrix pattern
#' follows the approach used by DiffBind (Stark & Brown, 2011).
#'
#' @section References:
#' \itemize{
#'   \item Church KW, Hanks P (1990) Computational Linguistics 16(1):22-29.
#'     "Word Association Norms, Mutual Information, and Lexicography."
#'     Pointwise mutual information framework adapted for co-binding analysis.
#'   \item Stark R, Brown GD (2011) DiffBind, Bioconductor.
#'     Presence/absence matrix approach for binding site overlap analysis.
#' }
#'
#' @param epiRomics_enhanceosome epiRomics class database containing enhanceosome calls,
#'   or an epiRomics database with ChIP data
#' @param epiRomics_dB epiRomics class database containing all data initially loaded
#' @param regions GRanges object defining regions to analyze. If NULL, uses all
#'   enhanceosome annotation regions.
#' @return list with components:
#'   \describe{
#'     \item{overlap_matrix}{matrix of pairwise overlap coefficients (asymmetric: fraction of row TF overlapping col TF)}
#'     \item{overlap_counts}{matrix of pairwise absolute overlap counts}
#'     \item{jaccard_matrix}{matrix of pairwise Jaccard indices (symmetric: intersection/union)}
#'     \item{unique_counts}{named integer vector of regions bound by ONLY that TF}
#'     \item{shared_all}{integer count of regions bound by ALL TFs}
#'     \item{n_tf_counts}{table of how many regions are bound by exactly 1, 2, 3... N TFs}
#'     \item{tf_names}{character vector of TF names analyzed}
#'     \item{summary}{data.frame with per-TF summary statistics}
#'   }
#' @export
#' @examples
#' db <- epiRomicsS4(annotations = GenomicRanges::GRanges(),
#'   meta = data.frame(name = character(), type = character(),
#'     file = character(), stringsAsFactors = FALSE),
#'   genome = "hg38")
#' tryCatch(epiRomics_tf_overlap(db, db), error = function(e) message(e$message))
#' \donttest{
#' overlap <- epiRomics_tf_overlap(enhanceosome, epiRomics_dB)
#' # Pairwise overlap matrix
#' overlap$overlap_matrix
#' # Regions unique to each TF
#' overlap$unique_counts
#' # Summary
#' overlap$summary
#' }
epiRomics_tf_overlap <- function(epiRomics_enhanceosome, epiRomics_dB, regions = NULL) {
  validate_epiRomics_dB(epiRomics_enhanceosome, "epiRomics_tf_overlap")
  validate_epiRomics_dB(epiRomics_dB, "epiRomics_tf_overlap")
  genome <- epiRomics_dB@genome

  # Get TF names from meta
  tf_names <- epiRomics_dB@meta[epiRomics_dB@meta$type == "chip", "name"]
  if (base::length(tf_names) == 0) {
    base::stop("No ChIP (TF) data found in database meta.")
  }
  n_tf <- base::length(tf_names)
  tf_db_names <- base::paste0(genome, "_custom_", tf_names)

  # Pre-split annotations by type for O(1) lookup
  annot_by_type <- base::split(
    epiRomics_dB@annotations,
    epiRomics_dB@annotations$type
  )

  # Define regions to analyze
  if (!base::is.null(regions) && !methods::is(regions, "GRanges")) {
    base::stop("regions must be a GRanges object")
  }
  if (base::is.null(regions)) {
    regions <- epiRomics_enhanceosome@annotations
  }
  n_regions <- base::length(regions)
  if (n_regions == 0) {
    base::warning("No regions to analyze.")
  }

  # Build overlap presence matrix: regions x TFs (parallel when available)
  subject_list <- stats::setNames(
    base::lapply(tf_db_names, function(key) annot_by_type[[key]]),
    tf_names
  )
  presence <- .parallel_count_overlaps(regions, subject_list)

  # Vectorized pairwise overlap via matrix algebra (replaces O(n_tf^2) loop)
  presence_num <- presence * 1L
  overlap_abs <- base::crossprod(presence_num)
  base::storage.mode(overlap_abs) <- "integer"
  base::rownames(overlap_abs) <- tf_names
  base::colnames(overlap_abs) <- tf_names

  # Overlap coefficient: fraction of TF_i overlapping TF_j (asymmetric)
  col_sums <- base::colSums(presence_num)
  overlap_frac <- base::round(
    base::sweep(overlap_abs, 1, base::ifelse(col_sums > 0, col_sums, 1), "/"), 4
  )
  base::rownames(overlap_frac) <- tf_names
  base::colnames(overlap_frac) <- tf_names

  # Jaccard index: intersection / union (symmetric)
  union_mat <- base::outer(col_sums, col_sums, "+") - overlap_abs
  jaccard <- base::ifelse(union_mat > 0, base::round(overlap_abs / union_mat, 4), 0)
  base::rownames(jaccard) <- tf_names
  base::colnames(jaccard) <- tf_names

  # Count TFs per region
  tfs_per_region <- base::rowSums(presence)

  # Unique: bound by ONLY that TF (and no others) -- vectorized
  solo_mask <- tfs_per_region == 1
  unique_counts <- base::colSums(presence[solo_mask, , drop = FALSE])
  base::names(unique_counts) <- tf_names

  # Shared by all TFs
  shared_all <- base::sum(tfs_per_region == n_tf)

  # Distribution of TF binding: how many regions bound by exactly k TFs
  n_tf_counts <- base::table(factor(tfs_per_region, levels = 0:n_tf))
  base::names(n_tf_counts) <- base::paste0(0:n_tf, "_TFs")

  # Per-TF summary
  total_bound <- base::colSums(presence)
  summary_df <- base::data.frame(
    tf = tf_names,
    total_regions = base::as.integer(total_bound),
    unique_regions = unique_counts,
    pct_of_all = if (n_regions > 0) base::round(total_bound / n_regions * 100, 1) else base::rep(0, n_tf),
    pct_unique = base::round(unique_counts / base::pmax(total_bound, 1) * 100, 1),
    stringsAsFactors = FALSE
  )
  base::rownames(summary_df) <- NULL

  return(list(
    overlap_matrix = overlap_frac,
    overlap_counts = overlap_abs,
    jaccard_matrix = jaccard,
    unique_counts = unique_counts,
    shared_all = shared_all,
    n_tf_counts = n_tf_counts,
    tf_names = tf_names,
    summary = summary_df
  ))
}
