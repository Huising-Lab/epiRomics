#' Analyze statistical significance of TF co-binding at enhanceosome regions
#'
#' Performs pairwise statistical testing of transcription factor co-occurrence
#' at enhanceosome regions using Fisher's exact test or permutation testing,
#' with odds ratios, Pointwise Mutual Information (PMI), and hierarchical
#' clustering.
#'
#' This replaces the previous decision-tree approach (\code{epiRomics_predictors})
#' with statistically rigorous co-binding analysis. For each pair of TFs, a 2x2
#' contingency table is constructed from the enhanceosome presence matrix.
#' P-values are corrected using Benjamini-Hochberg FDR.
#'
#' @section Statistical methods:
#' \itemize{
#'   \item \strong{Fisher's exact test} (method = "fisher"): Tests whether two
#'     TFs co-occur at enhanceosome regions more (or less) often than expected
#'     by chance. Assumes independence between regions. This is the default and
#'     is appropriate when regions are largely non-overlapping.
#'     Reference: Fisher, R.A. (1922) J Royal Stat Soc.
#'   \item \strong{Permutation test} (method = "permutation"): Shuffles TF_B
#'     binding labels across regions to generate a null distribution, accounting
#'     for spatial autocorrelation between nearby genomic regions. More
#'     conservative but robust to violations of independence.
#'     Reference: Gel et al. (2016) Bioinformatics 32(2):289-291. "regioneR:
#'     an R/Bioconductor package for the association analysis of genomic
#'     regions."
#'   \item \strong{Odds ratio}: Measures strength of association. OR > 1
#'     indicates co-occurrence; OR < 1 indicates mutual exclusion.
#'   \item \strong{PMI}: Pointwise Mutual Information quantifies the degree
#'     of association between two TFs: \code{PMI(A,B) = log2(P(A,B) / (P(A)*P(B)))}.
#'     PMI > 0 indicates co-occurrence; PMI < 0 indicates avoidance.
#'     Reference: Church & Hanks (1990) Computational Linguistics.
#'   \item \strong{BH-FDR}: Benjamini-Hochberg correction controls the false
#'     discovery rate across all pairwise tests.
#'     Reference: Benjamini & Hochberg (1995) J Royal Stat Soc B.
#' }
#'
#' @section Note on spatial autocorrelation:
#' Fisher's exact test assumes independence between observations (regions).
#' Nearby genomic regions may be spatially correlated (e.g., broad TF binding
#' domains), which can inflate significance. If your enhanceosome regions
#' contain many closely spaced or overlapping intervals, consider using
#' \code{method = "permutation"} for more conservative p-values.
#'
#' @param epiRomics_enhanceosome epiRomics class database containing enhanceosome calls
#' @param epiRomics_dB epiRomics class database containing all data initially loaded
#' @param fdr_threshold numeric, FDR threshold for significance (default: 0.05)
#' @param min_regions integer, minimum number of co-bound regions to report
#'   a pair (default: 5)
#' @param method character, statistical method: "fisher" (default) for Fisher's
#'   exact test or "permutation" for permutation-based testing that accounts
#'   for spatial autocorrelation.
#' @param n_permutations integer, number of permutations when
#'   \code{method = "permutation"} (default: 1000). Ignored for Fisher's test.
#' @return list with components:
#'   \describe{
#'     \item{pairwise}{data.frame with columns: tf1, tf2, n_both, n_tf1_only,
#'       n_tf2_only, n_neither, odds_ratio, pvalue, fdr, pmi, significant}
#'     \item{presence_matrix}{logical matrix (regions x TFs) of binding presence}
#'     \item{clustering}{hclust object from hierarchical clustering of TF
#'       co-occurrence (Jaccard distance, Ward.D2 linkage)}
#'     \item{tf_names}{character vector of TF names analyzed}
#'     \item{n_regions}{integer, total number of enhanceosome regions}
#'     \item{method}{character, statistical method used}
#'   }
#' @export
#' @seealso \code{\link{epiRomics_tf_overlap}} for overlap fractions without
#'   significance testing
#' @examples
#' db <- epiRomicsS4(annotations = GenomicRanges::GRanges(),
#'   meta = data.frame(name = character(), type = character(),
#'     file = character(), stringsAsFactors = FALSE),
#'   genome = "hg38")
#' tryCatch(epiRomics_tf_cobinding(db, db), error = function(e) message(e$message))
#' \donttest{
#' cobinding <- epiRomics_tf_cobinding(enhanceosome, epiRomics_dB)
#' # Significant co-binding pairs
#' sig <- cobinding$pairwise[cobinding$pairwise$significant, ]
#' sig[order(sig$odds_ratio, decreasing = TRUE), ]
#'
#' # Dendrogram of TF co-occurrence
#' plot(cobinding$clustering)
#'
#' # Permutation test (accounts for spatial autocorrelation)
#' cobinding_perm <- epiRomics_tf_cobinding(enhanceosome, epiRomics_dB,
#'   method = "permutation", n_permutations = 1000)
#' }
epiRomics_tf_cobinding <- function(epiRomics_enhanceosome, epiRomics_dB,
                                    fdr_threshold = 0.05, min_regions = 5L,
                                    method = c("fisher", "permutation"),
                                    n_permutations = 1000L) {
  validate_epiRomics_dB(epiRomics_enhanceosome, "epiRomics_tf_cobinding")
  validate_epiRomics_dB(epiRomics_dB, "epiRomics_tf_cobinding")
  method <- base::match.arg(method)

  if (!base::is.numeric(fdr_threshold) || base::length(fdr_threshold) != 1 ||
      fdr_threshold <= 0 || fdr_threshold > 1) {
    base::stop("fdr_threshold must be a number between 0 and 1")
  }
  if (!base::is.numeric(min_regions) || base::length(min_regions) != 1 || min_regions < 0) {
    base::stop("min_regions must be a non-negative integer")
  }
  if (method == "permutation" && (!base::is.numeric(n_permutations) ||
      base::length(n_permutations) != 1 || n_permutations < 100)) {
    base::stop("n_permutations must be >= 100 for permutation testing")
  }

  genome <- epiRomics_dB@genome
  tf_names <- epiRomics_dB@meta[epiRomics_dB@meta$type == "chip", "name"]
  if (base::length(tf_names) < 2) {
    base::stop("At least 2 ChIP (TF) datasets required for co-binding analysis.")
  }
  n_tf <- base::length(tf_names)
  tf_db_names <- base::paste0(genome, "_custom_", tf_names)

  # Pre-split annotations by type for O(1) lookup
  annot_by_type <- base::split(
    epiRomics_dB@annotations,
    epiRomics_dB@annotations$type
  )

  regions <- epiRomics_enhanceosome@annotations
  n_regions <- base::length(regions)
  if (n_regions == 0) {
    base::stop("No enhanceosome regions found.")
  }

  # Build presence matrix: regions x TFs (parallel when available)
  subject_list <- stats::setNames(
    base::lapply(tf_db_names, function(key) annot_by_type[[key]]),
    tf_names
  )
  presence <- .parallel_count_overlaps(regions, subject_list)

  # Pairwise testing — pre-compute co-occurrence via crossprod
  n_pairs <- base::as.integer(n_tf * (n_tf - 1) / 2)
  co_count <- base::crossprod(presence * 1L)  # all-pairs n_both matrix
  col_sums <- base::colSums(presence)

  results <- base::data.frame(
    tf1 = base::character(n_pairs),
    tf2 = base::character(n_pairs),
    n_both = base::integer(n_pairs),
    n_tf1_only = base::integer(n_pairs),
    n_tf2_only = base::integer(n_pairs),
    n_neither = base::integer(n_pairs),
    odds_ratio = base::numeric(n_pairs),
    pvalue = base::numeric(n_pairs),
    pmi = base::numeric(n_pairs),
    stringsAsFactors = FALSE
  )

  idx <- 1L
  for (i in base::seq_len(n_tf - 1)) {
    for (j in (i + 1):n_tf) {
      n_both <- co_count[i, j]
      n_a_only <- col_sums[i] - n_both
      n_b_only <- col_sums[j] - n_both
      n_neither <- n_regions - col_sums[i] - col_sums[j] + n_both

      results$tf1[idx] <- tf_names[i]
      results$tf2[idx] <- tf_names[j]
      results$n_both[idx] <- n_both
      results$n_tf1_only[idx] <- n_a_only
      results$n_tf2_only[idx] <- n_b_only
      results$n_neither[idx] <- n_neither

      # Odds ratio (always computed)
      contingency <- base::matrix(
        base::c(n_both, n_a_only, n_b_only, n_neither), nrow = 2)
      fisher_result <- stats::fisher.test(contingency)
      results$odds_ratio[idx] <- base::round(
        base::as.numeric(fisher_result$estimate), 4)

      if (method == "fisher") {
        results$pvalue[idx] <- fisher_result$p.value
      } else {
        # Permutation test: shuffle TF_B labels, count co-occurrences
        a <- presence[, i]
        b <- presence[, j]
        observed <- n_both
        perm_counts <- base::replicate(n_permutations, {
          base::sum(a & base::sample(b))
        })
        n_ge <- base::sum(perm_counts >= observed)
        results$pvalue[idx] <- (n_ge + 1) / (n_permutations + 1)
      }

      # Pointwise Mutual Information: log2(P(A,B) / (P(A) * P(B)))
      p_ab <- n_both / n_regions
      p_a <- col_sums[i] / n_regions
      p_b <- col_sums[j] / n_regions
      if (p_ab > 0 && p_a > 0 && p_b > 0) {
        results$pmi[idx] <- base::round(base::log2(p_ab / (p_a * p_b)), 4)
      } else {
        results$pmi[idx] <- NA_real_
      }

      idx <- idx + 1L
    }
  }

  # Benjamini-Hochberg FDR correction
  results$fdr <- stats::p.adjust(results$pvalue, method = "BH")
  results$pvalue <- round(results$pvalue, 6)
  results$fdr <- round(results$fdr, 6)

  # Significance flag
  results$significant <- results$fdr < fdr_threshold & results$n_both >= min_regions

  # Sort by odds ratio (descending) within significant pairs first
  results <- results[order(!results$significant, -results$odds_ratio), ]
  rownames(results) <- NULL

  # Hierarchical clustering on Jaccard distance (vectorized via crossprod)
  presence_num <- presence * 1L
  intersection_mat <- base::crossprod(presence_num)
  col_sums <- base::colSums(presence_num)
  union_mat <- base::outer(col_sums, col_sums, "+") - intersection_mat
  jaccard_dist <- base::ifelse(union_mat > 0, 1 - intersection_mat / union_mat, 1)
  base::rownames(jaccard_dist) <- tf_names
  base::colnames(jaccard_dist) <- tf_names

  clustering <- stats::hclust(stats::as.dist(jaccard_dist), method = "ward.D2")

  return(list(
    pairwise = results,
    presence_matrix = presence,
    clustering = clustering,
    tf_names = tf_names,
    n_regions = n_regions,
    method = method
  ))
}


#' Predict TF behavior at enhanceosome regions (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' Use \code{\link{epiRomics_tf_cobinding}} instead, which provides Fisher's
#' exact test, odds ratios, PMI, and hierarchical clustering for TF co-binding
#' analysis.
#'
#' @param epiRomics_putative_enhanceosome epiRomics class database containing
#'   putative enhanceosome calls
#' @return list (from epiRomics_tf_cobinding). Previously returned a ctree object.
#' @export
#' @examples
#' db <- epiRomicsS4(annotations = GenomicRanges::GRanges(),
#'   meta = data.frame(name = character(), type = character(),
#'     file = character(), stringsAsFactors = FALSE),
#'   genome = "hg38")
#' tryCatch(epiRomics_predictors(db), error = function(e) message(e$message))
#' \donttest{
#' # Deprecated: use epiRomics_tf_cobinding() instead
#' }
epiRomics_predictors <- function(epiRomics_putative_enhanceosome) {
  .Deprecated("epiRomics_tf_cobinding",
    msg = paste0(
      "epiRomics_predictors() is deprecated. Use epiRomics_tf_cobinding() instead.\n",
      "The decision tree approach has been replaced with statistically rigorous\n",
      "co-binding analysis (Fisher's exact test, odds ratios, PMI)."
    )
  )
  base::warning(
    "epiRomics_predictors() requires epiRomics_dB as second argument in the new API. ",
    "Returning NULL. Please call epiRomics_tf_cobinding(enhanceosome, dB) instead."
  )
  return(NULL)
}
