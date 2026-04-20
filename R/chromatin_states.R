#' Classify genomic regions by histone chromatin state with genomic context
#'
#' Given histone mark combinations in the epiRomics database, classifies
#' regions based on curated chromatin state definitions (ChromHMM/Roadmap
#' Epigenomics conventions). States are refined by TSS proximity so that
#' "promoter" labels are assigned only to regions near transcription start
#' sites (within \code{tss_window} bp). Regions with promoter-associated
#' marks (H3K4me3) that fall outside TSS windows are reclassified as
#' enhancers (e.g., "active_enhancer" instead of "active_promoter").
#'
#' Chromatin state definitions (6 simplified labels, priority order):
#' \itemize{
#'   \item \strong{repressed}: H3K27me3 + H3K9me3, or H3K9me3 alone, or
#'     H3K27me3 alone (Polycomb/heterochromatin)
#'   \item \strong{bivalent}: H3K4me3 + H3K27me3 (poised for activation)
#'   \item \strong{active}: H3K4me3 + H3K27ac; H3K4me1 + H3K27ac;
#'     H3K4me1 + H3K27ac + H3K36me3; H3K36me3 alone; H2A.Z + H3K27ac
#'   \item \strong{poised}: H3K4me1 + H3K27me3; H2A.Z + H3K27me3
#'   \item \strong{primed}: H3K4me1 only
#'   \item \strong{unmarked}: no marks, H2A.Z alone, or unclassifiable
#' }
#' Genomic context (promoter/gene_body/intergenic) is reported separately
#' so users can infer regulatory identity from position.
#'
#' @section TSS Refinement:
#' H3K4me3 is a promoter-specific mark that peaks at TSS regions
#' (Santos-Rosa et al. 2002, Bernstein et al. 2005). When H3K4me3-containing
#' states are observed outside TSS windows, they likely represent either:
#' (a) strong/broad enhancers that recruit H3K4me3 (Pekowska et al. 2011),
#' or (b) unannotated alternative promoters. By default, these are
#' reclassified as enhancer states. Set \code{refine_by_tss = FALSE} to
#' disable this behavior.
#'
#' @section Genomic Context:
#' A \code{genomic_context} column is added to the output indicating whether
#' each region is at a "promoter" (within \code{tss_window} of a TSS),
#' "gene_body" (overlapping a gene but not near TSS), or "intergenic"
#' (not overlapping any annotated gene).
#'
#' @section References:
#' \itemize{
#'   \item Ernst J, Kellis M (2012) Nature Methods 9(3):215-216. "ChromHMM:
#'     automating chromatin-state discovery and characterization."
#'   \item Kundaje A et al. (2015) Nature 518(7539):317-330. "Integrative
#'     analysis of 111 reference human epigenomes."
#'   \item Creyghton MP et al. (2010) PNAS 107(50):21931-21936. "Histone
#'     H3K27ac separates active from poised enhancers."
#'   \item Rada-Iglesias A et al. (2011) Nature 470(7333):279-283. "A unique
#'     chromatin signature uncovers early developmental enhancers in humans."
#'   \item Santos-Rosa H et al. (2002) Nature 419(6905):407-411. "Active genes
#'     are tri-methylated at K4 of histone H3." H3K4me3 TSS specificity.
#'   \item Pekowska A et al. (2011) EMBO J 30(20):4198-4210. H3K4me3 at
#'     strong enhancers.
#' }
#'
#' @param database epiRomics class database containing
#'   all data initially loaded
#' @param histone_marks character vector of histone mark
#'   names to use for classification.
#'   Must match names in \code{\link{meta}(database)}. If NULL,
#'   auto-detects from meta.
#' @param regions GRanges object of regions to classify. If NULL, uses all
#'   annotations in the database.
#' @param refine_by_tss logical. If TRUE (default), promoter states are assigned
#'   only to regions within \code{tss_window} of an annotated TSS. Regions with
#'   promoter marks (H3K4me3) outside TSS windows are reclassified as enhancers.
#' @param tss_window integer. Distance in bp around each TSS to define the
#'   promoter zone (default: 2000L). Regions within +/- \code{tss_window} of any
#'   annotated TSS are considered "promoter" context.
#' @return data.frame with columns: seqnames, start, end, chromatin_state,
#'   genomic_context ("promoter"/"gene_body"/"intergenic"),
#'   marks_present (comma-separated), n_marks, is_hotspot
#' @export
#' @examples
#' db <- make_example_database()
#' states <- classify_chromatin_states(db)
#' table(states$chromatin_state)
classify_chromatin_states <- function(
    database, histone_marks = NULL,
    regions = NULL, refine_by_tss = TRUE,
    tss_window = 2000L) {
  setup <- .chromatin_state_setup(database, histone_marks, regions,
                                   refine_by_tss, tss_window,
                                   caller = "classify_chromatin_states")
  if (base::is.data.frame(setup)) return(setup)  # empty regions case

  n_regions <- setup$n_regions

  # Vectorized chromatin state classification using data.table::fcase()
  # Priority order follows the Roadmap Epigenomics 18-state ChromHMM model
  # (Kundaje et al. 2015) extended with H3K36me3 gene body states and
  # H2A.Z histone rules (Giaimo et al. 2019; Lai & Pugh 2017).
  #
  # Output uses 6 simplified labels: active, poised, repressed, bivalent,
  # primed, unmarked. Genomic context (promoter/gene_body/intergenic) is
  # provided separately so users can infer regulatory identity from context.
  states <- data.table::fcase(
    # --- Standard histone rules (highest priority) ---
    setup$has_h3k27me3 &
      setup$has_h3k9me3, "repressed",
    setup$has_h3k4me3 &
      setup$has_h3k27me3, "bivalent",
    setup$has_h3k4me3 &
      setup$has_h3k27ac &
      !setup$has_h3k27me3, "active",
    setup$has_h3k4me1 &
      setup$has_h3k27me3, "poised",
    setup$has_h3k4me1 &
      setup$has_h3k27ac &
      setup$has_h3k36me3 &
      !setup$has_h3k27me3, "active",
    setup$has_h3k4me1 &
      setup$has_h3k27ac &
      !setup$has_h3k27me3, "active",
    # H3K36me3 + H3K27ac (no K4me1/K4me3)
    setup$has_h3k36me3 &
      setup$has_h3k27ac &
      !setup$has_h3k4me1 &
      !setup$has_h3k4me3 &
      !setup$has_h3k27me3, "active",
    # H3K36me3 alone: transcribed gene body
    setup$has_h3k36me3 &
      !setup$has_h3k27ac &
      !setup$has_h3k4me1 &
      !setup$has_h3k4me3 &
      !setup$has_h3k9me3 &
      !setup$has_h3k27me3, "unmarked",
    setup$has_h3k4me1 &
      !setup$has_h3k27ac &
      !setup$has_h3k27me3, "primed",
    # H3K27ac alone: active (Creyghton 2010)
    setup$has_h3k27ac &
      !setup$has_h3k4me1 &
      !setup$has_h3k4me3 &
      !setup$has_h3k36me3 &
      !setup$has_h3k27me3, "active",
    # --- H2A.Z rules ---
    setup$has_h2az &
      setup$has_h3k27ac &
      !setup$has_h3k4me3 &
      !setup$has_h3k27me3, "active",
    setup$has_h2az &
      setup$has_h3k4me3 &
      setup$has_h3k27ac, "active",
    setup$has_h2az &
      setup$has_h3k27me3 &
      !setup$has_h3k4me1, "poised",
    setup$has_h2az &
      !setup$any_standard, "unmarked",
    # --- Remaining standard rules ---
    setup$has_h3k9me3 &
      !setup$has_h3k27me3, "repressed",
    setup$has_h3k27me3 &
      !setup$has_h3k4me1 &
      !setup$has_h3k4me3, "repressed",
    # H3K4me3 alone: poised promoter-like
    setup$has_h3k4me3 &
      !setup$has_h3k27ac &
      !setup$has_h3k27me3, "primed",
    !setup$any_mark, "unmarked",
    base::rep(TRUE, n_regions), "unmarked"
  )

  .chromatin_state_result(setup, states)
}


#' Return chromatin state category definitions
#'
#' Variant of \code{\link{classify_chromatin_states}} that outputs
#' simplified single-state categories. Resolves dual-mark states
#' (bivalent = H3K4me3 + H3K27me3, poised = H3K4me1 + H3K27me3)
#' to the single most likely state, useful for applications that
#' require one-state-per-region classification.
#'
#' @details
#' State resolution rules:
#' \itemize{
#'   \item \strong{active_promoter}: H3K4me3 + H3K27ac (replaces "active" when
#'     K4me3 present; includes former bivalent with K27ac)
#'   \item \strong{active_enhancer}: H3K4me1 + H3K27ac (distal active elements)
#'   \item \strong{active}: H3K27ac alone (sufficient per Creyghton et al. 2010)
#'   \item \strong{primed_enhancer}: H3K4me1 alone (no K27ac, no K27me3)
#'   \item \strong{transcribed}: H3K36me3 alone (SETD2/Pol II gene body mark;
#'     standard epiRomics marks this as "unmarked")
#'   \item \strong{polycomb_repressed}: H3K27me3 (Polycomb Repressive Complex;
#'     absorbs former poised and bivalent-without-K27ac)
#'   \item \strong{heterochromatin}: H3K9me3 (constitutive heterochromatin)
#'   \item \strong{quiescent}: No histone marks detected
#'   \item \strong{low_signal}: Borderline marks or H2A.Z only
#' }
#'
#' Bivalent resolution: K4me3 + K27me3 + K27ac -> active_promoter
#' (K27ac sufficient for active). K4me3 + K27me3 without K27ac ->
#' polycomb_repressed (K27me3 dominates without activating mark).
#'
#' Poised resolution: K4me1 + K27me3 -> polycomb_repressed (K27me3
#' dominates in single-state HMM framework).
#'
#' @section References:
#' Ernst J, Kellis M (2012) Nature Methods 9(3):215-216. "ChromHMM:
#' automating chromatin-state discovery and characterization."
#'
#' @inheritParams classify_chromatin_states
#' @return data.frame with same structure as \code{classify_chromatin_states}
#'   but chromatin_state uses ChromHMM-compatible labels:
#'   active_promoter, active_enhancer, active, primed_enhancer,
#'   transcribed, polycomb_repressed, heterochromatin, quiescent, low_signal
#' @export
#' @seealso \code{\link{classify_chromatin_states}} for the standard
#'   epiRomics classification with bivalent/poised states.
#' @examples
#' db <- make_example_database()
#' cats <- chromatin_state_categories(db)
#' table(cats$chromatin_state)
chromatin_state_categories <- function(
    database, histone_marks = NULL, regions = NULL,
    refine_by_tss = TRUE, tss_window = 2000L) {

  setup <- .chromatin_state_setup(
    database, histone_marks, regions,
    refine_by_tss, tss_window,
    caller = "chromatin_state_categories"
  )
  if (base::is.data.frame(setup)) return(setup)  # empty regions case

  n_regions <- setup$n_regions

  # ChromHMM-compatible states: resolves bivalent/poised to single state
  states <- data.table::fcase(
    # --- Constitutive heterochromatin ---
    setup$has_h3k27me3 &
      setup$has_h3k9me3, "heterochromatin",

    # --- Former bivalent: resolved ---
    setup$has_h3k4me3 &
      setup$has_h3k27me3 &
      setup$has_h3k27ac, "active_promoter",
    setup$has_h3k4me3 &
      setup$has_h3k27me3, "polycomb_repressed",

    # --- Active promoter (K4me3+K27ac) ---
    setup$has_h3k4me3 &
      setup$has_h3k27ac &
      !setup$has_h3k27me3, "active_promoter",

    # --- Former poised: K27me3 dominates ---
    setup$has_h3k4me1 &
      setup$has_h3k27me3, "polycomb_repressed",

    # --- Active enhancer (K4me1+K27ac) ---
    setup$has_h3k4me1 &
      setup$has_h3k27ac &
      setup$has_h3k36me3 &
      !setup$has_h3k27me3, "active_enhancer",
    setup$has_h3k4me1 &
      setup$has_h3k27ac &
      !setup$has_h3k27me3, "active_enhancer",

    # --- Active (K36me3+K27ac, no K4) ---
    setup$has_h3k36me3 &
      setup$has_h3k27ac &
      !setup$has_h3k4me1 &
      !setup$has_h3k4me3 &
      !setup$has_h3k27me3, "active",

    # --- Transcribed gene body ---
    setup$has_h3k36me3 &
      !setup$has_h3k27ac &
      !setup$has_h3k4me1 &
      !setup$has_h3k4me3 &
      !setup$has_h3k9me3 &
      !setup$has_h3k27me3, "transcribed",

    # --- Primed enhancer (K4me1 alone) ---
    setup$has_h3k4me1 &
      !setup$has_h3k27ac &
      !setup$has_h3k27me3, "primed_enhancer",

    # --- Active (K27ac alone) ---
    setup$has_h3k27ac &
      !setup$has_h3k4me1 &
      !setup$has_h3k4me3 &
      !setup$has_h3k36me3 &
      !setup$has_h3k27me3, "active",

    # --- H2A.Z rules ---
    setup$has_h2az &
      setup$has_h3k27ac &
      !setup$has_h3k4me3 &
      !setup$has_h3k27me3, "active",
    setup$has_h2az &
      setup$has_h3k4me3 &
      setup$has_h3k27ac, "active_promoter",
    setup$has_h2az &
      setup$has_h3k27me3 &
      !setup$has_h3k4me1, "polycomb_repressed",
    setup$has_h2az &
      !setup$any_standard, "low_signal",

    # --- Remaining repressive states ---
    setup$has_h3k9me3 &
      !setup$has_h3k27me3, "heterochromatin",
    setup$has_h3k27me3 &
      !setup$has_h3k4me1 &
      !setup$has_h3k4me3, "polycomb_repressed",

    # --- K4me3 alone ---
    setup$has_h3k4me3 &
      !setup$has_h3k27ac &
      !setup$has_h3k27me3, "low_signal",

    # --- No marks ---
    !setup$any_mark, "quiescent",
    base::rep(TRUE, n_regions), "quiescent"
  )

  .chromatin_state_result(setup, states)
}


# ---- Internal helpers (not exported) ----

#' Auto-detect and validate histone marks from epiRomics database meta
#'
#' If histone_marks is NULL, extracts mark names from the database meta rows
#' with type == "histone". Validates that the result is a character vector.
#'
#' @param database epiRomics S4 database object
#' @param histone_marks character vector or NULL
#' @param caller character, calling function name for error messages
#' @return character vector of validated histone mark names
#' @noRd
.detect_histone_marks <- function(database, histone_marks, caller) {
  if (base::is.null(histone_marks)) {
    histone_marks <- database@meta[
      database@meta$type == "histone", "name"
    ]
    if (base::length(histone_marks) == 0) {
      base::stop(
        "No histone marks found in database ",
        "meta. Provide histone_marks parameter."
      )
    }
  }
  if (!base::is.character(histone_marks)) {
    base::stop("histone_marks must be a character vector")
  }
  base::return(histone_marks)
}


#' Resolve regions for chromatin state classification
#'
#' If regions is NULL, builds the union of all mark annotations. If regions
#' is provided, validates it as GRanges. Returns the resolved GRanges or
#' an empty data.frame if n_regions == 0.
#'
#' @param database epiRomics S4 database object
#' @param annot_by_type named list of GRanges, split by annotation type
#' @param mark_db_names character vector of database-prefixed mark names
#' @param histone_marks character vector of histone mark names
#' @param regions GRanges or NULL
#' @return GRanges of regions, or data.frame if empty (early exit)
#' @noRd
.resolve_classification_regions <- function(database, annot_by_type,
                                             mark_db_names, histone_marks,
                                             regions) {
  if (!base::is.null(regions) && !methods::is(regions, "GRanges")) {
    base::stop("regions must be a GRanges object")
  }

  if (base::is.null(regions)) {
    mark_ranges <- base::list()
    for (i in base::seq_along(mark_db_names)) {
      mark_data <- annot_by_type[[mark_db_names[i]]]
      if (!base::is.null(mark_data) && base::length(mark_data) > 0) {
        mark_ranges[[histone_marks[i]]] <- GenomicRanges::reduce(mark_data)
      }
    }
    if (base::length(mark_ranges) == 0) {
      base::stop("No histone mark annotations found in database.")
    }
    all_regions <- base::do.call(c, base::unname(mark_ranges))
    regions <- GenomicRanges::reduce(all_regions)
  }

  if (base::length(regions) == 0) {
    base::warning("No regions to classify. Returning empty data.frame.")
    return(base::data.frame(
      seqnames = base::character(0), start = base::integer(0),
      end = base::integer(0), width = base::integer(0),
      chromatin_state = base::character(0),
      marks_present = base::character(0), n_marks = base::integer(0),
      is_hotspot = base::logical(0), stringsAsFactors = FALSE
    ))
  }

  base::return(regions)
}


#' Build boolean overlap matrix and per-mark vectors
#'
#' Computes an n_regions x n_marks boolean matrix indicating which regions
#' overlap each histone mark, plus individual has_* boolean vectors and
#' aggregate any_standard / any_mark flags.
#'
#' @param regions GRanges of classification regions
#' @param n_regions integer, number of regions
#' @param histone_marks character vector of histone mark names
#' @param mark_db_names character vector of database-prefixed mark names
#' @param annot_by_type named list of GRanges, split by annotation type
#' @return named list with overlap_matrix, has_h3k4me1, has_h3k27ac,
#'   has_h3k4me3, has_h3k27me3, has_h3k9me3, has_h3k36me3, has_h2az,
#'   any_standard, any_mark
#' @noRd
.build_overlap_matrix <- function(regions, n_regions, histone_marks,
                                   mark_db_names, annot_by_type) {
  overlap_matrix <- base::matrix(FALSE, nrow = n_regions,
                                  ncol = base::length(histone_marks))
  base::colnames(overlap_matrix) <- histone_marks

  for (i in base::seq_along(histone_marks)) {
    mark_data <- annot_by_type[[mark_db_names[i]]]
    if (!base::is.null(mark_data) && base::length(mark_data) > 0) {
      overlaps <- GenomicRanges::countOverlaps(regions, mark_data)
      overlap_matrix[, i] <- overlaps > 0
    }
  }

  marks_lower <- base::tolower(histone_marks)
  has_mark_vec <- function(pattern) {
    idx <- base::grep(pattern, marks_lower)
    if (base::length(idx) == 0) return(base::rep(FALSE, n_regions))
    if (base::length(idx) == 1) return(overlap_matrix[, idx])
    return(base::rowSums(overlap_matrix[, idx, drop = FALSE]) > 0L)
  }

  has_h3k4me1 <- has_mark_vec("h3k4me1")
  has_h3k27ac <- has_mark_vec("h3k27ac")
  has_h3k4me3 <- has_mark_vec("h3k4me3")
  has_h3k27me3 <- has_mark_vec("h3k27me3")
  has_h3k9me3 <- has_mark_vec("h3k9me3")
  has_h3k36me3 <- has_mark_vec("h3k36me3")
  has_h2az <- has_mark_vec("h2a\\.?z")

  any_standard <- has_h3k4me1 | has_h3k27ac | has_h3k4me3 | has_h3k27me3 |
                  has_h3k9me3 | has_h3k36me3
  any_mark <- any_standard | has_h2az

  base::list(
    overlap_matrix = overlap_matrix,
    has_h3k4me1 = has_h3k4me1,
    has_h3k27ac = has_h3k27ac,
    has_h3k4me3 = has_h3k4me3,
    has_h3k27me3 = has_h3k27me3,
    has_h3k9me3 = has_h3k9me3,
    has_h3k36me3 = has_h3k36me3,
    has_h2az = has_h2az,
    any_standard = any_standard,
    any_mark = any_mark
  )
}


#' Compute genomic context via TSS proximity refinement
#'
#' Classifies each region as "promoter" (within tss_window of TSS),
#' "gene_body" (overlapping a gene), or "intergenic" (neither).
#'
#' @param regions GRanges of classification regions
#' @param n_regions integer, number of regions
#' @param database epiRomics S4 database object
#' @param refine_by_tss logical, whether to compute TSS proximity
#' @param tss_window integer, bp around TSS for promoter zone
#' @return character vector of genomic context labels
#' @noRd
.compute_tss_proximity <- function(regions, n_regions, database,
                                    refine_by_tss, tss_window) {
  genomic_context <- base::rep("intergenic", n_regions)

  if (refine_by_tss && base::length(database@txdb) > 0 &&
      base::nchar(database@txdb) > 0) {
    txdb_obj <- base::tryCatch(
      resolve_txdb(database@txdb),
      error = function(e) NULL)

    if (!base::is.null(txdb_obj)) {
      # Suppress informational message from GenomicFeatures::genes() about
      # sequences in the TxDb that have no annotated genes (dropped seqlevels).
      all_genes <- base::tryCatch(
        base::suppressMessages(GenomicFeatures::genes(txdb_obj)),
        error = function(e) NULL)

      if (!base::is.null(all_genes) && base::length(all_genes) > 0) {
        # Suppress GenomicRanges warning about strand-unaware resize on "*"
        # strand regions — expected when TxDb contains unstranded entries.
        tss <- GenomicRanges::trim(base::suppressWarnings(
          GenomicRanges::resize(all_genes, width = 1L, fix = "start")))
        # Suppress same upstream warning for TSS-window resize step.
        tss_windows <- GenomicRanges::trim(base::suppressWarnings(
          GenomicRanges::resize(tss, width = tss_window * 2L, fix = "center")))
        tss_windows <- GenomicRanges::reduce(tss_windows)

        near_tss <- GenomicRanges::countOverlaps(regions, tss_windows) > 0
        in_gene_body <- GenomicRanges::countOverlaps(regions, all_genes) > 0
        genomic_context <- base::ifelse(
          near_tss, "promoter",
          base::ifelse(in_gene_body, "gene_body", "intergenic"))
      }
    }
  }

  base::return(genomic_context)
}


#' Compute marks_present, n_marks, is_hotspot, and h2az_overlap columns
#'
#' Builds the marks_present string (comma-separated), counts n_marks per
#' region, flags hotspots (>= 3 marks), and detects H2A.Z overlap.
#'
#' @param overlap_matrix boolean matrix (n_regions x n_marks)
#' @param histone_marks character vector of histone mark names
#' @param n_regions integer, number of regions
#' @param database epiRomics S4 database object
#' @param genome character, genome assembly name
#' @param annot_by_type named list of GRanges, split by annotation type
#' @param regions GRanges of classification regions
#' @return named list with marks_present, n_marks, is_hotspot, h2az_overlap
#' @noRd
.compute_marks_present <- function(overlap_matrix, histone_marks, n_regions,
                                    database, genome, annot_by_type,
                                    regions) {
  # Vectorized marks_present — column-wise with separator-strip pattern
  marks_present <- base::rep("", n_regions)
  for (j in base::seq_along(histone_marks)) {
    has_this <- overlap_matrix[, j]
    if (base::any(has_this)) {
      marks_present[has_this] <- base::paste0(marks_present[has_this],
                                               ",", histone_marks[j])
    }
  }
  # Strip leading comma from all non-empty entries
  has_marks_str <- base::nchar(marks_present) > 0L
  if (base::any(has_marks_str)) {
    marks_present[has_marks_str] <- base::substr(
      marks_present[has_marks_str], 2L,
      base::nchar(marks_present[has_marks_str]))
  }

  n_marks <- base::rowSums(overlap_matrix)
  is_hotspot <- n_marks >= 3L

  # H2A.Z overlap detection
  h2az_overlap <- base::rep(FALSE, n_regions)
  h2az_names <- database@meta[
    database@meta$type == "histone" &
      base::grepl("h2a\\.?z", base::tolower(database@meta$name)), "name"]
  if (base::length(h2az_names) > 0) {
    for (hz in h2az_names) {
      hz_key <- base::paste0(genome, "_custom_", hz)
      hz_data <- annot_by_type[[hz_key]]
      if (!base::is.null(hz_data) && base::length(hz_data) > 0) {
        h2az_overlap <- h2az_overlap |
          (GenomicRanges::countOverlaps(regions, hz_data) > 0)
      }
    }
  }

  base::list(
    marks_present = marks_present,
    n_marks = n_marks,
    is_hotspot = is_hotspot,
    h2az_overlap = h2az_overlap
  )
}


#' Shared setup for chromatin state classification
#'
#' Validates inputs, computes overlap matrix, resolves mark boolean vectors,
#' determines genomic context via TSS refinement, and computes summary columns
#' (marks_present, n_marks, is_hotspot, h2az_overlap).
#' Both public classification functions call this
#' helper then apply their own fcase() rules.
#'
#' @param database epiRomics S4 database object
#' @param histone_marks character vector of histone mark
#'   names, or NULL for auto-detect
#' @param regions GRanges of regions to classify, or NULL for auto-detect
#' @param refine_by_tss logical, whether to compute genomic context
#' @param tss_window integer, bp around TSS for promoter zone
#' @param caller character, name of calling function for error messages
#' @return named list with all pre-computed vectors,
#'   or data.frame if empty regions
#' @noRd
.chromatin_state_setup <- function(database, histone_marks, regions,
                                    refine_by_tss, tss_window, caller) {
  validate_database(database, caller)
  genome <- database@genome

  histone_marks <- .detect_histone_marks(database, histone_marks, caller)
  mark_db_names <- base::paste0(genome, "_custom_", histone_marks)

  # Pre-split annotations by type for O(1) lookup
  annot_by_type <- base::split(
    database@annotations,
    database@annotations$type
  )

  regions <- .resolve_classification_regions(
    database, annot_by_type, mark_db_names, histone_marks, regions
  )
  if (base::is.data.frame(regions)) return(regions)

  n_regions <- base::length(regions)

  overlap_info <- .build_overlap_matrix(
    regions, n_regions, histone_marks, mark_db_names, annot_by_type
  )

  genomic_context <- .compute_tss_proximity(
    regions, n_regions, database, refine_by_tss, tss_window
  )

  summary_cols <- .compute_marks_present(
    overlap_info$overlap_matrix, histone_marks, n_regions,
    database, genome, annot_by_type, regions
  )

  base::c(
    base::list(
      regions = regions,
      n_regions = n_regions,
      histone_marks = histone_marks,
      genomic_context = genomic_context
    ),
    overlap_info,
    summary_cols
  )
}


#' Build result data.frame from setup and classified states
#' @param setup named list from .chromatin_state_setup()
#' @param states character vector of classified chromatin states
#' @return data.frame with chromatin state results
#' @noRd
.chromatin_state_result <- function(setup, states) {
  # Direct column extraction avoids slow as.data.frame(GRanges) dispatch
  regions <- setup$regions
  base::data.frame(
    seqnames = base::as.character(GenomeInfoDb::seqnames(regions)),
    start = BiocGenerics::start(regions),
    end = BiocGenerics::end(regions),
    width = BiocGenerics::width(regions),
    chromatin_state = states,
    genomic_context = setup$genomic_context,
    marks_present = setup$marks_present,
    n_marks = setup$n_marks,
    is_hotspot = setup$is_hotspot,
    h2az_overlap = setup$h2az_overlap,
    stringsAsFactors = FALSE
  )
}


