#' Parallel countOverlaps for multiple subject GRanges
#'
#' Computes boolean overlap columns in parallel when \code{parallel} package
#' is available. Falls back to sequential \code{for} loop otherwise.
#'
#' @param query GRanges of query regions
#' @param subject_list named list of GRanges objects to test against
#' @return logical matrix with nrow = length(query), ncol = length(subject_list)
#' @noRd
.parallel_count_overlaps <- function(query, subject_list) {
  n_q <- base::length(query)
  n_s <- base::length(subject_list)
  s_names <- base::names(subject_list)

  if (n_s == 0L) {
    m <- base::matrix(FALSE, nrow = n_q, ncol = 0L)
    return(m)
  }

  # Use parallel when available and >2 subjects (overhead not worth it for 1-2)
  use_par <- n_s > 2L &&
    base::requireNamespace("parallel", quietly = TRUE) &&
    .Platform$OS.type == "unix"

  if (use_par) {
    n_cores <- base::min(n_s, base::max(1L, parallel::detectCores() - 1L))
    cols <- parallel::mclapply(subject_list, function(s_gr) {
      if (base::is.null(s_gr) || base::length(s_gr) == 0L) {
        return(base::rep(FALSE, n_q))
      }
      GenomicRanges::countOverlaps(query, s_gr) > 0L
    }, mc.cores = n_cores)
  } else {
    cols <- base::lapply(subject_list, function(s_gr) {
      if (base::is.null(s_gr) || base::length(s_gr) == 0L) {
        return(base::rep(FALSE, n_q))
      }
      GenomicRanges::countOverlaps(query, s_gr) > 0L
    })
  }

  m <- base::do.call(base::cbind, cols)
  base::colnames(m) <- s_names
  m
}

#' Identify putative enhancer regions using rule-based histone logic
#'
#' Automatically scans all histone and histone variant marks in the epiRomics
#' database and applies ChromHMM-based classification rules to identify putative
#' enhancer regions.
#'
#' This function uses a multi-source approach:
#' \describe{
#'   \item{Chromatin states}{Leverages \code{\link{epiRomics_chromatin_states}}
#'     to classify all genomic regions covered by histone marks. Regions
#'     classified as enhancer-related states are included.}
#'   \item{Hi-C contacts}{If provided, Hi-C contact anchors are added as
#'     putative enhancers. Anchors are classified using the available histone
#'     data at each anchor; anchors with no histone coverage are labeled
#'     \code{"Unmarked"}.}
#'   \item{TF binding}{Regions bound by at least one TF (type = "chip") are
#'     included as putative enhancers. TF binding alone yields \code{"Unmarked"}
#'     chromatin state.}
#'   \item{H2A.Z regions}{H2A.Z-positive regions that were classified as
#'     \code{"unmarked"} by chromatin states are
#'     recovered as putative enhancers, since H2A.Z is
#'     enriched at regulatory elements (Giaimo et al. 2019; Lai & Pugh 2017).
#'     H2A.Z alone is insufficient for specific chromatin state assignment.}
#' }
#'
#' Unlike earlier versions that required the user to specify exactly two
#' histone marks, this function automatically uses ALL histone marks in
#' the database and applies the full set of classification rules.
#'
#' @param epiRomics_dB An epiRomics S4 database object.
#' @param chromatin_states data.frame or NULL. Pre-computed output from
#'   \code{\link{epiRomics_chromatin_states}}. If NULL (default), computed
#'   automatically from all available histone marks.
#' @param hic_contacts data.frame or NULL. Hi-C contacts in BEDPE format.
#'   Anchors are added as putative enhancers and classified using available
#'   histone data.
#' @param enhancer_states Character vector. Chromatin states to include as
#'   putative enhancers. Default includes all enhancer-related states.
#' @return A data.frame with columns:
#'   \describe{
#'     \item{putative_id}{Integer. Unique enhancer index (sorted by
#'       TF co-binding then histone marks).}
#'     \item{chr}{Character. Chromosome.}
#'     \item{start}{Integer. Start position.}
#'     \item{end}{Integer. End position.}
#'     \item{width}{Integer. Region width.}
#'     \item{source}{Character. Origin: \code{"histone"}, \code{"hic"},
#'       \code{"tf"}, or comma-separated if multiple sources contribute.}
#'     \item{chromatin_state}{Character. Broad state category:
#'       Active, Poised, Repressed, or Unmarked.}
#'     \item{chromatin_state_detail}{Character. Specific state from
#'       \code{\link{epiRomics_chromatin_states}} (e.g. active_enhancer,
#'       poised_enhancer, primed_enhancer).}
#'     \item{histone_marks}{Character. Comma-separated histone marks
#'       overlapping this region.}
#'     \item{n_histone_marks}{Integer. Number of histone marks.}
#'     \item{h2az}{Logical. Whether H2A.Z overlaps this region.}
#'     \item{tf_names}{Character. Comma-separated TF names with peaks
#'       overlapping this region (H2A.Z excluded from TF count).}
#'     \item{n_tfs}{Integer. Number of TFs with binding peaks.}
#'   }
#' @export
#' @examples
#' db <- epiRomicsS4(annotations = GenomicRanges::GRanges(),
#'   meta = data.frame(name = character(), type = character(),
#'     file = character(), stringsAsFactors = FALSE),
#'   genome = "hg38")
#' tryCatch(epiRomics_putative_enhancers(db),
#'   error = function(e) message(e$message))
#' \donttest{
#' ## Fully automatic -- uses all available histones
#' pe <- epiRomics_putative_enhancers(epiRomics_dB)
#'
#' ## With pre-computed chromatin states
#' cs <- epiRomics_chromatin_states(epiRomics_dB)
#' pe <- epiRomics_putative_enhancers(epiRomics_dB, chromatin_states = cs)
#'
#' ## With Hi-C contacts -- anchors become putative enhancers too
#' contacts <- read.delim("path/to/contacts.bedpe", header = FALSE)
#' pe <- epiRomics_putative_enhancers(epiRomics_dB, hic_contacts = contacts)
#'
#' ## Browse results
#' head(pe[pe$chromatin_state == "Active", ])
#' pe[pe$n_tfs >= 3, ]  # high TF co-binding
#' pe[pe$h2az, ]          # H2A.Z-positive enhancers
#' }
epiRomics_putative_enhancers <- function(epiRomics_dB,
                                          chromatin_states = NULL,
                                          hic_contacts = NULL,
                                          enhancer_states = base::c(
                                            "active",
                                            "poised",
                                            "primed")) {

  validate_epiRomics_dB(epiRomics_dB, "epiRomics_putative_enhancers")

  genome <- epiRomics_dB@genome
  db_prefix <- base::paste0(genome, "_custom_")
  source_gr_list <- base::list()
  source_label_list <- base::list()

  # Pre-split annotations by type for O(1) lookup (avoids repeated == scans)
  annot_by_type <- base::split(
    epiRomics_dB@annotations,
    epiRomics_dB@annotations$type
  )

  # ================================================================
  # SOURCE 1: Rule-based histone classification (automatic)
  # ================================================================
  if (base::is.null(chromatin_states)) {
    base::message("Auto-computing chromatin states from all histone marks...")
    chromatin_states <- epiRomics_chromatin_states(epiRomics_dB)
  }

  if (!base::is.data.frame(chromatin_states) ||
      base::nrow(chromatin_states) == 0) {
    base::stop("chromatin_states is empty or not a data.frame")
  }

  # Select enhancer-related regions
  enh_rows <- chromatin_states$chromatin_state %in% enhancer_states

  # Also include H2A.Z-positive regions NOT already classified as enhancer
  # H2A.Z alone = primed enhancer; H2A.Z can mark regulatory elements
  # even without traditional histone marks
  has_h2az_col <- "h2az_overlap" %in% base::names(chromatin_states)
  if (has_h2az_col) {
    # Include H2A.Z regions that aren't at promoters (use genomic_context)
    non_promoter <- !base::grepl("promoter",
      chromatin_states$genomic_context)
    h2az_extra <- chromatin_states$h2az_overlap & non_promoter & !enh_rows
    enh_rows <- enh_rows | h2az_extra

    n_h2az_extra <- base::sum(h2az_extra, na.rm = TRUE)
    if (n_h2az_extra > 0) {
      base::message(base::sprintf(
        "H2A.Z: recovered %d additional regions as putative enhancers (H2A.Z-positive, not classified as enhancer by standard histone rules)",
        n_h2az_extra))
    }
  }

  enh_data <- chromatin_states[enh_rows, , drop = FALSE]

  if (base::nrow(enh_data) > 0) {
    histone_gr <- GenomicRanges::GRanges(
      seqnames = enh_data$seqnames,
      ranges = IRanges::IRanges(start = enh_data$start, end = enh_data$end))
    histone_gr <- GenomicRanges::reduce(histone_gr)

    source_gr_list[["histone"]] <- histone_gr
    source_label_list[["histone"]] <- base::rep("histone",
      base::length(histone_gr))

    base::message(base::sprintf(
      "Histone-based enhancers: %d regions (from %d classified regions)",
      base::length(histone_gr), base::nrow(enh_data)))
    # State breakdown
    state_tab <- base::table(enh_data$chromatin_state)
    base::message(base::sprintf("  States: %s", base::paste(
      base::paste0(base::names(state_tab), "=", state_tab),
      collapse = ", ")))
  } else {
    base::message("No enhancer-like regions found from histone marks")
  }

  # ================================================================
  # SOURCE 2: Hi-C contact anchors
  # ================================================================
  if (!base::is.null(hic_contacts) && base::is.data.frame(hic_contacts) &&
      base::nrow(hic_contacts) > 0) {

    anchor1 <- GenomicRanges::GRanges(
      seqnames = hic_contacts$chr1,
      ranges = IRanges::IRanges(start = hic_contacts$x1,
        end = hic_contacts$x2))
    anchor2 <- GenomicRanges::GRanges(
      seqnames = hic_contacts$chr2,
      ranges = IRanges::IRanges(start = hic_contacts$y1,
        end = hic_contacts$y2))
    hic_anchors <- base::c(anchor1, anchor2)
    hic_anchors <- GenomicRanges::reduce(hic_anchors)

    source_gr_list[["hic"]] <- hic_anchors
    source_label_list[["hic"]] <- base::rep("hic",
      base::length(hic_anchors))

    base::message(base::sprintf(
      "Hi-C source: %d unique anchor regions from %d contacts",
      base::length(hic_anchors), base::nrow(hic_contacts)))
  }

  # ================================================================
  # SOURCE 3: TF binding regions (ChIP-seq peaks, excluding H2A.Z)
  # ================================================================
  # Any region with at least one TF bound is a putative regulatory element.
  # TF binding alone does not define chromatin state (classified as "Unmarked").
  tf_meta <- epiRomics_dB@meta[epiRomics_dB@meta$type == "chip", , drop = FALSE]
  tf_chip_names <- tf_meta$name[!base::grepl("h2a\\.?z",
    base::tolower(tf_meta$name))]

  if (base::length(tf_chip_names) > 0) {
    tf_gr_list <- base::list()
    for (tn in tf_chip_names) {
      tk <- base::paste0(db_prefix, tn)
      tf_data <- annot_by_type[[tk]]
      if (!base::is.null(tf_data) && base::length(tf_data) > 0) {
        tf_gr_list[[tn]] <- tf_data
      }
    }
    if (base::length(tf_gr_list) > 0) {
      tf_all <- base::do.call(c, base::unname(tf_gr_list))
      tf_reduced <- GenomicRanges::reduce(tf_all)
      source_gr_list[["tf"]] <- tf_reduced
      source_label_list[["tf"]] <- base::rep("tf",
        base::length(tf_reduced))
      base::message(base::sprintf(
        "TF binding source: %d regions from %d TFs (%s)",
        base::length(tf_reduced), base::length(tf_chip_names),
        base::paste(tf_chip_names, collapse = ", ")))
    }
  }

  # ================================================================
  # UNION all regions
  # ================================================================
  if (base::length(source_gr_list) == 0) {
    base::stop("No putative enhancer regions found from any source")
  }

  # Concatenate all GRanges from the list (single dispatch)
  all_gr <- base::do.call(c, base::unname(source_gr_list))
  all_labels <- base::do.call(base::c, base::unname(source_label_list))

  merged <- GenomicRanges::reduce(all_gr)
  n_merged <- base::length(merged)

  base::message(base::sprintf(
    "Union: %d non-overlapping putative enhancer regions", n_merged))

  # ================================================================
  # ANNOTATE: source attribution (vectorized)
  # ================================================================
  source_sets <- base::split(base::seq_along(all_gr), all_labels)
  # Build a boolean matrix: merged[i] overlaps source[j]?
  src_names <- base::names(source_sets)
  src_hits <- base::matrix(FALSE, nrow = n_merged,
    ncol = base::length(src_names))
  base::colnames(src_hits) <- src_names
  for (j in base::seq_along(src_names)) {
    src_gr <- all_gr[source_sets[[src_names[j]]]]
    src_hits[, j] <- GenomicRanges::countOverlaps(merged, src_gr) > 0
  }
  # Column-wise source string building
  sorted_srcs <- base::sort(src_names)
  merged_sources <- base::rep("", n_merged)
  for (sn in sorted_srcs) {
    idx <- base::which(src_hits[, sn])
    if (base::length(idx) > 0) {
      merged_sources[idx] <- base::ifelse(
        merged_sources[idx] == "", sn,
        base::paste0(merged_sources[idx], ",", sn))
    }
  }

  # ================================================================
  # ANNOTATE: chromatin state (vectorized via findOverlaps)
  # ================================================================
  state_to_broad <- base::c(
    active_enhancer = "Active", strong_enhancer = "Active",
    active_promoter = "Active",
    poised_enhancer = "Poised", primed_enhancer = "Poised",
    bivalent_promoter = "Poised",
    repressed = "Repressed", heterochromatin = "Repressed",
    quiescent = "Unmarked", unclassified = "Unmarked",
    active_gene_body = "Unmarked"
  )
  priority_order <- base::c(Active = 4, Poised = 3, Repressed = 2, Unmarked = 1)

  cs_gr <- GenomicRanges::GRanges(
    seqnames = chromatin_states$seqnames,
    ranges = IRanges::IRanges(
      start = chromatin_states$start, end = chromatin_states$end))
  cs_detail <- chromatin_states$chromatin_state
  cs_broad <- state_to_broad[cs_detail]
  cs_broad[base::is.na(cs_broad)] <- "Unmarked"

  broad_col <- base::rep("Unmarked", n_merged)
  detail_col <- base::rep("unmarked", n_merged)

  cs_hits <- GenomicRanges::findOverlaps(merged, cs_gr)
  q_idx <- S4Vectors::queryHits(cs_hits)
  s_idx <- S4Vectors::subjectHits(cs_hits)

  # For each merged region, find best broad + detail state (data.table grouping)
  hit_dt <- data.table::data.table(
    q = q_idx,
    broad = cs_broad[s_idx],
    detail = cs_detail[s_idx],
    broad_score = priority_order[cs_broad[s_idx]])
  best_dt <- hit_dt[, {
    best_b <- broad[base::which.max(broad_score)]
    enh_d <- detail[detail %in% enhancer_states]
    best_d <- if (base::length(enh_d) > 0) enh_d[1L] else detail[1L]
    base::list(best_broad = best_b, best_detail = best_d)
  }, by = q]
  broad_col[best_dt$q] <- best_dt$best_broad
  detail_col[best_dt$q] <- best_dt$best_detail

  # ================================================================
  # ANNOTATE: histone marks + H2A.Z (vectorized)
  # ================================================================
  hist_names <- epiRomics_dB@meta$name[epiRomics_dB@meta$type == "histone"]

  # Build boolean matrix: merged[i] overlaps histone[j]? (parallel when available)
  hist_subjects <- stats::setNames(base::lapply(hist_names, function(hn) {
    annot_by_type[[base::paste0(db_prefix, hn)]]
  }), hist_names)
  hist_matrix <- .parallel_count_overlaps(merged, hist_subjects)

  # Column-wise string building (much faster than row-wise apply on 180k+ rows)
  sorted_hist <- base::sort(hist_names)
  hist_overlap_names <- base::rep("", n_merged)
  for (hn in sorted_hist) {
    idx <- base::which(hist_matrix[, hn])
    if (base::length(idx) > 0) {
      hist_overlap_names[idx] <- base::ifelse(
        hist_overlap_names[idx] == "", hn,
        base::paste0(hist_overlap_names[idx], ",", hn))
    }
  }
  hist_overlap_counts <- base::rowSums(hist_matrix)

  # H2A.Z detection
  h2az_idx <- base::grep("h2a\\.?z", base::tolower(hist_names))
  if (base::length(h2az_idx) > 0) {
    h2az_col <- base::rowSums(hist_matrix[, h2az_idx, drop = FALSE]) > 0
  } else {
    h2az_col <- base::rep(FALSE, n_merged)
  }

  # ================================================================
  # ANNOTATE: TF co-binding (vectorized, H2A.Z excluded)
  # ================================================================
  tf_names <- epiRomics_dB@meta$name[epiRomics_dB@meta$type == "chip"]
  tf_names <- tf_names[!base::grepl("h2a\\.?z", base::tolower(tf_names))]

  # Build boolean matrix: merged[i] overlaps TF[j]? (parallel when available)
  if (base::length(tf_names) > 0) {
    tf_subjects <- stats::setNames(base::lapply(tf_names, function(tn) {
      annot_by_type[[base::paste0(db_prefix, tn)]]
    }), tf_names)
    tf_matrix <- .parallel_count_overlaps(merged, tf_subjects)
  } else {
    tf_matrix <- base::matrix(FALSE, nrow = n_merged, ncol = 0L)
  }

  if (base::length(tf_names) > 0) {
    sorted_tfs <- base::sort(tf_names)
    tf_overlap_names <- base::rep("", n_merged)
    for (tn in sorted_tfs) {
      idx <- base::which(tf_matrix[, tn])
      if (base::length(idx) > 0) {
        tf_overlap_names[idx] <- base::ifelse(
          tf_overlap_names[idx] == "", tn,
          base::paste0(tf_overlap_names[idx], ",", tn))
      }
    }
    tf_overlap_counts <- base::rowSums(tf_matrix)
  } else {
    tf_overlap_names <- base::character(n_merged)
    tf_overlap_counts <- base::integer(n_merged)
  }

  # ================================================================
  # BUILD OUTPUT
  # ================================================================
  result <- base::data.frame(
    putative_id = base::seq_len(n_merged),
    chr = base::as.character(GenomeInfoDb::seqnames(merged)),
    start = BiocGenerics::start(merged),
    end = BiocGenerics::end(merged),
    width = GenomicRanges::width(merged),
    source = merged_sources,
    chromatin_state = broad_col,
    chromatin_state_detail = detail_col,
    histone_marks = hist_overlap_names,
    n_histone_marks = hist_overlap_counts,
    h2az = h2az_col,
    tf_names = tf_overlap_names,
    n_tfs = tf_overlap_counts,
    stringsAsFactors = FALSE
  )

  # Sort: Active > Poised > Unmarked, then by TF co-binding
  state_sort <- base::c(Active = 1, Poised = 2, Repressed = 3, Unmarked = 4)
  result$sort_key <- state_sort[result$chromatin_state]
  result$sort_key[base::is.na(result$sort_key)] <- 5
  result <- result[base::order(result$sort_key, -result$n_tfs,
    -result$n_histone_marks), , drop = FALSE]
  result$sort_key <- NULL
  result$putative_id <- base::seq_len(base::nrow(result))
  base::rownames(result) <- NULL

  base::message(base::sprintf(
    "Putative enhancers: %d total", base::nrow(result)))
  base::message(base::sprintf(
    "  Active: %d | Poised: %d | Unmarked: %d | Repressed: %d",
    base::sum(result$chromatin_state == "Active"),
    base::sum(result$chromatin_state == "Poised"),
    base::sum(result$chromatin_state == "Unmarked"),
    base::sum(result$chromatin_state == "Repressed")))
  base::message(base::sprintf(
    "  H2A.Z+: %d | TF-bound: %d | High co-binding (>=3 TFs): %d",
    base::sum(result$h2az),
    base::sum(result$n_tfs > 0),
    base::sum(result$n_tfs >= 3)))

  # Source breakdown
  src_counts <- base::table(result$source)
  base::message(base::sprintf(
    "  Sources: %s", base::paste(
      base::paste0(base::names(src_counts), "=", src_counts),
      collapse = ", ")))

  result
}


#' Annotate putative enhancers with functional database overlaps
#'
#' Cross-references putative enhancers against all functional databases
#' (FANTOM5, UCNEs, Regulome Active, Regulome Super, etc.) loaded in the
#' epiRomics database. Adds boolean overlap columns for each database, a
#' count of overlapping databases, and a \code{novel} flag for enhancers
#' not found in any known database.
#'
#' @param putative_enhancers data.frame. Output from
#'   \code{\link{epiRomics_putative_enhancers}}.
#' @param epiRomics_dB An epiRomics S4 database object.
#' @return The input data.frame with additional columns:
#'   \describe{
#'     \item{overlap_<name>}{Logical. TRUE if putative enhancer overlaps the
#'       named functional database (one column per database).}
#'     \item{n_databases}{Integer. Count of functional databases overlapping.}
#'     \item{novel}{Logical. TRUE if no functional database overlap (completely
#'       novel enhancer call).}
#'   }
#' @export
#' @examples
#' db <- epiRomicsS4(annotations = GenomicRanges::GRanges(),
#'   meta = data.frame(name = character(), type = character(),
#'     file = character(), stringsAsFactors = FALSE),
#'   genome = "hg38")
#' pe_df <- data.frame(chr = "chr1", start = 1000L, end = 2000L)
#' tryCatch(epiRomics_annotate_putative(pe_df, db),
#'   error = function(e) message(e$message))
#' \donttest{
#' pe <- epiRomics_putative_enhancers(epiRomics_dB)
#' pe_annot <- epiRomics_annotate_putative(pe, epiRomics_dB)
#' table(pe_annot$novel)       # novel vs validated
#' pe_annot[pe_annot$n_databases >= 2, ]  # multi-database validated
#' }
epiRomics_annotate_putative <- function(putative_enhancers, epiRomics_dB) {

  if (!base::is.data.frame(putative_enhancers) ||
      base::nrow(putative_enhancers) == 0) {
    base::stop("putative_enhancers must be a non-empty data.frame")
  }
  validate_epiRomics_dB(epiRomics_dB, "epiRomics_annotate_putative")

  genome <- epiRomics_dB@genome
  db_prefix <- base::paste0(genome, "_custom_")
  meta <- epiRomics_dB@meta

  # Pre-split annotations by type for O(1) lookup
  annot_by_type <- base::split(
    epiRomics_dB@annotations,
    epiRomics_dB@annotations$type
  )

  # Build GRanges from putative enhancers
  pe_gr <- GenomicRanges::GRanges(
    seqnames = putative_enhancers$chr,
    ranges = IRanges::IRanges(
      start = putative_enhancers$start,
      end = putative_enhancers$end))

  # Get all functional database entries
  func_meta <- meta[meta$type == "functional", , drop = FALSE]

  if (base::nrow(func_meta) == 0) {
    base::message("No functional databases found in epiRomics_dB")
    putative_enhancers$n_databases <- 0L
    putative_enhancers$novel <- TRUE
    return(putative_enhancers)
  }

  n_pe <- base::nrow(putative_enhancers)
  overlap_cols <- base::vector("list", base::nrow(func_meta))
  col_names <- base::character(base::nrow(func_meta))

  for (i in base::seq_len(base::nrow(func_meta))) {
    db_name <- func_meta$name[i]
    db_key <- base::paste0(db_prefix, db_name)
    db_data <- annot_by_type[[db_key]]

    col_name <- base::paste0("overlap_", db_name)
    col_names[i] <- col_name
    if (!base::is.null(db_data) && base::length(db_data) > 0) {
      hits <- GenomicRanges::countOverlaps(pe_gr, db_data) > 0
      overlap_cols[[i]] <- hits
      n_hit <- base::sum(hits)
      base::message(base::sprintf(
        "  %s: %d / %d (%s%%) overlap",
        db_name, n_hit, n_pe, base::round(n_hit / n_pe * 100, 1)))
    } else {
      overlap_cols[[i]] <- base::rep(FALSE, n_pe)
      base::message(base::sprintf("  %s: 0 annotations in DB", db_name))
    }
  }
  base::names(overlap_cols) <- col_names

  # Add overlap columns to data.frame (single cbind avoids repeated copy)
  for (cn in col_names) {
    putative_enhancers[[cn]] <- overlap_cols[[cn]]
  }

  # Count databases overlapping and novel flag
  overlap_matrix <- base::do.call(base::cbind, overlap_cols)
  putative_enhancers$n_databases <- base::rowSums(overlap_matrix)
  putative_enhancers$novel <- putative_enhancers$n_databases == 0L

  # Summary
  n_novel <- base::sum(putative_enhancers$novel)
  n_validated <- n_pe - n_novel
  base::message(base::sprintf(
    "Functional annotation: %d validated (%s%%), %d novel (%s%%)",
    n_validated, base::round(n_validated / n_pe * 100, 1),
    n_novel, base::round(n_novel / n_pe * 100, 1)))
  db_dist <- base::table(putative_enhancers$n_databases)
  base::message(base::sprintf(
    "  Database overlap distribution: %s", base::paste(
      base::paste0(base::names(db_dist), " DBs=", db_dist),
      collapse = ", ")))

  putative_enhancers
}


#' Filter putative enhancers by chromatin accessibility evidence
#'
#' Multi-mode accessibility filter for putative enhancers. Supports four
#' complementary evidence types that can be used independently or combined.
#'
#' @section Modes:
#' \describe{
#'   \item{signal}{Import ATAC-seq/DNase-seq BigWig signal over each region.
#'     Regions with mean signal above a z-score threshold are flagged accessible.
#'     Requires \code{epiRomics_track_connection} with ATAC/DNase tracks.}
#'   \item{bed}{Overlap with an external accessibility BED file (e.g., ENCODE
#'     peaks, DHS hotspots). Any region overlapping a BED entry is flagged.
#'     Requires \code{bed_path}.}
#'   \item{genelist}{Retain enhancers linked to expressed genes. Regions whose
#'     \code{SYMBOL} column matches a gene in \code{gene_list} are flagged.
#'     Requires \code{gene_list}.}
#'   \item{combined}{Union of all available evidence. A region is retained if
#'     ANY mode flags it as accessible.}
#' }
#'
#' @section Scope:
#' The \code{scope} parameter controls which regions are evaluated:
#' \describe{
#'   \item{filter_distal}{Only distal (non-promoter) enhancers are filtered;
#'     promoter-proximal regions are always retained. (Default)}
#'   \item{filter_all}{All regions are subject to filtering, including
#'     promoter-proximal ones.}
#' }
#'
#' @param putative_enhancers data.frame. Output from
#'   \code{\link{epiRomics_putative_enhancers}}. Must contain columns
#'   \code{chr}, \code{start}, \code{end}. For \code{genelist} mode, also
#'   requires a \code{SYMBOL} column.
#' @param epiRomics_track_connection data.frame or NULL. BigWig track
#'   connection sheet. Required for \code{mode = "signal"}.
#' @param mode Character. Filtering mode: \code{"signal"} (default),
#'   \code{"bed"}, \code{"genelist"}, or \code{"combined"}.
#' @param scope Character. Filtering scope: \code{"filter_distal"} (default)
#'   or \code{"filter_all"}.
#' @param signal_threshold Numeric. Z-score threshold for signal mode
#'   (default 2).
#' @param bed_path Character or NULL. Path to a BED file for \code{bed} mode.
#' @param gene_list Character vector or NULL. Expressed gene symbols for
#'   \code{genelist} mode.
#' @param promoter_distance Integer. Distance from TSS to classify as
#'   promoter-proximal (default 2000 bp). Only used when
#'   \code{scope = "filter_distal"}.
#' @return The input data.frame with an \code{atac_accessible} logical column.
#'   For \code{signal} mode, also includes per-sample signal and accessibility
#'   columns.
#' @export
#' @examples
#' pe_df <- data.frame(chr = "chr1", start = 1000L, end = 2000L,
#'   SYMBOL = "GENE1", stringsAsFactors = FALSE)
#' tryCatch(epiRomics_filter_accessible(pe_df, mode = "genelist",
#'   gene_list = c("GENE1")),
#'   error = function(e) message(e$message))
#' \donttest{
#' # Mode 1: BigWig signal
#' pe_signal <- epiRomics_filter_accessible(pe, track_connection,
#'   mode = "signal", scope = "filter_distal")
#'
#' # Mode 2: BED file overlap
#' pe_bed <- epiRomics_filter_accessible(pe, mode = "bed",
#'   bed_path = "ENCODE_DHS_peaks.bed", scope = "filter_all")
#'
#' # Mode 3: Gene list (from RNA-seq)
#' pe_genes <- epiRomics_filter_accessible(pe, mode = "genelist",
#'   gene_list = c("INS", "GCG", "SST", "PPY"))
#'
#' # Mode 4: Combined — accessible if ANY evidence
#' pe_combined <- epiRomics_filter_accessible(pe, track_connection,
#'   mode = "combined", bed_path = "peaks.bed",
#'   gene_list = expressed_genes)
#' }
epiRomics_filter_accessible <- function(putative_enhancers,
                                         epiRomics_track_connection = NULL,
                                         mode = "signal",
                                         scope = "filter_distal",
                                         signal_threshold = 2,
                                         bed_path = NULL,
                                         gene_list = NULL,
                                         promoter_distance = 2000L) {

  # --- Validation ---
  if (!base::is.data.frame(putative_enhancers) ||
      base::nrow(putative_enhancers) == 0) {
    base::stop("putative_enhancers must be a non-empty data.frame")
  }
  mode <- base::match.arg(mode, base::c("signal", "bed", "genelist", "combined"))
  scope <- base::match.arg(scope, base::c("filter_distal", "filter_all"))

  required_cols <- base::c("chr", "start", "end")
  missing_cols <- required_cols[!required_cols %in% base::names(putative_enhancers)]
  if (base::length(missing_cols) > 0) {
    base::stop("putative_enhancers missing required columns: ",
      base::paste(missing_cols, collapse = ", "))
  }

  n_pe <- base::nrow(putative_enhancers)
  pe_gr <- GenomicRanges::GRanges(
    seqnames = putative_enhancers$chr,
    ranges = IRanges::IRanges(
      start = putative_enhancers$start,
      end = putative_enhancers$end))

  # --- Scope: identify promoter-proximal regions ---
  is_promoter <- base::rep(FALSE, n_pe)
  if (scope == "filter_distal" &&
      "distanceToTSS" %in% base::names(putative_enhancers)) {
    is_promoter <- base::abs(putative_enhancers$distanceToTSS) <=
      promoter_distance
  }

  # --- Mode: signal ---
  signal_accessible <- base::rep(FALSE, n_pe)
  if (mode %in% base::c("signal", "combined")) {
    if (base::is.null(epiRomics_track_connection) ||
        !base::is.data.frame(epiRomics_track_connection)) {
      if (mode == "signal") {
        base::stop("epiRomics_track_connection required for mode='signal'")
      }
    } else {
      tc_access <- epiRomics_track_connection[
        epiRomics_track_connection$type %in% base::c("atac", "dnase"), ,
        drop = FALSE]

      if (base::nrow(tc_access) > 0) {
        base::message(base::sprintf(
          "Screening %d regions against %d ATAC/DNase track(s)...",
          n_pe, base::nrow(tc_access)))

        for (i in base::seq_len(base::nrow(tc_access))) {
          bw_path <- tc_access[i, 1]
          sample_name <- tc_access[i, 2]
          sample_type <- tc_access[i, 4]
          col_safe <- base::gsub("[^a-zA-Z0-9_]", "_",
            base::tolower(sample_name))

          base::message(base::sprintf(
            "  Importing %s (%s)...", sample_name, sample_type))

          signal_vals <- base::tryCatch({
            .fast_bw_signal(bw_path, pe_gr, type = "mean")
          }, error = function(e) {
            base::message(base::sprintf(
              "    Warning: failed to import %s - %s", bw_path, e$message))
            base::rep(0, n_pe)
          })

          putative_enhancers[[base::paste0("atac_", col_safe)]] <- signal_vals

          sig_mean <- base::mean(signal_vals)
          sig_sd <- stats::sd(signal_vals)
          threshold_val <- sig_mean + signal_threshold * sig_sd
          is_accessible <- signal_vals > threshold_val

          putative_enhancers[[base::paste0("accessible_", col_safe)]] <-
            is_accessible
          signal_accessible <- signal_accessible | is_accessible

          n_acc <- base::sum(is_accessible)
          base::message(base::sprintf(
            "    %s: %d / %d accessible (%s%%, threshold=%s)",
            sample_name, n_acc, n_pe,
            base::round(n_acc / n_pe * 100, 1),
            base::round(threshold_val, 2)))
        }
      } else if (mode == "signal") {
        base::message("No ATAC/DNase tracks found. Returning input unchanged.")
        putative_enhancers$atac_accessible <- NA
        return(putative_enhancers)
      }
    }
  }

  # --- Mode: bed ---
  bed_accessible <- base::rep(FALSE, n_pe)
  if (mode %in% base::c("bed", "combined")) {
    if (!base::is.null(bed_path)) {
      if (!base::file.exists(bed_path)) {
        base::stop("BED file not found: ", bed_path)
      }
      bed_gr <- rtracklayer::import(bed_path)
      hits <- GenomicRanges::findOverlaps(pe_gr, bed_gr)
      bed_accessible[S4Vectors::queryHits(hits)] <- TRUE
      n_bed <- base::sum(bed_accessible)
      base::message(base::sprintf(
        "BED overlap: %d / %d regions overlap BED peaks (%s%%)",
        n_bed, n_pe, base::round(n_bed / n_pe * 100, 1)))
    } else if (mode == "bed") {
      base::stop("bed_path required for mode='bed'")
    }
  }

  # --- Mode: genelist ---
  gene_accessible <- base::rep(FALSE, n_pe)
  if (mode %in% base::c("genelist", "combined")) {
    if (!base::is.null(gene_list)) {
      if (!"SYMBOL" %in% base::names(putative_enhancers)) {
        base::stop("putative_enhancers must have a 'SYMBOL' column for ",
          "mode='genelist'")
      }
      gene_accessible <- putative_enhancers$SYMBOL %in% gene_list
      n_gene <- base::sum(gene_accessible)
      base::message(base::sprintf(
        "Gene list: %d / %d regions linked to expressed genes (%s%%)",
        n_gene, n_pe, base::round(n_gene / n_pe * 100, 1)))
    } else if (mode == "genelist") {
      base::stop("gene_list required for mode='genelist'")
    }
  }

  # --- Combine evidence ---
  if (mode == "signal") {
    final_accessible <- signal_accessible
  } else if (mode == "bed") {
    final_accessible <- bed_accessible
  } else if (mode == "genelist") {
    final_accessible <- gene_accessible
  } else {
    # combined: union of all available evidence
    final_accessible <- signal_accessible | bed_accessible | gene_accessible
  }

  # Apply scope: promoter-proximal always retained
  if (scope == "filter_distal") {
    final_accessible <- final_accessible | is_promoter
  }

  putative_enhancers$atac_accessible <- final_accessible

  n_pass <- base::sum(final_accessible)
  base::message(base::sprintf(
    "Accessibility filter (mode=%s, scope=%s): %d / %d regions pass (%s%%)",
    mode, scope, n_pass, n_pe, base::round(n_pass / n_pe * 100, 1)))

  putative_enhancers
}
