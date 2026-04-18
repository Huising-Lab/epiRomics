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
    n_cores <- .detect_cores(max_cores = n_s)
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

#' Collect histone-based putative enhancer regions
#'
#' Auto-computes chromatin states if NULL, selects enhancer-related rows
#' including H2A.Z recovery, builds histone GRanges from classified regions.
#'
#' @param chromatin_states data.frame or NULL. Pre-computed chromatin states.
#' @param database An epiRomics S4 database object.
#' @param enhancer_states Character vector of enhancer state names.
#' @return list with gr (GRanges or NULL) and labels (character or NULL)
#'   and chromatin_states (data.frame, possibly auto-computed)
#' @noRd
.collect_histone_source <- function(chromatin_states, database,
                                     enhancer_states) {
  if (base::is.null(chromatin_states)) {
    base::message("Auto-computing chromatin states from all histone marks...")
    chromatin_states <- classify_chromatin_states(database)
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
      fmt <- base::paste0(
        "H2A.Z: recovered %d additional ",
        "regions as putative enhancers ",
        "(H2A.Z-positive, not classified ",
        "as enhancer by standard histone ",
        "rules)"
      )
      base::message(base::sprintf(fmt, n_h2az_extra))
    }
  }

  enh_data <- chromatin_states[enh_rows, , drop = FALSE]

  if (base::nrow(enh_data) > 0) {
    histone_gr <- GenomicRanges::GRanges(
      seqnames = enh_data$seqnames,
      ranges = IRanges::IRanges(start = enh_data$start, end = enh_data$end))
    histone_gr <- GenomicRanges::reduce(histone_gr)

    base::message(base::sprintf(
      "Histone-based enhancers: %d regions (from %d classified regions)",
      base::length(histone_gr), base::nrow(enh_data)))
    # State breakdown
    state_tab <- base::table(enh_data$chromatin_state)
    base::message(base::sprintf("  States: %s", base::paste(
      base::paste0(base::names(state_tab), "=", state_tab),
      collapse = ", ")))

    base::return(base::list(gr = histone_gr,
      labels = base::rep("histone", base::length(histone_gr)),
      chromatin_states = chromatin_states))
  }

  base::message("No enhancer-like regions found from histone marks")
  base::list(gr = NULL, labels = NULL, chromatin_states = chromatin_states)
}

#' Collect Hi-C contact anchor regions
#'
#' Builds GRanges from Hi-C contact anchors for putative enhancer sources.
#'
#' @param hic_contacts data.frame or NULL. Hi-C contacts in BEDPE format.
#' @return list with gr (GRanges or NULL) and labels (character or NULL)
#' @noRd
.collect_hic_source <- function(hic_contacts) {
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

    base::message(base::sprintf(
      "Hi-C source: %d unique anchor regions from %d contacts",
      base::length(hic_anchors), base::nrow(hic_contacts)))

    base::return(base::list(gr = hic_anchors,
      labels = base::rep("hic", base::length(hic_anchors))))
  }

  base::list(gr = NULL, labels = NULL)
}

#' Collect TF binding regions excluding H2A.Z
#'
#' Collects TF binding regions from ChIP-seq peaks, excluding H2A.Z.
#'
#' @param database An epiRomics S4 database object.
#' @param annot_by_type Named list of GRanges, split by annotation type.
#' @param db_prefix Character. Genome-specific prefix for annotation keys.
#' @return list with gr (GRanges or NULL) and labels (character or NULL)
#' @noRd
.collect_tf_source <- function(database, annot_by_type, db_prefix) {
  tf_meta <- database@meta[database@meta$type == "chip", , drop = FALSE]
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
      base::message(base::sprintf(
        "TF binding source: %d regions from %d TFs (%s)",
        base::length(tf_reduced), base::length(tf_chip_names),
        base::paste(tf_chip_names, collapse = ", ")))
      base::return(base::list(gr = tf_reduced,
        labels = base::rep("tf", base::length(tf_reduced))))
    }
  }

  base::list(gr = NULL, labels = NULL)
}

#' Union and annotate putative enhancer source regions
#'
#' Concatenates all GRanges sources, reduces overlaps, and builds
#' source attribution strings per merged region.
#'
#' @param source_gr_list Named list of GRanges from each source.
#' @param source_label_list Named list of character label vectors.
#' @return list with merged (GRanges), merged_sources (character),
#'   all_gr (concatenated GRanges), all_labels (character)
#' @noRd
.union_and_annotate_sources <- function(source_gr_list, source_label_list) {
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

  # Source attribution (vectorized)
  source_sets <- base::split(base::seq_along(all_gr), all_labels)
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

  base::list(merged = merged, merged_sources = merged_sources,
    all_gr = all_gr, all_labels = all_labels)
}

#' Annotate merged regions with chromatin state
#'
#' Vectorized chromatin state annotation using findOverlaps and
#' data.table grouping to assign broad and detail states.
#'
#' @param merged GRanges of merged putative enhancer regions.
#' @param chromatin_states data.frame of chromatin state classifications.
#' @param enhancer_states Character vector of enhancer state names.
#' @return list with broad_col (character) and detail_col (character)
#' @noRd
.annotate_chromatin_state <- function(merged, chromatin_states,
                                       enhancer_states) {
  n_merged <- base::length(merged)

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

  base::list(broad_col = broad_col, detail_col = detail_col)
}

#' Annotate merged regions with histone mark overlaps
#'
#' Builds histone overlap matrix, computes overlap name strings,
#' counts, and H2A.Z detection.
#'
#' @param merged GRanges of merged putative enhancer regions.
#' @param database An epiRomics S4 database object.
#' @param annot_by_type Named list of GRanges, split by annotation type.
#' @param db_prefix Character. Genome-specific prefix for annotation keys.
#' @return list with hist_overlap_names (character), hist_overlap_counts
#'   (integer), h2az_col (logical)
#' @noRd
.annotate_histone_marks <- function(merged, database, annot_by_type,
                                     db_prefix) {
  n_merged <- base::length(merged)
  hist_names <- database@meta$name[database@meta$type == "histone"]

  # Build boolean matrix: merged[i] overlaps histone[j]?
  # (parallel when available)
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

  base::list(hist_overlap_names = hist_overlap_names,
    hist_overlap_counts = hist_overlap_counts, h2az_col = h2az_col)
}

#' Annotate merged regions with TF co-binding overlaps
#'
#' Builds TF overlap matrix, computes TF overlap name strings and counts.
#' H2A.Z is excluded from TF counting.
#'
#' @param merged GRanges of merged putative enhancer regions.
#' @param database An epiRomics S4 database object.
#' @param annot_by_type Named list of GRanges, split by annotation type.
#' @param db_prefix Character. Genome-specific prefix for annotation keys.
#' @return list with tf_overlap_names (character) and tf_overlap_counts
#'   (integer)
#' @noRd
.annotate_tf_cobinding <- function(merged, database, annot_by_type,
                                    db_prefix) {
  n_merged <- base::length(merged)
  tf_names <- database@meta$name[database@meta$type == "chip"]
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

  base::list(tf_overlap_names = tf_overlap_names,
    tf_overlap_counts = tf_overlap_counts)
}

#' Build putative enhancer result data.frame
#'
#' Assembles the final output data.frame from merged regions and all
#' annotation columns, sorts by state priority and co-binding, and
#' prints summary messages.
#'
#' @param merged GRanges of merged putative enhancer regions.
#' @param merged_sources Character vector of source attribution strings.
#' @param broad_col Character vector of broad chromatin state labels.
#' @param detail_col Character vector of detail chromatin state labels.
#' @param hist_overlap_names Character vector of histone mark name strings.
#' @param hist_overlap_counts Integer vector of histone mark counts.
#' @param h2az_col Logical vector of H2A.Z overlap flags.
#' @param tf_overlap_names Character vector of TF name strings.
#' @param tf_overlap_counts Integer vector of TF counts.
#' @return data.frame with putative enhancer annotations, sorted.
#' @noRd
.build_putative_result <- function(merged, merged_sources, broad_col,
                                    detail_col, hist_overlap_names,
                                    hist_overlap_counts, h2az_col,
                                    tf_overlap_names, tf_overlap_counts) {
  n_merged <- base::length(merged)

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

#' Identify putative enhancer regions using rule-based histone logic
#'
#' Automatically scans all histone and histone variant marks in the epiRomics
#' database and applies ChromHMM-based classification rules to identify putative
#' enhancer regions.
#'
#' This function uses a multi-source approach:
#' \describe{
#'   \item{Chromatin states}{Leverages \code{\link{classify_chromatin_states}}
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
#' @param database An epiRomics S4 database object.
#' @param chromatin_states data.frame or NULL. Pre-computed output from
#'   \code{\link{classify_chromatin_states}}. If NULL (default), computed
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
#'       \code{\link{classify_chromatin_states}} (e.g. active_enhancer,
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
#' db <- make_example_database()
#' pe <- find_putative_enhancers(db)
#' head(pe[, c("chr", "start", "end", "chromatin_state")])
find_putative_enhancers <- function(database,
                                          chromatin_states = NULL,
                                          hic_contacts = NULL,
                                          enhancer_states = base::c(
                                            "active",
                                            "poised",
                                            "primed")) {

  validate_database(database, "find_putative_enhancers")

  genome <- database@genome
  db_prefix <- base::paste0(genome, "_custom_")
  source_gr_list <- base::list()
  source_label_list <- base::list()

  # Pre-split annotations by type for O(1) lookup
  annot_by_type <- base::split(
    database@annotations,
    database@annotations$type
  )

  # SOURCE 1: Rule-based histone classification (automatic)
  hist_result <- .collect_histone_source(chromatin_states, database,
    enhancer_states)
  chromatin_states <- hist_result$chromatin_states
  if (!base::is.null(hist_result$gr)) {
    source_gr_list[["histone"]] <- hist_result$gr
    source_label_list[["histone"]] <- hist_result$labels
  }

  # SOURCE 2: Hi-C contact anchors
  hic_result <- .collect_hic_source(hic_contacts)
  if (!base::is.null(hic_result$gr)) {
    source_gr_list[["hic"]] <- hic_result$gr
    source_label_list[["hic"]] <- hic_result$labels
  }

  # SOURCE 3: TF binding regions (ChIP-seq peaks, excluding H2A.Z)
  tf_result <- .collect_tf_source(database, annot_by_type, db_prefix)
  if (!base::is.null(tf_result$gr)) {
    source_gr_list[["tf"]] <- tf_result$gr
    source_label_list[["tf"]] <- tf_result$labels
  }

  # UNION all regions and annotate sources
  union_result <- .union_and_annotate_sources(source_gr_list,
    source_label_list)
  merged <- union_result$merged
  merged_sources <- union_result$merged_sources

  # ANNOTATE: chromatin state
  cs_annot <- .annotate_chromatin_state(merged, chromatin_states,
    enhancer_states)

  # ANNOTATE: histone marks + H2A.Z
  hist_annot <- .annotate_histone_marks(merged, database, annot_by_type,
    db_prefix)

  # ANNOTATE: TF co-binding
  tf_annot <- .annotate_tf_cobinding(merged, database, annot_by_type,
    db_prefix)

  # BUILD OUTPUT
  .build_putative_result(merged, merged_sources,
    cs_annot$broad_col, cs_annot$detail_col,
    hist_annot$hist_overlap_names, hist_annot$hist_overlap_counts,
    hist_annot$h2az_col,
    tf_annot$tf_overlap_names, tf_annot$tf_overlap_counts)
}


#' Compute functional database overlap columns
#'
#' Loops over functional databases computing boolean overlap columns
#' for each database against the putative enhancer GRanges.
#'
#' @param pe_gr GRanges of putative enhancer regions.
#' @param n_pe Integer. Number of putative enhancers.
#' @param func_meta data.frame. Functional database metadata rows.
#' @param annot_by_type Named list of GRanges, split by annotation type.
#' @param db_prefix Character. Genome-specific prefix for annotation keys.
#' @return list with overlap_cols (named list of logical vectors) and
#'   col_names (character vector of column names)
#' @noRd
.compute_database_overlaps <- function(pe_gr, n_pe, func_meta, annot_by_type,
                                        db_prefix) {
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

  base::list(overlap_cols = overlap_cols, col_names = col_names)
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
#'   \code{\link{find_putative_enhancers}}.
#' @param database An epiRomics S4 database object.
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
#' db <- make_example_database()
#' pe <- find_putative_enhancers(db)
#' pe_annot <- annotate_enhancers(pe, db)
#' table(pe_annot$novel)
annotate_enhancers <- function(putative_enhancers, database) {

  if (!base::is.data.frame(putative_enhancers) ||
      base::nrow(putative_enhancers) == 0) {
    base::stop("putative_enhancers must be a non-empty data.frame")
  }
  validate_database(database, "annotate_enhancers")

  genome <- database@genome
  db_prefix <- base::paste0(genome, "_custom_")
  meta <- database@meta

  # Pre-split annotations by type for O(1) lookup
  annot_by_type <- base::split(
    database@annotations,
    database@annotations$type
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
    base::message("No functional databases found in database")
    putative_enhancers$n_databases <- 0L
    putative_enhancers$novel <- TRUE
    return(putative_enhancers)
  }

  n_pe <- base::nrow(putative_enhancers)

  # Compute overlaps via helper
  db_result <- .compute_database_overlaps(pe_gr, n_pe, func_meta,
    annot_by_type, db_prefix)
  overlap_cols <- db_result$overlap_cols
  col_names <- db_result$col_names

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


#' Filter putative enhancers by signal-based accessibility
#'
#' Imports ATAC-seq/DNase-seq BigWig signal over each region and flags
#' regions with mean signal above a z-score threshold as accessible.
#' Modifies putative_enhancers by adding per-sample signal and
#' accessibility columns.
#'
#' @param putative_enhancers data.frame of putative enhancers.
#' @param pe_gr GRanges of putative enhancer regions.
#' @param n_pe Integer. Number of putative enhancers.
#' @param track_connection data.frame or NULL. BigWig track sheet.
#' @param mode Character. Current filtering mode.
#' @param signal_threshold Numeric. Z-score threshold.
#' @return list with putative_enhancers (updated data.frame),
#'   signal_accessible (logical vector), and early_return (data.frame or NULL)
#' @noRd
.filter_by_signal <- function(putative_enhancers, pe_gr, n_pe,
                               track_connection, mode,
                               signal_threshold) {
  signal_accessible <- base::rep(FALSE, n_pe)

  if (!mode %in% base::c("signal", "combined")) {
    base::return(base::list(putative_enhancers = putative_enhancers,
      signal_accessible = signal_accessible, early_return = NULL))
  }

  if (base::is.null(track_connection) ||
      !base::is.data.frame(track_connection)) {
    if (mode == "signal") {
      base::stop("track_connection required for mode='signal'")
    }
    base::return(base::list(putative_enhancers = putative_enhancers,
      signal_accessible = signal_accessible, early_return = NULL))
  }

  tc_access <- track_connection[
    track_connection$type %in% base::c("atac", "dnase"), ,
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
        base::warning(base::sprintf(
          "failed to import %s - %s", bw_path, e$message))
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
    base::return(base::list(putative_enhancers = putative_enhancers,
      signal_accessible = signal_accessible,
      early_return = putative_enhancers))
  }

  base::list(putative_enhancers = putative_enhancers,
    signal_accessible = signal_accessible, early_return = NULL)
}

#' Filter putative enhancers by BED overlap accessibility
#'
#' Overlaps putative enhancer GRanges with an external BED file and
#' returns a logical vector of accessible regions.
#'
#' @param pe_gr GRanges of putative enhancer regions.
#' @param n_pe Integer. Number of putative enhancers.
#' @param bed_path Character or NULL. Path to a BED file.
#' @param mode Character. Current filtering mode.
#' @return Logical vector of BED-accessible regions.
#' @noRd
.filter_by_bed <- function(pe_gr, n_pe, bed_path, mode) {
  bed_accessible <- base::rep(FALSE, n_pe)
  if (!mode %in% base::c("bed", "combined")) {
    base::return(bed_accessible)
  }

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

  bed_accessible
}

#' Filter putative enhancers by gene list accessibility
#'
#' Checks putative enhancer SYMBOL column against a gene list and
#' returns a logical vector of gene-list-accessible regions.
#'
#' @param putative_enhancers data.frame of putative enhancers.
#' @param n_pe Integer. Number of putative enhancers.
#' @param gene_list Character vector or NULL. Expressed gene symbols.
#' @param mode Character. Current filtering mode.
#' @return Logical vector of gene-list-accessible regions.
#' @noRd
.filter_by_genelist <- function(putative_enhancers, n_pe, gene_list, mode) {
  gene_accessible <- base::rep(FALSE, n_pe)
  if (!mode %in% base::c("genelist", "combined")) {
    base::return(gene_accessible)
  }

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

  gene_accessible
}

#' Combine accessibility evidence and apply scope
#'
#' Combines signal, BED, and genelist accessibility evidence, applies
#' scope (filter_distal retains promoters), adds atac_accessible column,
#' and prints summary.
#'
#' @param signal_accessible Logical vector of signal-accessible regions.
#' @param bed_accessible Logical vector of BED-accessible regions.
#' @param gene_accessible Logical vector of gene-list-accessible regions.
#' @param is_promoter Logical vector of promoter-proximal regions.
#' @param mode Character. Current filtering mode.
#' @param scope Character. Current filtering scope.
#' @param putative_enhancers data.frame of putative enhancers.
#' @param n_pe Integer. Number of putative enhancers.
#' @return Updated putative_enhancers with atac_accessible column.
#' @noRd
.combine_accessibility <- function(signal_accessible, bed_accessible,
                                    gene_accessible, is_promoter, mode,
                                    scope, putative_enhancers, n_pe) {
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

#' Filter putative enhancers by chromatin accessibility evidence
#'
#' Multi-mode accessibility filter for putative enhancers. Supports four
#' complementary evidence types that can be used independently or combined.
#'
#' @section Modes:
#' \describe{
#'   \item{signal}{Import ATAC-seq/DNase-seq BigWig signal over each region.
#'     Regions with mean signal above a z-score threshold
#'     are flagged accessible.
#'     Requires \code{track_connection} with ATAC/DNase tracks.}
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
#'   \code{\link{find_putative_enhancers}}. Must contain columns
#'   \code{chr}, \code{start}, \code{end}. For \code{genelist} mode, also
#'   requires a \code{SYMBOL} column.
#' @param track_connection data.frame or NULL. BigWig track
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
#' db <- make_example_database()
#' pe <- find_putative_enhancers(db)
#' ## Gene-list mode: attach synthetic SYMBOL column, then filter
#' pe$SYMBOL <- paste0("GENE", seq_len(nrow(pe)))
#' pe_genes <- filter_accessible_regions(
#'   pe, mode = "genelist", gene_list = c("GENE1", "GENE2")
#' )
#' sum(pe_genes$atac_accessible)
filter_accessible_regions <- function(putative_enhancers,
                                         track_connection = NULL,
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
  mode <- base::match.arg(
    mode,
    base::c("signal", "bed", "genelist", "combined"))
  scope <- base::match.arg(scope, base::c("filter_distal", "filter_all"))

  required_cols <- base::c("chr", "start", "end")
  missing_cols <- required_cols[
    !required_cols %in% base::names(putative_enhancers)]
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

  # --- Filter by each mode ---
  sig_result <- .filter_by_signal(putative_enhancers, pe_gr, n_pe,
    track_connection, mode, signal_threshold)
  putative_enhancers <- sig_result$putative_enhancers
  if (!base::is.null(sig_result$early_return)) {
    base::return(sig_result$early_return)
  }

  bed_accessible <- .filter_by_bed(pe_gr, n_pe, bed_path, mode)
  gene_accessible <- .filter_by_genelist(putative_enhancers, n_pe,
    gene_list, mode)

  # --- Combine evidence and apply scope ---
  .combine_accessibility(sig_result$signal_accessible, bed_accessible,
    gene_accessible, is_promoter, mode, scope, putative_enhancers, n_pe)
}
