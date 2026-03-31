#' Define regions of interest and filter enhanceosome by overlap
#'
#' Evaluates whether regions of interest derived from external experiments
#' (ATAC-seq, DBA, gene lists, BED files) correspond with enhanceosome regions.
#' Supports multiple input types for flexible region definition.
#'
#' @param epiRomics_putative_enhanceosome epiRomics class database containing
#'   putative enhanceosome calls. Must have non-empty \code{@@annotations}.
#' @param epiRomics_test_regions GRanges or NULL. Direct GRanges regions of
#'   interest (original interface). If provided, \code{input_type} is ignored
#'   and this is used directly for overlap. Default NULL.
#' @param input_type Character. How to construct test regions when
#'   \code{epiRomics_test_regions} is NULL. One of:
#'   \describe{
#'     \item{\code{"granges"}}{Use \code{epiRomics_test_regions} directly
#'       (default, backward-compatible).}
#'     \item{\code{"bed"}}{Import regions from a BED file via
#'       \code{bed_path}.}
#'     \item{\code{"genelist"}}{Select enhanceosome regions whose SYMBOL
#'       column matches genes in \code{gene_list}.}
#'     \item{\code{"combined"}}{Union of all provided evidence sources
#'       (GRanges, BED, genelist).}
#'   }
#' @param bed_path Character or NULL. Path to a BED file for
#'   \code{input_type = "bed"} or \code{"combined"}.
#' @param gene_list Character vector or NULL. Gene symbols for
#'   \code{input_type = "genelist"} or \code{"combined"}.
#' @return Variable of class epiRomics with enhanceosome regions overlapping
#'   with regions of interest.
#' @export
#' @examples
#' db <- epiRomicsS4(annotations = GenomicRanges::GRanges(),
#'   meta = data.frame(name = character(), type = character(),
#'     file = character(), stringsAsFactors = FALSE),
#'   genome = "hg38")
#' tryCatch(
#'   epiRomics_regions_of_interest(db),
#'   error = function(e) message(e$message))
#' \donttest{
#' # Original GRanges interface (backward-compatible)
#' roi <- epiRomics_regions_of_interest(dB, test_gr)
#'
#' # BED file input
#' roi <- epiRomics_regions_of_interest(dB, input_type = "bed",
#'   bed_path = "peaks.bed")
#'
#' # Gene list input
#' roi <- epiRomics_regions_of_interest(dB, input_type = "genelist",
#'   gene_list = c("INS", "GCG", "SST"))
#'
#' # Combined: union of BED + genelist
#' roi <- epiRomics_regions_of_interest(dB, input_type = "combined",
#'   bed_path = "peaks.bed", gene_list = c("INS", "GCG"))
#' }
epiRomics_regions_of_interest <- function(epiRomics_putative_enhanceosome,
                                           epiRomics_test_regions = NULL,
                                           input_type = c("granges", "bed",
                                             "genelist", "combined"),
                                           bed_path = NULL,
                                           gene_list = NULL) {

  # --- Core validation ---
  validate_epiRomics_dB(epiRomics_putative_enhanceosome,
    "epiRomics_regions_of_interest")

  if (base::length(epiRomics_putative_enhanceosome@annotations) == 0) {
    base::stop("epiRomics_putative_enhanceosome has no annotations")
  }

  annotations <- epiRomics_putative_enhanceosome@annotations

  # --- Backward-compatible: if GRanges passed directly, use it ---
  if (!base::is.null(epiRomics_test_regions)) {
    validate_genomic_ranges(epiRomics_test_regions,
      "epiRomics_regions_of_interest")
    return(.roi_filter(epiRomics_putative_enhanceosome,
      epiRomics_test_regions))
  }

  # --- New multi-mode logic ---
  input_type <- base::match.arg(input_type)

  if (input_type == "granges") {
    base::stop(
      "epiRomics_test_regions (GRanges) required for input_type='granges'"
    )
  }

  if (input_type == "bed") {
    if (base::is.null(bed_path)) {
      base::stop("bed_path required for input_type='bed'")
    }
    test_gr <- .roi_from_bed(bed_path)
    return(.roi_filter(epiRomics_putative_enhanceosome, test_gr))
  }

  if (input_type == "genelist") {
    if (base::is.null(gene_list) || base::length(gene_list) == 0) {
      base::stop("gene_list required for input_type='genelist'")
    }
    return(.roi_from_genelist(epiRomics_putative_enhanceosome, gene_list))
  }

  if (input_type == "combined") {
    return(.roi_combined_filter(epiRomics_putative_enhanceosome, annotations,
      bed_path, gene_list))
  }
}

#' Combined-mode ROI filter: union of BED and genelist evidence
#'
#' Applies BED overlap and/or genelist matching to the enhanceosome
#' annotations. Any region matching at least one evidence source is retained.
#'
#' @param epiRomics_putative_enhanceosome epiRomicsS4 object.
#' @param annotations GRanges of enhanceosome annotations.
#' @param bed_path Character or NULL. Path to a BED file.
#' @param gene_list Character vector or NULL. Gene symbols.
#' @return Filtered epiRomicsS4 object with combined ROI.
#' @noRd
.roi_combined_filter <- function(epiRomics_putative_enhanceosome, annotations,
                                  bed_path, gene_list) {
  keep <- base::rep(FALSE, base::length(annotations))

  # BED evidence
  if (!base::is.null(bed_path)) {
    bed_gr <- .roi_from_bed(bed_path)
    hits <- GenomicRanges::findOverlaps(annotations, bed_gr)
    keep[S4Vectors::queryHits(hits)] <- TRUE
  }

  # Genelist evidence
  if (!base::is.null(gene_list) && base::length(gene_list) > 0) {
    if ("SYMBOL" %in% base::names(S4Vectors::mcols(annotations))) {
      keep <- keep | (annotations$SYMBOL %in% gene_list)
    }
  }

  if (!base::any(keep)) {
    base::warning("No overlapping regions found for any combined evidence")
  }

  result <- epiRomics_putative_enhanceosome
  result@annotations <- annotations[keep]
  base::message(base::sprintf("Combined ROI: %d of %d regions retained",
    base::sum(keep), base::length(annotations)))
  base::return(result)
}

#' Filter enhanceosome annotations by GRanges overlap
#' @param dB epiRomicsS4 object
#' @param test_gr GRanges of test regions
#' @return Filtered epiRomicsS4 object
#' @noRd
.roi_filter <- function(dB, test_gr) {
  result <- dB
  result@annotations <- IRanges::subsetByOverlaps(dB@annotations, test_gr)

  if (base::length(result@annotations) == 0) {
    base::warning(
      "No overlapping regions found between enhanceosome and test regions"
    )
  }

  n_kept <- base::length(result@annotations)
  n_total <- base::length(dB@annotations)
  base::message(base::sprintf("ROI filter: %d of %d regions retained",
    n_kept, n_total))
  return(result)
}

#' Import GRanges from a BED file
#' @param bed_path Path to BED file
#' @return GRanges
#' @noRd
.roi_from_bed <- function(bed_path) {
  if (!base::file.exists(bed_path)) {
    base::stop("BED file not found: ", bed_path)
  }
  rtracklayer::import(bed_path)
}

#' Filter enhanceosome by gene list matching SYMBOL column
#' @param dB epiRomicsS4 object
#' @param gene_list Character vector of gene symbols
#' @return Filtered epiRomicsS4 object
#' @noRd
.roi_from_genelist <- function(dB, gene_list) {
  annotations <- dB@annotations

  if (!"SYMBOL" %in% base::names(S4Vectors::mcols(annotations))) {
    base::stop("Annotations must have a SYMBOL column for genelist mode")
  }

  keep <- annotations$SYMBOL %in% gene_list
  result <- dB
  result@annotations <- annotations[keep]

  if (!base::any(keep)) {
    base::warning("No annotations matched any gene in gene_list")
  }

  n_matched <- base::sum(keep)
  base::message(base::sprintf("Genelist ROI: %d of %d regions matched (%s)",
    n_matched, base::length(annotations),
    base::paste(base::intersect(gene_list, annotations$SYMBOL),
      collapse = ", ")))
  return(result)
}
