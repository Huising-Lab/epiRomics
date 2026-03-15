#' Gviz-based multi-track genomic visualization (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function has been superseded by \code{\link{epiRomics_track_layer}},
#' which uses base R graphics for ~20x faster rendering. This Gviz-based
#' version is retained for users who need IdeogramTrack rendering or
#' Gviz-specific customization.
#'
#' Generates a multi-track genomic visualization for a putative enhanceosome
#' region, including gene model, chromatin accessibility (ATAC), RNA-seq signal,
#' enhancer annotation, and ChIP-seq peak tracks. Genome is automatically
#' detected from the database object.
#'
#' @inheritParams epiRomics_track_layer
#' @return GViz plot
#' @export
#' @examples
#' db <- epiRomicsS4(annotations = GenomicRanges::GRanges(),
#'   meta = data.frame(name = character(), type = character(),
#'     file = character(), stringsAsFactors = FALSE),
#'   genome = "hg38")
#' tryCatch(epiRomics_track_layer_deprecated(db, 1, db,
#'   data.frame(file = character(), name = character(),
#'     color = character(), type = character())),
#'   error = function(e) message(e$message))
#' \donttest{
#' epiRomics_track_layer_deprecated(enhanceosome, 1, epiRomics_dB,
#'   epiRomics_track_connection)
#' }
epiRomics_track_layer_deprecated <- function(epiRomics_putative_enhanceosome, epiRomics_index, epiRomics_dB, epiRomics_track_connection, epiRomics_keep_epitracks = TRUE, chromatin_states = NULL) {
  if (!base::requireNamespace("Gviz", quietly = TRUE)) {
    base::stop("Package 'Gviz' is required for this deprecated function. ",
               "Install it with: BiocManager::install('Gviz')\n",
               "Or use epiRomics_track_layer() instead (no Gviz dependency).",
               call. = FALSE)
  }
  validate_epiRomics_dB(epiRomics_putative_enhanceosome, "epiRomics_track_layer_deprecated")
  validate_numeric_param(epiRomics_index, "epiRomics_index", "epiRomics_track_layer_deprecated", min_val = 1)
  validate_epiRomics_dB(epiRomics_dB, "epiRomics_track_layer_deprecated")
  validate_logical_param(epiRomics_keep_epitracks, "epiRomics_keep_epitracks", "epiRomics_track_layer_deprecated")

  if (!base::is.data.frame(epiRomics_track_connection)) {
    base::stop("epiRomics_track_layer_deprecated: epiRomics_track_connection must be a data frame")
  }
  if (base::nrow(epiRomics_track_connection) == 0) {
    base::stop("epiRomics_track_layer_deprecated: epiRomics_track_connection is empty")
  }
  if (epiRomics_index > base::length(epiRomics_putative_enhanceosome@annotations)) {
    base::stop(base::sprintf(
      "epiRomics_track_layer_deprecated: epiRomics_index (%d) exceeds number of annotations (%d)",
      epiRomics_index, base::length(epiRomics_putative_enhanceosome@annotations)
    ))
  }
  tryCatch(
    {
      txdb_obj <- resolve_txdb(epiRomics_dB@txdb)

      genome <- epiRomics_dB@genome

      epiRomics_class_save <- epiRomics_putative_enhanceosome
      epiRomics_putative_enhanceosome <- epiRomics_putative_enhanceosome@annotations[epiRomics_index, ]
      chr <- base::as.character(BiocGenerics::unique(GenomeInfoDb::seqnames(epiRomics_putative_enhanceosome)))
      gtrack <- Gviz::GenomeAxisTrack()
      txTr <- Gviz::GeneRegionTrack(txdb_obj,
        start = epiRomics_putative_enhanceosome$geneStart,
        end = epiRomics_putative_enhanceosome$geneEnd,
        chromosome = epiRomics_putative_enhanceosome$geneChr,
        symbol = epiRomics_putative_enhanceosome$SYMBOL,
        name = epiRomics_putative_enhanceosome$SYMBOL,
        transcriptAnnotation = "symbol",
        geneSymbol = TRUE, howId = TRUE, min.height = 10, size = 10
      )

      # Determine chromatin state label for the enhancer index
      enhancer_fill <- "#7B2D8E"  # strong violet --highlights the region of interest
      enhancer_state_label <- "Putative Enhancer"
      if (!base::is.null(chromatin_states) && base::is.data.frame(chromatin_states)) {
        enh_start <- BiocGenerics::start(epiRomics_putative_enhanceosome)
        enh_end <- BiocGenerics::end(epiRomics_putative_enhanceosome)
        enh_chr <- base::as.character(GenomeInfoDb::seqnames(epiRomics_putative_enhanceosome))
        # Find overlapping chromatin state for label only
        match_idx <- base::which(
          chromatin_states$seqnames == enh_chr &
          chromatin_states$start <= enh_end &
          chromatin_states$end >= enh_start
        )
        if (base::length(match_idx) > 0) {
          state <- chromatin_states$chromatin_state[match_idx[1]]
          enhancer_state_label <- base::gsub("_", " ", state)
        }
      }

      atrack <- Gviz::AnnotationTrack(epiRomics_putative_enhanceosome,
        name = base::paste0("Enhancer (", enhancer_state_label, ")"),
        fill = enhancer_fill,
        col = enhancer_fill,
        min.height = 15, size = 15
      )
      itrack <- get_ideogram_track(genome, chr)

      if (base::min(atrack) < base::min(txTr)) {
        min_track <- base::min(atrack)
      } else {
        min_track <- base::min(txTr)
      }
      if (base::max(atrack) < base::max(txTr)) {
        max_track <- base::max(txTr)
      } else {
        max_track <- base::max(atrack)
      }
      min_track <- min_track - 1000
      max_track <- max_track + 1000
      region_track <- GenomicRanges::GRanges(seqnames = chr, ranges = IRanges::IRanges(start = min_track, end = max_track))

      # Build ChIP/TF, histone variant, and ncRNA annotation tracks
      # Color scheme: dull blue for background peaks, strong violet for peaks
      # overlapping the enhancer index (region of interest)
      epiRomics_db_data_table <- epiRomics_class_save@meta
      dull_blue <- "#7F8FA6"
      highlight_violet <- "#7B2D8E"
      enh_region <- epiRomics_putative_enhanceosome  # the single enhancer index

      # Pre-filter ALL annotations to viewing region ONCE (major speed boost)
      # This avoids repeated subsetByOverlaps on the full annotation database
      region_annots <- IRanges::subsetByOverlaps(epiRomics_dB@annotations, region_track)
      region_types <- region_annots$type

      # Helper: build colored annotation tracks from pre-filtered data
      .build_annot_tracks <- function(names_vec, db_prefix, base_color = dull_blue,
                                       min_h = 10, sz = 10,
                                       display_transform = base::identity) {
        tracks_out <- base::list()
        for (i in base::seq_along(names_vec)) {
          db_key <- base::paste0(db_prefix, names_vec[i])
          annot_data <- region_annots[region_types == db_key]
          if (base::length(annot_data) > 0) {
            overlaps_enh <- GenomicRanges::countOverlaps(annot_data, enh_region) > 0
            peak_colors <- base::ifelse(overlaps_enh, highlight_violet, base_color)
            tracks_out <- base::c(tracks_out, Gviz::AnnotationTrack(
              annot_data,
              name = display_transform(names_vec[i]),
              fill = peak_colors,
              col = peak_colors,
              min.height = min_h, size = sz
            ))
          }
        }
        tracks_out
      }

      db_prefix <- base::paste0(genome, "_custom_")

      # TF tracks (type == "chip") --dull blue base
      epiRomics_chips <- epiRomics_db_data_table[epiRomics_db_data_table$type == "chip", "name"]
      chip_tracks <- .build_annot_tracks(epiRomics_chips, db_prefix)

      # Histone variant peaks (e.g., H2A.Z) are NOT shown as separate tracks.
      # Their information is captured through the chromatin state classification
      # (e.g., active enhancer, poised enhancer) which integrates all histone
      # marks. Showing individual histone peaks would be redundant with the
      # chromatin state tracks.

      # ncRNA tracks (auto-detected from database) --dull blue base
      ncrna_meta <- epiRomics_db_data_table[epiRomics_db_data_table$type == "ncrna", , drop = FALSE]
      ncrna_display <- function(nm) base::toupper(base::gsub("_", " ", nm))
      ncrna_tracks <- .build_annot_tracks(
        ncrna_meta$name, db_prefix, min_h = 8, sz = 8,
        display_transform = ncrna_display
      )

      # Build per-state chromatin annotation tracks if states provided
      state_tracks <- base::list()
      if (!base::is.null(chromatin_states) && base::is.data.frame(chromatin_states) &&
          base::nrow(chromatin_states) > 0) {
        region_states <- chromatin_states[
          chromatin_states$seqnames == chr &
          chromatin_states$end >= min_track &
          chromatin_states$start <= max_track, , drop = FALSE
        ]
        if (base::nrow(region_states) > 0) {
          palette <- .chromatin_state_palette()
          unique_states <- base::unique(region_states$chromatin_state)
          for (st in unique_states) {
            st_rows <- region_states[region_states$chromatin_state == st, , drop = FALSE]
            st_gr <- GenomicRanges::GRanges(
              seqnames = st_rows$seqnames,
              ranges = IRanges::IRanges(start = st_rows$start, end = st_rows$end)
            )
            st_color <- if (st %in% base::names(palette)) palette[st] else "#AAAAAA"
            st_label <- base::gsub("_", " ", st)
            state_tracks <- base::c(state_tracks, Gviz::AnnotationTrack(
              st_gr,
              name = st_label,
              fill = st_color,
              col = st_color,
              min.height = 8, size = 8
            ))
          }
        }
      }

      # Pre-import ATAC BigWig data and build tracks (single-pass I/O)
      epiRomics_track_connection_chromatin <- epiRomics_track_connection[epiRomics_track_connection$type == "atac", ]
      atac_result <- .build_data_tracks(
        bw_paths = epiRomics_track_connection_chromatin[, 1],
        region = region_track,
        sample_names = epiRomics_track_connection_chromatin[, 2],
        colors = epiRomics_track_connection_chromatin[, 3],
        genome = genome,
        chr = chr
      )
      chromatin_tracks <- atac_result$tracks

      # Pre-import RNA BigWig data and build tracks (single-pass I/O)
      rna_tracks <- base::list()
      if (base::any(epiRomics_track_connection$type == "rna")) {
        epiRomics_track_connection_rna <- epiRomics_track_connection[epiRomics_track_connection$type == "rna", ]
        rna_result <- .build_data_tracks(
          bw_paths = epiRomics_track_connection_rna[, 1],
          region = region_track,
          sample_names = epiRomics_track_connection_rna[, 2],
          colors = epiRomics_track_connection_rna[, 3],
          genome = genome,
          chr = chr
        )
        rna_tracks <- rna_result$tracks
      }
      # Pre-import histone BigWig data and build tracks (single-pass I/O)
      # Supports type = "histone" in epiRomics_track_connection CSV
      histone_signal_tracks <- base::list()
      if (base::any(epiRomics_track_connection$type == "histone")) {
        tc_histone <- epiRomics_track_connection[
          epiRomics_track_connection$type == "histone", , drop = FALSE]
        histone_result <- .build_data_tracks(
          bw_paths = tc_histone[, 1],
          region = region_track,
          sample_names = tc_histone[, 2],
          colors = tc_histone[, 3],
          genome = genome,
          chr = chr
        )
        histone_signal_tracks <- histone_result$tracks
      }

      # Wrap data tracks in a HighlightTrack to draw a vertical shaded
      # column over the enhancer index region of interest.
      # Track order: signal data -> enhancer index -> chromatin context ->
      #   TF binding -> ncRNA
      highlight_region <- epiRomics_putative_enhanceosome
      all_data_tracks <- base::c(
        chromatin_tracks,          # ATAC signal (continuous)
        rna_tracks,                # RNA signal (continuous)
        histone_signal_tracks,     # Histone mark signal (continuous)
        atrack,                    # Enhancer index (region of interest)
        state_tracks,              # Chromatin state annotations
        chip_tracks,               # TF binding tracks
        ncrna_tracks               # ncRNA annotation tracks
      )
      ht <- Gviz::HighlightTrack(
        trackList = all_data_tracks,
        range = highlight_region,
        col = "#D8BFD8",    # light violet border
        fill = "#F0E6F6",   # very light violet background
        alpha = 0.3,
        inBackground = TRUE
      )

      if (epiRomics_keep_epitracks) {
        Gviz::plotTracks(base::c(itrack, gtrack, txTr, ht), from = min_track, to = max_track)
      } else {
        Gviz::plotTracks(base::c(itrack, gtrack, txTr, chromatin_tracks, rna_tracks), from = min_track, to = max_track)
      }
    },
    error = function(e) {
      base::stop(base::sprintf("epiRomics_track_layer_deprecated: %s", e$message))
    }
  )
}
