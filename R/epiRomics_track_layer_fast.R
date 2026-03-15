#' Import BigWig signal data for a track group
#'
#' Imports BigWig files for a set of tracks within a genomic region.
#' Returns a list with imported data, shared max score, names, and colors.
#' Used internally by \code{epiRomics_track_layer()} to load ATAC, RNA,
#' and histone signal data.
#'
#' @param tc_sub data.frame subset of track_connection (path, name, color columns)
#' @param region_gr GRanges object defining the genomic region
#' @param xlims numeric vector of length 2 with plot x-axis limits
#' @return list with elements: data (list of data.frames), max (numeric),
#'   names (character), colors (character); or NULL if tc_sub has 0 rows
#' @noRd
.import_bw_group_data <- function(tc_sub, region_gr, xlims) {
  n <- base::nrow(tc_sub)
  if (n == 0) return(NULL)
  imported <- base::vector("list", n)
  max_scores <- base::numeric(n)
  for (j in base::seq_len(n)) {
    bw <- rtracklayer::import(tc_sub[j, 1], format = "BigWig",
      selection = rtracklayer::BigWigSelection(region_gr))
    ovlp <- IRanges::subsetByOverlaps(bw, region_gr)
    if (base::length(ovlp) > 0) {
      imported[[j]] <- base::data.frame(
        pos = (BiocGenerics::start(ovlp) + BiocGenerics::end(ovlp)) / 2,
        score = ovlp$score)
      # Use 98th percentile to cap extreme outlier peaks --standard
      # genome browser practice (IGV auto-scale). Prevents ultra-high
      # loci like INS promoter (88K ATAC) from compressing all signal.
      # Peaks above the cap are clipped (flat-topped) by base R graphics.
      if (base::length(ovlp$score) > 100) {
        max_scores[j] <- stats::quantile(ovlp$score, 0.98,
          names = FALSE)
      } else {
        max_scores[j] <- base::max(ovlp$score)
      }
    } else {
      imported[[j]] <- base::data.frame(pos = xlims, score = base::c(0, 0))
    }
  }
  shared_max <- base::ceiling((base::max(max_scores) * 1.1) / 10) * 10
  if (shared_max <= 0) shared_max <- 1
  base::list(data = imported, max = shared_max,
    names = tc_sub[, 2], colors = tc_sub[, 3])
}

#' Multi-track genomic visualization (base R graphics)
#'
#' Renders a stacked multi-panel genome browser view using base R graphics.
#' Includes gene model, BigWig signal tracks (ATAC, RNA, histone), enhancer
#' index, chromatin state bars, FANTOM/UCNE annotations,
#' TF binding peaks, and ncRNA annotations. A translucent
#' violet highlight column marks the enhancer region of interest across all
#' panels. Typically renders in 1-2 seconds.
#'
#' For the original Gviz-based version (slower, but supports IdeogramTrack),
#' see \code{\link{epiRomics_track_layer_deprecated}}.
#'
#' @param epiRomics_putative_enhanceosome epiRomics class database containing putative enhanceosome calls
#' @param epiRomics_index numeric of row value from epiRomics_putative_enhanceosome to visualize
#' @param epiRomics_dB epiRomics class database containing all data initially loaded
#' @param epiRomics_track_connection data frame containing bigwig track locations and their names
#' @param epiRomics_keep_epitracks logical indicating whether to show enhancer and chip tracks, default is TRUE
#' @param chromatin_states data.frame, optional output from
#'   \code{\link{epiRomics_chromatin_states}}. When provided, the enhancer
#'   track is colored by chromatin state and separate per-state annotation
#'   tracks are added.
#' @param gene_lookup logical. When TRUE, omits the enhancer index bar and
#'   violet highlight column. Used internally by
#'   \code{\link{epiRomics_track_layer_gene}} to display a gene locus without
#'   enhancer-specific visual elements. Default is FALSE.
#' @param show_bigwig logical. Show BigWig signal tracks. Default TRUE.
#' @param show_chromatin logical. Show chromatin state color tracks. Default TRUE.
#' @param show_annotations logical. Show BED annotation tracks. Default TRUE.
#' @param show_gene_model logical. Show gene model panel. Default TRUE.
#' @param show_enhancer_highlight logical. Show enhancer index highlight.
#'   Default TRUE.
#' @param mirror logical. Use symmetric mirrored axes for paired ATAC/RNA
#'   signals. Default TRUE.
#' @param signal_style character. Signal rendering style: \code{"line"}
#'   (default) draws vertical bars at each position (IGV/UCSC browser style),
#'   \code{"polygon"} draws filled area charts.
#' @param signal_layout character. Signal rendering layout: \code{"auto"}
#'   detects paired signals and mirrors, \code{"stacked"} renders one row per
#'   signal, \code{"mirrored"} always mirrors first two. Default \code{"auto"}.
#' @param cex_cell_label numeric. Font size for cell type labels (e.g. Alpha,
#'   Beta). Default 1.4.
#' @param cex_axis numeric. Font size for y-axis tick labels on signal tracks.
#'   Default 1.2.
#' @param cex_coord numeric. Font size for chromosome, start, stop, genome
#'   build text. Default 1.3.
#' @param cex_annotation numeric. Font size for BED annotation labels.
#'   Default 1.1.
#' @param cex_gene numeric. Font size for gene name labels in gene model.
#'   Default 1.2.
#' @param cex_title numeric. Font size for main plot title. Default 1.5.
#' @param cex_signal_label numeric. Font size for ATAC/RNA type indicator
#'   text. Default 1.2.
#' @param export character or NULL. File path to auto-export the plot (e.g.
#'   \code{"plot.pdf"}, \code{"plot.eps"}, \code{"plot.png"}). When NULL
#'   (default), the plot is drawn on the current device only.
#' @param width numeric. Export width in inches. Default 10.
#' @param height numeric. Export height in inches. Default 8.
#' @return Invisible \code{NULL}. The function produces a base R multi-panel
#'   plot on the current graphics device.
#' @export
#' @examples
#' db <- epiRomicsS4(annotations = GenomicRanges::GRanges(),
#'   meta = data.frame(name = character(), type = character(),
#'     file = character(), stringsAsFactors = FALSE),
#'   genome = "hg38")
#' tryCatch(epiRomics_track_layer(db, 1, db,
#'   data.frame(file = character(), name = character(),
#'     color = character(), type = character())),
#'   error = function(e) message(e$message))
#' \donttest{
#' epiRomics_track_layer(enhanceosome, 1, epiRomics_dB,
#'   epiRomics_track_connection)
#' # With font size and export customization:
#' epiRomics_track_layer(enhanceosome, 1, epiRomics_dB,
#'   epiRomics_track_connection, cex_title = 2.0,
#'   mirror = FALSE, export = "figure.pdf")
#' }
epiRomics_track_layer <- function(epiRomics_putative_enhanceosome,
                                        epiRomics_index,
                                        epiRomics_dB,
                                        epiRomics_track_connection,
                                        epiRomics_keep_epitracks = TRUE,
                                        chromatin_states = NULL,
                                        gene_lookup = FALSE,
                                        show_bigwig = TRUE,
                                        show_chromatin = TRUE,
                                        show_annotations = TRUE,
                                        show_gene_model = TRUE,
                                        show_enhancer_highlight = TRUE,
                                        mirror = TRUE,
                                        signal_style = c("line", "polygon"),
                                        signal_layout = "auto",
                                        cex_cell_label = 1.4,
                                        cex_axis = 1.2,
                                        cex_coord = 1.3,
                                        cex_annotation = 1.1,
                                        cex_gene = 1.2,
                                        cex_title = 1.5,
                                        cex_signal_label = 1.2,
                                        export = NULL,
                                        width = 10,
                                        height = 8) {
  validate_epiRomics_dB(epiRomics_putative_enhanceosome, "epiRomics_track_layer")
  validate_numeric_param(epiRomics_index, "epiRomics_index",
    "epiRomics_track_layer", min_val = 1)
  validate_epiRomics_dB(epiRomics_dB, "epiRomics_track_layer")
  validate_logical_param(epiRomics_keep_epitracks, "epiRomics_keep_epitracks",
    "epiRomics_track_layer")

  if (!base::is.data.frame(epiRomics_track_connection)) {
    base::stop("epiRomics_track_connection must be a data frame")
  }
  if (epiRomics_index > base::length(epiRomics_putative_enhanceosome@annotations)) {
    base::stop(base::sprintf("epiRomics_index (%d) exceeds number of annotations (%d)",
      epiRomics_index,
      base::length(epiRomics_putative_enhanceosome@annotations)))
  }
  # Validate new visibility/layout params
  signal_style <- base::match.arg(signal_style,
    base::c("line", "polygon"))
  signal_layout <- base::match.arg(signal_layout,
    base::c("auto", "stacked", "mirrored"))
  if (!base::is.null(export) && !base::is.character(export)) {
    base::stop("export must be NULL or a file path string")
  }

  # Open export device if requested
  if (!base::is.null(export)) {
    ext <- base::tolower(tools::file_ext(export))
    if (ext == "pdf") {
      grDevices::pdf(export, width = width, height = height)
    } else if (ext == "eps" || ext == "ps") {
      grDevices::postscript(export, width = width, height = height,
        horizontal = FALSE, paper = "special")
    } else if (ext == "png") {
      grDevices::png(export, width = width, height = height,
        units = "in", res = 300)
    } else if (ext == "tiff" || ext == "tif") {
      grDevices::tiff(export, width = width, height = height,
        units = "in", res = 300)
    } else {
      base::stop("Unsupported export format: ", ext,
        ". Use pdf, eps, png, or tiff.")
    }
    base::on.exit(grDevices::dev.off(), add = TRUE)
  }

  # --- Extract focal region ---
  enh_save <- epiRomics_putative_enhanceosome
  enh <- epiRomics_putative_enhanceosome@annotations[epiRomics_index]
  enh_chr <- base::as.character(GenomeInfoDb::seqnames(enh))
  enh_s <- BiocGenerics::start(enh)
  enh_e <- BiocGenerics::end(enh)
  gene_symbol <- enh$SYMBOL
  genome <- epiRomics_dB@genome

  # Compute viewing window:
  # Gene lookup mode: gene body +/- 1000bp (tight focus on the gene itself).
  # Enhancer mode: include both enhancer and gene with 5kb padding.
  if (gene_lookup) {
    pad <- 1000L
    min_track <- base::floor((enh$geneStart - pad) / 100) * 100
    max_track <- base::ceiling((enh$geneEnd + pad) / 100) * 100
  } else {
    min_track <- base::floor((base::min(enh_s, enh$geneStart) - 5000L) / 100) * 100
    max_track <- base::ceiling((base::max(enh_e, enh$geneEnd) + 5000L) / 100) * 100
  }
  region_gr <- GenomicRanges::GRanges(
    seqnames = enh_chr,
    ranges = IRanges::IRanges(start = min_track, end = max_track)
  )
  xlims <- base::c(min_track, max_track)
  # Get the actual enhanceosome index ID (the name) for the label

  enh_id <- base::names(epiRomics_putative_enhanceosome@annotations)[epiRomics_index]

  # Colors
  hl_violet <- "#7B2D8E"
  dull_blue <- "#7F8FA6"
  enh_fill <- "#F0E6F680"   # translucent violet
  enh_border <- "#D8BFD8"

  # Helper: draw enhancer highlight that extends into margins for contiguous
  # column across panels. Uses xpd=TRUE to draw beyond the plot region into
  # the figure margins, so adjacent panels' highlights merge seamlessly.
  # In gene_lookup mode, the highlight is skipped entirely since we are

  # just viewing a gene locus, not highlighting a specific enhancer.
  .draw_highlight <- function() {
    if (gene_lookup) return(base::invisible(NULL))
    usr <- graphics::par("usr")
    y_ext <- (usr[4] - usr[3]) * 20
    graphics::rect(enh_s, usr[3] - y_ext, enh_e, usr[4] + y_ext,
      col = enh_fill, border = NA, xpd = TRUE)
    graphics::segments(enh_s, usr[3] - y_ext, enh_s, usr[4] + y_ext,
      col = enh_border, lwd = 0.5, xpd = TRUE)
    graphics::segments(enh_e, usr[3] - y_ext, enh_e, usr[4] + y_ext,
      col = enh_border, lwd = 0.5, xpd = TRUE)
  }

  # ====== DATA PREPARATION ======

  # 1. Collapsed gene model: one row per gene (union of all transcript exons)
  txdb_obj <- resolve_txdb(epiRomics_dB@txdb)
  orgdb <- base::tryCatch(base::get(epiRomics_dB@organism), error = function(e) NULL)

  # Get all genes overlapping the region
  all_txdb_genes <- GenomicFeatures::genes(txdb_obj)
  region_genes <- IRanges::subsetByOverlaps(all_txdb_genes, region_gr)

  gene_models <- base::list()
  if (base::length(region_genes) > 0 && !base::is.null(orgdb)) {
    gene_ids <- base::names(region_genes)

    # Map ENTREZID -> SYMBOL
    gene_syms <- base::tryCatch(
      AnnotationDbi::mapIds(orgdb, keys = gene_ids,
        keytype = "ENTREZID", column = "SYMBOL", multiVals = "first"),
      error = function(e) stats::setNames(base::rep(NA_character_,
        base::length(gene_ids)), gene_ids))

    # Secondary focal gene lookup: ensure focal gene is included
    if (!base::is.na(gene_symbol) && base::nzchar(gene_symbol) &&
        !base::any(!base::is.na(gene_syms) & gene_syms == gene_symbol)) {
      focal_eid <- base::tryCatch(
        AnnotationDbi::mapIds(orgdb, keys = gene_symbol,
          keytype = "SYMBOL", column = "ENTREZID",
          multiVals = "first"),
        error = function(e) NA_character_)
      if (!base::is.na(focal_eid) && focal_eid %in% gene_ids) {
        gene_syms[focal_eid] <- gene_symbol
      }
    }

    # Get exons grouped by gene for union model
    exons_by_gene <- GenomicFeatures::exonsBy(txdb_obj, by = "gene")

    gene_models <- base::vector("list", base::length(gene_ids))
    for (i in base::seq_along(gene_ids)) {
      gid <- gene_ids[i]
      sym <- gene_syms[i]
      if (base::is.na(sym)) sym <- gid

      # Merge overlapping exons into union gene model
      if (gid %in% base::names(exons_by_gene)) {
        gene_exons <- exons_by_gene[[gid]]
        merged <- GenomicRanges::reduce(gene_exons)
        gene_models[[i]] <- base::list(
          symbol = sym,
          gene_id = gid,
          start = BiocGenerics::start(region_genes[i]),
          end = BiocGenerics::end(region_genes[i]),
          strand = base::as.character(BiocGenerics::strand(region_genes[i])),
          ex_starts = BiocGenerics::start(merged),
          ex_ends = BiocGenerics::end(merged)
        )
      } else {
        # Gene record exists but no exons -- draw gene body only
        gene_models[[i]] <- base::list(
          symbol = sym,
          gene_id = gid,
          start = BiocGenerics::start(region_genes[i]),
          end = BiocGenerics::end(region_genes[i]),
          strand = base::as.character(BiocGenerics::strand(region_genes[i])),
          ex_starts = base::integer(0),
          ex_ends = base::integer(0)
        )
      }
    }
  }

  # Sort: focal gene last (drawn on top), others by gene width
  if (base::length(gene_models) > 0) {
    is_focal <- base::vapply(gene_models, function(g) {
      !base::is.na(g$symbol) && g$symbol == gene_symbol
    }, logical(1))
    widths <- base::vapply(gene_models, function(g) g$end - g$start, numeric(1))
    sort_order <- base::order(is_focal, widths)
    gene_models <- gene_models[sort_order]

    # Limit to 8 genes max (prioritizing focal gene)
    if (base::length(gene_models) > 8) {
      is_focal_sorted <- base::vapply(gene_models, function(g) {
        !base::is.na(g$symbol) && g$symbol == gene_symbol
      }, logical(1))
      focal_idx <- base::which(is_focal_sorted)
      other_idx <- base::setdiff(base::seq_along(gene_models), focal_idx)
      keep_idx <- base::sort(base::c(
        utils::head(other_idx, 8 - base::length(focal_idx)), focal_idx))
      gene_models <- gene_models[keep_idx]
    }

    # In gene_lookup mode, only show the focal gene (not readthrough
    # transcripts or overlapping genes like INS-IGF2)
    if (gene_lookup) {
      focal_only <- base::Filter(function(g) {
        !base::is.na(g$symbol) && g$symbol == gene_symbol
      }, gene_models)
      if (base::length(focal_only) > 0) gene_models <- focal_only
    }

    # In enhancer mode, expand window if gene models extend beyond
    if (!gene_lookup) {
      all_gm_min <- base::min(base::vapply(gene_models,
        function(g) g$start, numeric(1)))
      all_gm_max <- base::max(base::vapply(gene_models,
        function(g) g$end, numeric(1)))
      if (all_gm_min < min_track + 5000L) {
        min_track <- base::floor((all_gm_min - 5000L) / 100) * 100
      }
      if (all_gm_max > max_track - 5000L) {
        max_track <- base::ceiling((all_gm_max + 5000L) / 100) * 100
      }
      xlims <- base::c(min_track, max_track)
      region_gr <- GenomicRanges::GRanges(
        seqnames = enh_chr,
        ranges = IRanges::IRanges(start = min_track, end = max_track))
    }
  }
  n_genes <- base::max(1L, base::length(gene_models))

  # 2. BigWig signal data (single-pass import for each type)
  tc_atac <- epiRomics_track_connection[epiRomics_track_connection$type == "atac", , drop = FALSE]
  tc_rna <- epiRomics_track_connection[epiRomics_track_connection$type == "rna", , drop = FALSE]
  tc_hist <- epiRomics_track_connection[epiRomics_track_connection$type == "histone", , drop = FALSE]

  atac_data <- .import_bw_group_data(tc_atac, region_gr, xlims)
  rna_data <- .import_bw_group_data(tc_rna, region_gr, xlims)
  hist_data <- .import_bw_group_data(tc_hist, region_gr, xlims)

  # 3. Try to pair matched ATAC/RNA samples into mirrored panels
  #    Publication-style: positive Y = ATAC, negative Y = RNA
  mirrored_pairs <- base::list()
  unpaired_atac <- base::seq_len(base::nrow(tc_atac))
  unpaired_rna <- base::seq_len(base::nrow(tc_rna))
  if (!base::is.null(atac_data) && !base::is.null(rna_data)) {
    atac_prefix <- base::sub("[_.]?(ATAC|atac)$", "", atac_data$names)
    rna_prefix <- base::sub("[_.]?(RNA|rna)$", "", rna_data$names)
    used_atac <- base::logical(base::length(atac_prefix))
    used_rna <- base::logical(base::length(rna_prefix))
    for (ai in base::seq_along(atac_prefix)) {
      ri <- base::match(atac_prefix[ai], rna_prefix)
      if (!base::is.na(ri) && !used_rna[ri]) {
        mirrored_pairs <- base::c(mirrored_pairs, base::list(base::list(
          atac_idx = ai, rna_idx = ri,
          label = atac_prefix[ai])))
        used_atac[ai] <- TRUE
        used_rna[ri] <- TRUE
      }
    }
    unpaired_atac <- base::which(!used_atac)
    unpaired_rna <- base::which(!used_rna)
  }
  n_mirrored <- base::length(mirrored_pairs)

  # Apply mirror / signal_layout overrides
  if (signal_layout == "stacked" || !mirror) {
    # Force all signals into separate rows — no mirroring
    if (n_mirrored > 0) {
      for (mp in mirrored_pairs) {
        unpaired_atac <- base::c(unpaired_atac, mp$atac_idx)
        unpaired_rna <- base::c(unpaired_rna, mp$rna_idx)
      }
      mirrored_pairs <- base::list()
      n_mirrored <- 0L
    }
  } else if (signal_layout == "mirrored") {
    # Force mirroring of first ATAC + first RNA even if names don't match
    if (n_mirrored == 0 && !base::is.null(atac_data) &&
        !base::is.null(rna_data) && base::nrow(tc_atac) > 0 &&
        base::nrow(tc_rna) > 0) {
      mirrored_pairs <- base::list(base::list(
        atac_idx = 1L, rna_idx = 1L,
        label = atac_data$names[1]))
      unpaired_atac <- unpaired_atac[unpaired_atac != 1L]
      unpaired_rna <- unpaired_rna[unpaired_rna != 1L]
      n_mirrored <- 1L
    }
  }
  # signal_layout == "auto" keeps the default pairing behavior above

  # Apply show_bigwig toggle
  if (!show_bigwig) {
    mirrored_pairs <- base::list()
    n_mirrored <- 0L
    unpaired_atac <- base::integer(0)
    unpaired_rna <- base::integer(0)
    hist_data <- NULL
  }

  # Count total panels
  n_hist <- if (!base::is.null(hist_data)) base::length(hist_data$data) else 0L
  n_signal <- n_mirrored + base::length(unpaired_atac) +
    base::length(unpaired_rna) + n_hist
  n_annot <- 0
  annot_panels <- base::list()  # will collect {gr, label, color} entries

  if (epiRomics_keep_epitracks) {
    # Pre-filter annotations to viewing region
    region_annots <- IRanges::subsetByOverlaps(epiRomics_dB@annotations, region_gr)
    region_types <- region_annots$type
    db_prefix <- base::paste0(genome, "_custom_")
    meta <- enh_save@meta
    all_types <- base::unique(epiRomics_dB@annotations$type)

    # ---- Track assembly (ordered top-to-bottom) ----
    # Order: Enhancer index (if not gene_lookup) → Chromatin states
    #   (Active/Poised/Repressed/Unmarked) → FANTOM5 → UCNEs →
    #   Regulome Active → Regulome Super → T2D GWAS → TFs (ALL) → ncRNAs

    # 1. Enhancer index -- SKIP in gene_lookup mode or if toggled off
    if (!gene_lookup && show_enhancer_highlight) {
      enh_label <- if (!base::is.null(enh_id) && !base::is.na(enh_id)) {
        base::paste0("Enhancer #", enh_id)
      } else {
        base::paste0("Enhancer #", epiRomics_index)
      }
      annot_panels <- base::c(annot_panels, base::list(base::list(
        gr = enh, label = enh_label, color = hl_violet, force_color = TRUE)))
      n_annot <- n_annot + 1
    }

    # 2. Chromatin states -- show ALL specific states (not just broad categories)
    # Each state gets its own color and annotation track.
    cs_palette <- .chromatin_state_palette()
    # Display labels: 6 simplified chromatin states
    cs_labels <- base::c(
      active    = "Active",
      bivalent  = "Bivalent",
      poised    = "Poised",
      primed    = "Primed",
      repressed = "Repressed",
      unmarked  = "Unmarked"
    )
    # Fixed display order
    cs_display_order <- base::c(
      "active", "bivalent", "poised", "primed", "repressed", "unmarked")
    if (show_chromatin && !base::is.null(chromatin_states) &&
        base::is.data.frame(chromatin_states) &&
        base::nrow(chromatin_states) > 0) {
      region_cs <- chromatin_states[
        chromatin_states$seqnames == enh_chr &
        chromatin_states$end >= min_track &
        chromatin_states$start <= max_track, , drop = FALSE
      ]
      # ALWAYS show ALL states in fixed order (even if no data in this region)
      for (st in cs_display_order) {
        st_rows <- region_cs[region_cs$chromatin_state == st, , drop = FALSE]
        st_gr <- if (base::nrow(st_rows) > 0) {
          GenomicRanges::GRanges(
            seqnames = st_rows$seqnames,
            ranges = IRanges::IRanges(
              start = st_rows$start, end = st_rows$end))
        } else {
          GenomicRanges::GRanges()  # Empty GRanges for states with no data
        }
        st_label <- if (st %in% base::names(cs_labels)) cs_labels[st] else st
        st_color <- if (st %in% base::names(cs_palette)) cs_palette[st] else "#AAAAAA"
        annot_panels <- base::c(annot_panels, base::list(base::list(
          gr = st_gr, label = st_label,
          color = st_color,
          force_color = TRUE)))
        n_annot <- n_annot + 1
      }
    }

    # 2b. Context-specific overlays (promoter/gene body annotations)
    # These provide positional specificity on top of the 6-state classification.
    # CRITICAL: chromatin state regions can be wider than actual gene features.
    # We INTERSECT with actual TxDb gene boundaries to clip overlays:
    #   - "Active Gene Body"  → active regions ∩ gene bodies (excluding TSS ±2kb)
    if (show_chromatin && !base::is.null(chromatin_states) &&
        base::is.data.frame(chromatin_states) &&
        base::nrow(chromatin_states) > 0 &&
        "genomic_context" %in% base::names(chromatin_states) &&
        base::length(region_genes) > 0) {

      # Build TSS ±2kb windows from actual gene models in the view
      ctx_tss <- GenomicRanges::trim(base::suppressWarnings(
        GenomicRanges::resize(region_genes, width = 1L, fix = "start")))
      ctx_tss_win <- GenomicRanges::trim(base::suppressWarnings(
        GenomicRanges::resize(ctx_tss, width = 4000L, fix = "center")))
      ctx_tss_win <- GenomicRanges::reduce(ctx_tss_win)

      # Gene body = gene model MINUS TSS window
      ctx_gene_body <- GenomicRanges::setdiff(
        GenomicRanges::reduce(region_genes), ctx_tss_win)

      ctx_defs <- base::list(
        base::list(state = "active", clip_gr = ctx_gene_body,
          label = "Active Gene Body", color = "#66A61E")
      )
      for (ctx in ctx_defs) {
        ctx_rows <- region_cs[
          region_cs$chromatin_state == ctx$state, , drop = FALSE]
        if (base::nrow(ctx_rows) > 0) {
          ctx_gr <- GenomicRanges::GRanges(
            seqnames = ctx_rows$seqnames,
            ranges = IRanges::IRanges(
              start = ctx_rows$start, end = ctx_rows$end))
          # Intersect chromatin state regions with clipping boundary
          ctx_gr <- GenomicRanges::intersect(ctx_gr, ctx$clip_gr)
        } else {
          ctx_gr <- GenomicRanges::GRanges()
        }
        annot_panels <- base::c(annot_panels, base::list(base::list(
          gr = ctx_gr, label = ctx$label,
          color = ctx$color, force_color = TRUE)))
        n_annot <- n_annot + 1
      }
    }

    # 3-10: BED annotation tracks (gated by show_annotations)
    if (show_annotations) {

    # 3. FANTOM5 enhancers (reference database)
    fantom_key <- base::paste0(db_prefix, "fantom")
    if (fantom_key %in% all_types) {
      fantom_gr <- region_annots[region_types == fantom_key]
      annot_panels <- base::c(annot_panels, base::list(base::list(
        gr = fantom_gr, label = "FANTOM5", color = "#F39C12",
        force_color = TRUE)))
      n_annot <- n_annot + 1
    }

    # 4. UCNEs (ultra-conserved non-coding elements)
    ucne_key <- base::paste0(db_prefix, "ucnes")
    if (ucne_key %in% all_types) {
      ucne_gr <- region_annots[region_types == ucne_key]
      annot_panels <- base::c(annot_panels, base::list(base::list(
        gr = ucne_gr, label = "UCNEs", color = "#8E44AD",
        force_color = TRUE)))
      n_annot <- n_annot + 1
    }

    # 5. Regulome Active enhancers (Human Islet Active Enhancers)
    regulome_active_key <- base::paste0(db_prefix, "regulome_active")
    if (regulome_active_key %in% all_types) {
      ra_gr <- region_annots[region_types == regulome_active_key]
      annot_panels <- base::c(annot_panels, base::list(base::list(
        gr = ra_gr, label = "Regulome Active",
        color = "#1ABC9C", force_color = TRUE)))
      n_annot <- n_annot + 1
    }

    # 6. Regulome Super-enhancers (Human Islet Super-Enhancers)
    regulome_super_key <- base::paste0(db_prefix, "regulome_super")
    if (regulome_super_key %in% all_types) {
      rs_gr <- region_annots[region_types == regulome_super_key]
      annot_panels <- base::c(annot_panels, base::list(base::list(
        gr = rs_gr, label = "Regulome Super",
        color = "#008080", force_color = TRUE)))
      n_annot <- n_annot + 1
    }

    # 7. T2D GWAS and other GWAS tracks (diamond markers)
    gwas_key <- base::paste0(db_prefix, "t2d_gwas")
    if (gwas_key %in% all_types) {
      gwas_gr <- region_annots[region_types == gwas_key]
      annot_panels <- base::c(annot_panels, base::list(base::list(
        gr = gwas_gr, label = "T2D GWAS", color = "#E67E22",
        force_color = TRUE, point_style = TRUE)))
      n_annot <- n_annot + 1
    }
    gwas_names <- meta[meta$type == "gwas", "name"]
    gwas_names <- gwas_names[gwas_names != "t2d_gwas"]
    for (gwas_nm in gwas_names) {
      gk <- base::paste0(db_prefix, gwas_nm)
      if (gk %in% all_types) {
        gwas_gr <- region_annots[region_types == gk]
        annot_panels <- base::c(annot_panels, base::list(base::list(
          gr = gwas_gr,
          label = base::toupper(base::gsub("_", " ", gwas_nm)),
          color = "#E67E22", force_color = TRUE, point_style = TRUE)))
        n_annot <- n_annot + 1
      }
    }

    # H2A.Z is a histone (NOT a TF). It is NOT shown as a separate track.
    # Instead, it is used in chromatin state classification to help identify
    # regulatory elements that lack canonical histone marks.

    # 8. TF/ChIP tracks -- ALL TFs shown (empty bar if no peaks in window)
    tf_names <- meta[meta$type == "chip", "name"]
    tf_names <- tf_names[!base::grepl("h2a\\.?z", base::tolower(tf_names))]
    for (nm in tf_names) {
      pk <- region_annots[region_types == base::paste0(db_prefix, nm)]
      # Always add TF track even if empty (for visual consistency)
      if (base::length(pk) == 0) pk <- GenomicRanges::GRanges()
      annot_panels <- base::c(annot_panels, base::list(base::list(
        gr = pk, label = base::toupper(nm), color = dull_blue)))
      n_annot <- n_annot + 1
    }

    # 10. ncRNA tracks
    for (nm in meta[meta$type == "ncrna", "name"]) {
      nc <- region_annots[region_types == base::paste0(db_prefix, nm)]
      if (base::length(nc) > 0) {
        annot_panels <- base::c(annot_panels, base::list(base::list(
          gr = nc, label = base::toupper(base::gsub("_", " ", nm)),
          color = dull_blue)))
        n_annot <- n_annot + 1
      }
    }

    } # end if (show_annotations)
  }

  # Compute panel count (conditionally include gene model)
  n_gene_panels <- if (show_gene_model) 1L else 0L
  # gene + coord_bar + signals + annotations
  n_panels <- n_gene_panels + 1 + n_signal + n_annot

  # ====== DRAWING ======

  # Compute longest annotation label for dynamic left margin
  annot_labels <- base::vapply(annot_panels, function(ap) ap$label, "")
  max_label_nchar <- if (base::length(annot_labels) > 0) {
    base::max(base::nchar(annot_labels))
  } else {
    8L
  }
  # Margin width: ~0.55 lines per character for horizontal text at cex 0.65
  left_mar <- base::max(8, base::ceiling(max_label_nchar * 0.55) + 2)

  # Panel heights: gene=dynamic, coord=0.5, mirrored=3.5, signal=2, annotation=0.7
  n_unpaired_signal <- base::length(unpaired_atac) + base::length(unpaired_rna) +
    n_hist
  gene_panel_height <- base::max(1.5, n_genes * 0.8)
  panel_heights <- base::c(
    if (show_gene_model) gene_panel_height,  # gene model (1 row per gene)
    0.8,                                  # coordinate bar
    base::rep(3.5, n_mirrored),           # mirrored ATAC/RNA panels
    base::rep(2, n_unpaired_signal),      # individual signal tracks
    base::rep(0.7, n_annot)              # annotation bars (compact)
  )

  old_par <- graphics::par(no.readonly = TRUE)
  base::on.exit(graphics::par(old_par), add = TRUE)

  graphics::layout(base::matrix(base::seq_len(n_panels), ncol = 1),
    heights = panel_heights)
  graphics::par(mar = base::c(0.3, left_mar, 0.3, 1), oma = base::c(1, 0, 2, 0),
    bg = "white")

  # Panel 1: Gene model (collapsed: one row per gene, union of all exons)
  if (show_gene_model) {
  graphics::plot(1, type = "n", xlim = xlims,
    ylim = base::c(0.5, n_genes + 0.5),
    xaxt = "n", yaxt = "n", xlab = "", ylab = gene_symbol, bty = "n")
  .draw_highlight()
  if (base::length(gene_models) > 0) {
    view_span <- max_track - min_track
    arrow_half <- view_span * 0.008
    min_intron_width <- view_span * 0.015

    for (i in base::seq_along(gene_models)) {
      gm <- gene_models[[i]]
      is_focal_gene <- !base::is.na(gm$symbol) && gm$symbol == gene_symbol

      line_col <- if (is_focal_gene) "#1A1A1A" else "#999999"
      exon_col <- if (is_focal_gene) "#2C3E50" else "#BBBBBB"
      arrow_lwd <- if (is_focal_gene) 3.0 else 2.0
      end_lwd <- if (is_focal_gene) 3.5 else 2.5

      # Gene body intron line
      graphics::segments(gm$start, i, gm$end, i,
        col = line_col, lwd = if (is_focal_gene) 1 else 0.6)

      arrow_dir <- if (gm$strand == "+") 1 else -1

      # Direction arrows between exons (intron arrows)
      if (base::length(gm$ex_starts) >= 2) {
        sorted_idx <- base::order(gm$ex_starts)
        s_ex_s <- gm$ex_starts[sorted_idx]
        s_ex_e <- gm$ex_ends[sorted_idx]
        for (j in base::seq_len(base::length(s_ex_e) - 1)) {
          intron_mid <- (s_ex_e[j] + s_ex_s[j + 1]) / 2
          intron_w <- s_ex_s[j + 1] - s_ex_e[j]
          if (intron_w > min_intron_width) {
            graphics::arrows(
              intron_mid - arrow_dir * arrow_half, i,
              intron_mid + arrow_dir * arrow_half, i,
              length = 0.10, lwd = arrow_lwd, col = line_col, code = 2)
          }
        }
      }
      # Terminal direction arrow at 3' end
      end_x <- if (gm$strand == "+") gm$end else gm$start
      graphics::arrows(
        end_x - arrow_dir * arrow_half * 1.5, i,
        end_x + arrow_dir * arrow_half * 1.5, i,
        length = 0.12, lwd = end_lwd, col = line_col, code = 2)

      # Merged exon blocks (union of all transcript exons)
      if (base::length(gm$ex_starts) > 0) {
        graphics::rect(gm$ex_starts, i - 0.35, gm$ex_ends, i + 0.35,
          col = exon_col, border = NA)
      }

      # Label non-focal genes with their symbol
      if (!is_focal_gene && !base::is.na(gm$symbol)) {
        graphics::text(gm$start, i + 0.4, gm$symbol,
          cex = cex_gene, col = "#999999", adj = base::c(0, 0))
      }
    }
  }
  } # end if (show_gene_model)

  # Panel 2: Coordinate bar (chr, start, end) right below transcripts
  graphics::plot(1, type = "n", xlim = xlims, ylim = base::c(0, 1),
    xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
  .draw_highlight()
  # Draw a line spanning the full extent (raised to upper area)
  graphics::segments(xlims[1], 0.75, xlims[2], 0.75, lwd = 1.5, col = "#333333")
  # Tick marks at endpoints
  graphics::segments(xlims[1], 0.60, xlims[1], 0.90, lwd = 1.5, col = "#333333")
  graphics::segments(xlims[2], 0.60, xlims[2], 0.90, lwd = 1.5, col = "#333333")
  # Coordinate labels BELOW the line bar (well padded, xpd to prevent clipping)
  graphics::text(xlims[1], 0.15,
    base::format(xlims[1], big.mark = ",", scientific = FALSE),
    adj = base::c(0, 0.5), cex = cex_coord, col = "#333333", xpd = TRUE)
  graphics::text(xlims[2], 0.15,
    base::format(xlims[2], big.mark = ",", scientific = FALSE),
    adj = base::c(1, 0.5), cex = cex_coord, col = "#333333", xpd = TRUE)
  # Chromosome and genome build in the center, BELOW the line bar
  graphics::text(base::mean(xlims), 0.15,
    base::paste0(enh_chr, " (", genome, ")"),
    adj = base::c(0.5, 0.5), cex = cex_coord, col = "#333333", font = 2, xpd = TRUE)

  # Signal panels
  .draw_signal_panel <- function(df, ymax, col, label) {
    graphics::plot(1, type = "n", xlim = xlims,
      ylim = base::c(0, ymax),
      xaxt = "n", ylab = label, xlab = "", bty = "n")
    .draw_highlight()
    if (base::nrow(df) > 0) {
      ord <- base::order(df$pos)
      px <- df$pos[ord]
      py <- df$score[ord]
      if (signal_style == "polygon") {
        graphics::polygon(
          base::c(px[1], px, px[base::length(px)]),
          base::c(0, py, 0),
          col = grDevices::adjustcolor(col, alpha.f = 0.85),
          border = col, lwd = 0.3)
      } else {
        graphics::segments(px, 0, px, py, col = col, lwd = 1.5)
      }
    }
  }

  # Mirrored ATAC/RNA panels (publication style)
  # Both halves fill their side equally --signal is scaled to mirror_max
  # so that e.g. ATAC peak at 120 and RNA peak at 60 each fill to the edge.
  # Axis labels show actual data values (not the rescaled values).
  # SHARED max: all panels use the same ATAC max and same RNA max across

  # all samples so that panels are directly comparable (e.g. if Beta RNA
  # max=20 and Alpha RNA max=10, both panels scale to 20).
  if (n_mirrored > 0) {
    # Use the already-capped shared max from .import_bw_group_data()
    # (includes 98th percentile clipping for extreme outlier peaks)
    shared_atac_max <- atac_data$max
    shared_rna_max <- rna_data$max

    for (mp in mirrored_pairs) {
      atac_df <- atac_data$data[[mp$atac_idx]]
      rna_df <- rna_data$data[[mp$rna_idx]]
      # Use shared max across all samples for cross-panel comparability
      atac_max <- shared_atac_max
      rna_max <- shared_rna_max
      # Equal-length axes: use the larger of the two so both sides match
      mirror_max <- base::max(atac_max, rna_max)

      # Scale factors so each assay fills its half of the axis
      atac_scale <- if (atac_max > 0) mirror_max / atac_max else 1
      rna_scale <- if (rna_max > 0) mirror_max / rna_max else 1

      graphics::plot(1, type = "n", xlim = xlims,
        ylim = base::c(-mirror_max, mirror_max),
        xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
      .draw_highlight()

      # ATAC on positive Y
      if (base::nrow(atac_df) > 0) {
        ord_a <- base::order(atac_df$pos)
        ax <- atac_df$pos[ord_a]
        ay <- atac_df$score[ord_a] * atac_scale
        acol <- atac_data$colors[mp$atac_idx]
        if (signal_style == "polygon") {
          graphics::polygon(
            base::c(ax[1], ax, ax[base::length(ax)]),
            base::c(0, ay, 0),
            col = grDevices::adjustcolor(acol, alpha.f = 0.85),
            border = acol, lwd = 0.3)
        } else {
          graphics::segments(ax, 0, ax, ay, col = acol, lwd = 1.5)
        }
      }
      # RNA on negative Y
      if (base::nrow(rna_df) > 0) {
        ord_r <- base::order(rna_df$pos)
        rx <- rna_df$pos[ord_r]
        ry <- -(rna_df$score[ord_r] * rna_scale)
        rcol <- rna_data$colors[mp$rna_idx]
        if (signal_style == "polygon") {
          graphics::polygon(
            base::c(rx[1], rx, rx[base::length(rx)]),
            base::c(0, ry, 0),
            col = grDevices::adjustcolor(rcol, alpha.f = 0.85),
            border = rcol, lwd = 0.3)
        } else {
          graphics::segments(rx, 0, rx, ry, col = rcol, lwd = 1.5)
        }
      }
      # Zero line
      graphics::abline(h = 0, col = "#AAAAAA", lwd = 0.5)

      # Custom y-axis: show actual data max values at the symmetric extent
      graphics::axis(2, at = base::c(-mirror_max, 0, mirror_max),
        labels = base::c(rna_max, 0, atac_max), las = 1, cex.axis = cex_axis)
      # Label in margin
      graphics::mtext(mp$label, side = 2, line = 3, cex = cex_cell_label, las = 1)
      # Inline ATAC/RNA indicators
      usr <- graphics::par("usr")
      graphics::text(usr[2], mirror_max * 0.8, "ATAC", adj = base::c(1.1, 0.5),
        cex = cex_signal_label, col = atac_data$colors[mp$atac_idx], font = 2)
      graphics::text(usr[2], -mirror_max * 0.8, "RNA", adj = base::c(1.1, 0.5),
        cex = cex_signal_label, col = rna_data$colors[mp$rna_idx], font = 2)
    }
  }

  # Unpaired ATAC panels
  if (!base::is.null(atac_data) && base::length(unpaired_atac) > 0) {
    for (j in unpaired_atac) {
      .draw_signal_panel(atac_data$data[[j]], atac_data$max,
        atac_data$colors[j], atac_data$names[j])
    }
  }
  # Unpaired RNA panels
  if (!base::is.null(rna_data) && base::length(unpaired_rna) > 0) {
    for (j in unpaired_rna) {
      .draw_signal_panel(rna_data$data[[j]], rna_data$max,
        rna_data$colors[j], rna_data$names[j])
    }
  }
  # Histone panels
  if (!base::is.null(hist_data)) {
    for (j in base::seq_along(hist_data$data)) {
      .draw_signal_panel(hist_data$data[[j]], hist_data$max,
        hist_data$colors[j], hist_data$names[j])
    }
  }

  # Annotation panels --use horizontal mtext labels to avoid truncation
  if (n_annot > 0) {
    for (k in base::seq_along(annot_panels)) {
      ap <- annot_panels[[k]]
      is_last <- (k == base::length(annot_panels) && n_signal == 0) ||
                 (k == base::length(annot_panels))

      graphics::plot(1, type = "n", xlim = xlims, ylim = base::c(0, 1),
        xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
      # Horizontal label in left margin (las=1) --avoids vertical clipping
      graphics::mtext(ap$label, side = 2, line = 0.5, cex = cex_annotation, las = 1)
      .draw_highlight()

      gr <- ap$gr
      if (base::length(gr) > 0) {
        if (base::isTRUE(ap$point_style)) {
          # GWAS variants --diamond markers for single-nucleotide annotations
          gr_mid <- (BiocGenerics::start(gr) + BiocGenerics::end(gr)) / 2
          graphics::points(gr_mid, base::rep(0.5, base::length(gr_mid)),
            pch = 18, col = ap$color, cex = 1.5)
        } else if (base::isTRUE(ap$force_color)) {
          # Chromatin state / enhancer index --use assigned color directly
          graphics::rect(BiocGenerics::start(gr), 0, BiocGenerics::end(gr), 1,
            col = ap$color, border = NA)
        } else {
          # TF/ncRNA: draw ALL peaks in base color first (dull blue)
          gr_starts <- BiocGenerics::start(gr)
          gr_ends <- BiocGenerics::end(gr)
          graphics::rect(gr_starts, 0, gr_ends, 1,
            col = ap$color, border = NA)

          # Then overlay violet ONLY on the intersection with the enhancer
          # region --don't recolor the entire peak, just the overlapping slice.
          # Skip in gene_lookup mode: the synthetic 1-bp "enhancer" at gene
          # midpoint creates a thin vertical line through peaks.
          if (!gene_lookup) {
            clip_s <- base::pmax(gr_starts, enh_s)
            clip_e <- base::pmin(gr_ends, enh_e)
            has_overlap <- clip_s < clip_e
            if (base::any(has_overlap)) {
              graphics::rect(clip_s[has_overlap], 0, clip_e[has_overlap], 1,
                col = hl_violet, border = NA)
            }
          }
        }
      }
    }
  }

  # Title (coordinates are now in the dedicated bar below transcripts)
  title_text <- if (gene_lookup) {
    gene_symbol
  } else {
    base::paste0(gene_symbol, " enhanceosome region")
  }
  graphics::mtext(title_text,
    side = 3, outer = TRUE, line = 0.3, cex = cex_title, font = 2)

  base::invisible(NULL)
}

#' @rdname epiRomics_track_layer
#' @export
epiRomics_track_layer_fast <- epiRomics_track_layer


#' Visualize any gene locus with multi-track BigWig overlay
#'
#' A convenience wrapper around \code{\link{epiRomics_track_layer}} that
#' looks up a gene by symbol and creates a multi-track view without
#' requiring the full enhanceosome pipeline.
#'
#' The function queries the TxDb for the official gene body coordinates,
#' builds a synthetic single-region epiRomics object, and passes it to
#' \code{epiRomics_track_layer} for rendering. The viewing window spans
#' the gene body plus \code{padding} on each side.
#'
#' @param gene_symbol Character. HGNC gene symbol (e.g. \code{"INS"},
#'   \code{"GCG"}, \code{"PDX1"}).
#' @param epiRomics_dB An epiRomics S4 database object.
#' @param epiRomics_track_connection A data.frame from the BigWig CSV
#'   sheet (columns: path, name, color, type).
#' @param chromatin_states Optional. Output of
#'   \code{\link{epiRomics_chromatin_states}}.
#' @param padding Integer. Base pairs of padding around the gene body
#'   (default 1000).
#' @param show_bigwig Logical. Show BigWig signal tracks (default TRUE).
#' @param show_chromatin Logical. Show chromatin state tracks (default TRUE).
#' @param show_annotations Logical. Show BED annotation tracks (default TRUE).
#' @param show_gene_model Logical. Show gene model panel (default TRUE).
#' @param show_enhancer_highlight Logical. Show enhancer highlight (default TRUE).
#' @param mirror Logical. Enable mirrored ATAC/RNA layout (default TRUE).
#' @param signal_style Character. Signal rendering style: \code{"line"}
#'   (default) draws vertical bars at each position (IGV/UCSC browser style),
#'   \code{"polygon"} draws filled area charts.
#' @param signal_layout Character. Signal layout mode: \code{"auto"},
#'   \code{"stacked"}, or \code{"mirrored"}.
#' @param cex_cell_label Numeric. Font size for cell type labels (default 1.4).
#' @param cex_axis Numeric. Font size for axis labels (default 1.2).
#' @param cex_coord Numeric. Font size for coordinate bar (default 1.3).
#' @param cex_annotation Numeric. Font size for annotation labels (default 1.1).
#' @param cex_gene Numeric. Font size for gene model labels (default 1.2).
#' @param cex_title Numeric. Font size for plot title (default 1.5).
#' @param cex_signal_label Numeric. Font size for signal indicators (default 1.2).
#' @param export Character or NULL. File path for export (pdf, eps, png, tiff).
#' @param width Numeric. Export width in inches (default 10).
#' @param height Numeric. Export height in inches (default 8).
#' @return Invisible NULL. A multi-panel base-R plot is drawn on the
#'   current graphics device.
#' @export
#' @examples
#' db <- epiRomicsS4(annotations = GenomicRanges::GRanges(),
#'   meta = data.frame(name = character(), type = character(),
#'     file = character(), stringsAsFactors = FALSE),
#'   genome = "hg38")
#' tc <- data.frame(path = character(), name = character(),
#'   color = character(), type = character(), stringsAsFactors = FALSE)
#' tryCatch(epiRomics_track_layer_gene("INS", db, tc),
#'   error = function(e) message(e$message))
#' \donttest{
#' epiRomics_track_layer_gene("INS", epiRomics_dB, track_connection)
#' epiRomics_track_layer_gene("GCG", epiRomics_dB, track_connection,
#'   chromatin_states = cs, padding = 10000, export = "gcg_locus.pdf")
#' }
epiRomics_track_layer_gene <- function(gene_symbol,
                                        epiRomics_dB,
                                        epiRomics_track_connection,
                                        chromatin_states = NULL,
                                        padding = 1000L,
                                        show_bigwig = TRUE,
                                        show_chromatin = TRUE,
                                        show_annotations = TRUE,
                                        show_gene_model = TRUE,
                                        show_enhancer_highlight = TRUE,
                                        mirror = TRUE,
                                        signal_style = c("line", "polygon"),
                                        signal_layout = "auto",
                                        cex_cell_label = 1.4,
                                        cex_axis = 1.2,
                                        cex_coord = 1.3,
                                        cex_annotation = 1.1,
                                        cex_gene = 1.2,
                                        cex_title = 1.5,
                                        cex_signal_label = 1.2,
                                        export = NULL,
                                        width = 10,
                                        height = 8) {
  if (!base::is.character(gene_symbol) || base::length(gene_symbol) != 1) {
    base::stop("gene_symbol must be a single character string")
  }
  validate_epiRomics_dB(epiRomics_dB, "epiRomics_track_layer_gene")
  if (!base::is.data.frame(epiRomics_track_connection)) {
    base::stop("epiRomics_track_connection must be a data frame")
  }

  # Resolve TxDb
  txdb_obj <- resolve_txdb(epiRomics_dB@txdb)
  orgdb <- base::tryCatch(
    base::get(epiRomics_dB@organism),
    error = function(e) base::stop("Cannot resolve organism database: ",
      epiRomics_dB@organism)
  )

  # Look up gene -> ENTREZID
  gene_id <- base::tryCatch({
    res <- AnnotationDbi::select(orgdb, keys = gene_symbol,
      keytype = "SYMBOL", columns = "ENTREZID")
    res$ENTREZID[!base::is.na(res$ENTREZID)][1]
  }, error = function(e) NULL)

  if (base::is.null(gene_id) || base::is.na(gene_id)) {
    base::stop(base::sprintf("Gene symbol '%s' not found in organism database",
      gene_symbol))
  }

  # Use genes() for official gene body coordinates (not transcript extents
  # which can include readthrough transcripts like INS-IGF2)
  gene_gr <- GenomicFeatures::genes(txdb_obj,
    filter = base::list(gene_id = gene_id))
  if (base::length(gene_gr) == 0) {
    base::stop(base::sprintf("No gene record found for '%s'", gene_symbol))
  }

  gene_chr <- base::as.character(GenomeInfoDb::seqnames(gene_gr)[1])
  gene_start <- BiocGenerics::start(gene_gr)[1]
  gene_end <- BiocGenerics::end(gene_gr)[1]

  # Build a synthetic single-region epiRomics object
  # Place a 1-bp "enhancer" at the gene midpoint so the track_layer

  # centers the view on the gene
  mid <- base::as.integer((gene_start + gene_end) / 2)
  synth_gr <- GenomicRanges::GRanges(
    seqnames = gene_chr,
    ranges = IRanges::IRanges(start = mid, end = mid + 1L)
  )
  synth_gr$SYMBOL <- gene_symbol
  synth_gr$geneStart <- gene_start
  synth_gr$geneEnd <- gene_end
  synth_gr$geneChr <- base::sub("^chr", "", gene_chr)
  base::names(synth_gr) <- "gene_lookup"

  # Build a minimal epiRomics S4 wrapper
  synth_obj <- methods::new("epiRomicsS4")
  synth_obj@annotations <- synth_gr
  synth_obj@meta <- epiRomics_dB@meta
  synth_obj@txdb <- epiRomics_dB@txdb
  synth_obj@organism <- epiRomics_dB@organism
  synth_obj@genome <- epiRomics_dB@genome

  epiRomics_track_layer(
    synth_obj,
    epiRomics_index = 1L,
    epiRomics_dB = epiRomics_dB,
    epiRomics_track_connection = epiRomics_track_connection,
    chromatin_states = chromatin_states,
    gene_lookup = TRUE,
    show_bigwig = show_bigwig,
    show_chromatin = show_chromatin,
    show_annotations = show_annotations,
    show_gene_model = show_gene_model,
    show_enhancer_highlight = show_enhancer_highlight,
    mirror = mirror,
    signal_style = signal_style,
    signal_layout = signal_layout,
    cex_cell_label = cex_cell_label,
    cex_axis = cex_axis,
    cex_coord = cex_coord,
    cex_annotation = cex_annotation,
    cex_gene = cex_gene,
    cex_title = cex_title,
    cex_signal_label = cex_signal_label,
    export = export,
    width = width,
    height = height
  )
}
