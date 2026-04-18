#' Validate quick_view input parameters
#'
#' Checks gene/region exclusivity, bw_paths type and names,
#' labels length, file existence, and color vector length.
#' Returns updated bw_paths with label names applied if
#' provided.
#'
#' @param gene character or NULL
#' @param region list or NULL
#' @param bw_paths named character vector of BigWig file paths
#' @param labels character vector or NULL
#' @param colors character vector or NULL
#' @return character vector (bw_paths, possibly with updated
#'   names)
#' @noRd
.validate_quick_view_inputs <- function(gene, region,
                                        bw_paths, labels,
                                        colors) {
  if (base::is.null(gene) && base::is.null(region)) {
    base::stop(
      "Exactly one of 'gene' or 'region' must be provided"
    )
  }
  if (!base::is.null(gene) && !base::is.null(region)) {
    base::stop("Provide either 'gene' or 'region', not both")
  }
  if (!base::is.character(bw_paths) ||
      base::length(bw_paths) == 0) {
    base::stop(
      "bw_paths must be a non-empty named character vector of BigWig paths"
    )
  }
  if (base::is.null(base::names(bw_paths))) {
    base::stop(
      "bw_paths must be named (names are used as sample labels)"
    )
  }

  if (!base::is.null(labels)) {
    if (base::length(labels) != base::length(bw_paths)) {
      base::stop(
        "labels must have the same length as bw_paths"
      )
    }
    base::names(bw_paths) <- labels
  }

  missing <- bw_paths[!base::file.exists(bw_paths)]
  if (base::length(missing) > 0) {
    base::stop("BigWig file(s) not found: ",
      base::paste(missing, collapse = ", "))
  }

  base::return(bw_paths)
}


#' Build track connection data.frame for quick_view
#'
#' Infers signal types from sample name patterns (RNA vs ATAC)
#' and resolves colors from provided vector or a default
#' colorblind-friendly palette.
#'
#' @param bw_paths named character vector of BigWig file paths
#' @param colors character vector or NULL
#' @return data.frame with columns: path, name, color, type
#' @noRd
.build_quick_view_track_connection <- function(bw_paths,
                                               colors) {
  sample_names <- base::names(bw_paths)
  types <- base::ifelse(
    base::grepl(
      "rna|rnaseq|expression", sample_names,
      ignore.case = TRUE
    ),
    "rna", "atac"
  )

  n_bw <- base::length(bw_paths)
  if (!base::is.null(colors)) {
    if (base::length(colors) != n_bw) {
      base::stop(
        "colors must have the same length as bw_paths"
      )
    }
  } else {
    default_colors <- base::c(
      "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
      "#FF7F00", "#A65628", "#F781BF", "#999999"
    )
    colors <- default_colors[
      ((base::seq_len(n_bw) - 1L) %%
        base::length(default_colors)) + 1L
    ]
  }

  base::data.frame(
    path = base::unname(bw_paths),
    name = sample_names,
    color = colors,
    type = types,
    stringsAsFactors = FALSE
  )
}


#' Build a minimal epiRomicsS4 database for quick_view
#'
#' Creates a minimal S4 object with TxDb, organism, genome,
#' dummy meta, and a single-range annotation to satisfy
#' validate_database.
#'
#' @param genome character, genome assembly name
#' @param txdb character, TxDb package::object string
#' @param organism_pkg character, organism annotation package
#'   name
#' @return epiRomicsS4 object
#' @noRd
.build_minimal_dB <- function(genome, txdb, organism_pkg) {
  minimal_dB <- methods::new("epiRomicsS4")
  minimal_dB@txdb <- txdb
  minimal_dB@organism <- organism_pkg
  minimal_dB@genome <- genome
  minimal_dB@meta <- base::data.frame(
    name = "quick_view", type = "functional",
    file = "none", stringsAsFactors = FALSE
  )
  minimal_dB@annotations <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1L, end = 1L)
  )
  base::return(minimal_dB)
}


#' Build synthetic enhanceosome object for region mode
#'
#' Validates region format and constructs a synthetic
#' epiRomicsS4 with a single-range enhanceosome for
#' track_layer compatibility.
#'
#' @param region list with chr, start, end elements
#' @param genome character, genome assembly name
#' @param txdb character, TxDb package::object string
#' @param organism_pkg character, organism annotation package
#'   name
#' @param minimal_dB epiRomicsS4, minimal database for meta
#'   reference
#' @return epiRomicsS4 synthetic enhanceosome object
#' @noRd
.build_region_synth_obj <- function(region, genome, txdb,
                                     organism_pkg,
                                     minimal_dB) {
  required <- base::c("chr", "start", "end")
  if (!base::is.list(region) ||
      !base::all(required %in% base::names(region))) {
    base::stop(
      "region must be a list with 'chr', 'start', 'end' elements"
    )
  }
  chr <- base::as.character(region$chr)
  r_start <- base::as.integer(region$start)
  r_end <- base::as.integer(region$end)

  mid <- base::as.integer((r_start + r_end) / 2)
  synth_gr <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges = IRanges::IRanges(
      start = mid, end = mid + 1L
    )
  )
  synth_gr$SYMBOL <- base::paste0(
    chr, ":", r_start, "-", r_end
  )
  synth_gr$geneStart <- r_start
  synth_gr$geneEnd <- r_end
  synth_gr$geneChr <- base::sub("^chr", "", chr)
  base::names(synth_gr) <- "region_lookup"

  synth_obj <- methods::new("epiRomicsS4")
  synth_obj@annotations <- synth_gr
  synth_obj@meta <- minimal_dB@meta
  synth_obj@txdb <- txdb
  synth_obj@organism <- organism_pkg
  synth_obj@genome <- genome
  base::return(synth_obj)
}


#' Quick standalone gene/region visualization from BigWig
#' files
#'
#' A zero-configuration entry point for visualizing any
#' genomic locus with BigWig signal tracks. Requires only a
#' gene symbol (or genomic coordinates) and one or more BigWig
#' file paths. No database setup, no track connection CSV, no
#' chromatin states --- just signal overlaid on the gene model.
#'
#' This is the simplest way to use epiRomics for exploratory
#' visualization. Specify a gene by symbol (e.g.,
#' \code{"INS"}) or a region by coordinates. The function
#' auto-resolves the TxDb annotation database for the
#' specified genome assembly.
#'
#' @param gene Character or NULL. HGNC gene symbol (e.g.,
#'   \code{"INS"}, \code{"GCG"}). Exactly one of \code{gene}
#'   or \code{region} must be provided.
#' @param region List or NULL. Genomic coordinates as
#'   \code{list(chr = "chr11", start = 2159779,
#'   end = 2161209)}. Exactly one of \code{gene} or
#'   \code{region} must be provided.
#' @param bw_paths Named character vector of BigWig file
#'   paths. Names are used as sample labels. If names contain
#'   both \code{"atac"} and \code{"rna"} (case-insensitive),
#'   mirror pairing is attempted.
#' @param labels Character vector or NULL. Override sample
#'   labels. If NULL, uses \code{names(bw_paths)}.
#' @param colors Character vector or NULL. Override sample
#'   colors. If NULL, uses a default colorblind-friendly
#'   palette. Must have the same length as \code{bw_paths}.
#' @param genome Character. String name of the genome
#'   assembly associated with \code{bw_paths} (e.g.
#'   \code{"hg38"}, \code{"mm10"}, \code{"rn6"},
#'   \code{"dm6"}). The value must match the assembly
#'   recorded in the supplied \code{txdb}; a validator
#'   compares the string against
#'   \code{GenomeInfoDb::genome(txdb)}. Defaults to
#'   \code{"hg38"} for convenience with built-in fallbacks.
#' @param txdb Character or NULL. TxDb specification in
#'   \code{"package::object"} form (e.g.
#'   \code{"TxDb.Rnorvegicus.UCSC.rn6.refGene::
#'   TxDb.Rnorvegicus.UCSC.rn6.refGene"}).
#'   If NULL, auto-resolved for the built-in convenience
#'   genomes \code{"hg38"} and \code{"mm10"}; for any
#'   other genome, \code{txdb} must be supplied
#'   explicitly.
#' @param organism Character or NULL. org.db specification
#'   (e.g. \code{"org.Rn.eg.db"}). If NULL, auto-resolved
#'   for the built-in convenience genomes \code{"hg38"}
#'   and \code{"mm10"}; for any other genome,
#'   \code{organism} must be supplied explicitly.
#' @param padding Integer. Base pairs of padding around
#'   gene/region (default 5000).
#' @param mirror Logical. Enable mirrored ATAC/RNA layout
#'   (default FALSE).
#' @param signal_style Character. Signal rendering style:
#'   \code{"line"} (default) draws vertical bars at each
#'   position (IGV/UCSC browser style), \code{"polygon"}
#'   draws filled area charts.
#' @param signal_layout Character. Signal layout mode:
#'   \code{"auto"}, \code{"stacked"}, or \code{"mirrored"}.
#' @param cex_cell_label Numeric. Font size for cell type
#'   labels (default 1.4).
#' @param cex_axis Numeric. Font size for axis labels
#'   (default 1.2).
#' @param cex_coord Numeric. Font size for coordinate bar
#'   (default 1.3).
#' @param cex_annotation Numeric. Font size for annotation
#'   labels (default 1.1).
#' @param cex_gene Numeric. Font size for gene model labels
#'   (default 1.2).
#' @param cex_title Numeric. Font size for plot title
#'   (default 1.5).
#' @param cex_signal_label Numeric. Font size for signal
#'   indicators (default 1.2).
#' @param quantile_cap numeric. Percentile for capping
#'   extreme signal peaks (default 0.98). Peaks above this
#'   percentile are clipped to prevent axis compression.
#' @param scale_factor numeric. Y-axis headroom multiplier
#'   (default 1.1). Values above 1.0 add whitespace above
#'   the tallest signal peak.
#' @param export Character or NULL. File path for export
#'   (pdf, eps, png, tiff).
#' @param width Numeric. Export width in inches (default 10).
#' @param height Numeric. Export height in inches
#'   (default 8).
#' @return Invisible NULL. A multi-panel base-R plot is drawn
#'   on the current graphics device (or exported to file if
#'   \code{export} is set).
#' @export
#' @examples
#' ## plot_quick_view requires a TxDb (Suggests) and a real BigWig file.
#' ## This example confirms that input validation fires correctly.
#' if (requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
#'   bw_path <- make_example_bigwig()
#'   tryCatch(
#'     plot_quick_view(
#'       region = list(chr = "chr1", start = 1000L, end = 51000L),
#'       bw_paths = c(Sample = bw_path),
#'       genome   = "hg38"
#'     ),
#'     error = function(e) message(e$message)
#'   )
#'   file.remove(bw_path)
#' }
plot_quick_view <- function(gene = NULL,
                                  region = NULL,
                                  bw_paths,
                                  labels = NULL,
                                  colors = NULL,
                                  genome = "hg38",
                                  txdb = NULL,
                                  organism = NULL,
                                  padding = 5000L,
                                  mirror = FALSE,
                                  signal_style = c(
                                    "line", "polygon"
                                  ),
                                  signal_layout = "auto",
                                  cex_cell_label = 1.4,
                                  cex_axis = 1.2,
                                  cex_coord = 1.3,
                                  cex_annotation = 1.1,
                                  cex_gene = 1.2,
                                  cex_title = 1.5,
                                  cex_signal_label = 1.2,
                                  quantile_cap = 0.98,
                                  scale_factor = 1.1,
                                  export = NULL,
                                  width = 10,
                                  height = 8) {

  bw_paths <- .validate_quick_view_inputs(
    gene, region, bw_paths, labels, colors
  )

  # --- Validate genome and resolve TxDb / org.db ---
  # genome is organism-agnostic. For the built-in convenience genomes
  # ("hg38", "mm10") we auto-resolve TxDb and org.db if the user omitted
  # them. For any other genome the user must supply both explicitly.
  # When a txdb is supplied, a consistency check via
  # GenomeInfoDb::genome() ensures the declared genome string matches
  # the TxDb's recorded assembly.
  validate_character_param(
    genome, "genome", "plot_quick_view"
  )
  if (base::is.null(txdb)) {
    txdb <- base::switch(
      genome,
      "hg38" = base::paste0(
        "TxDb.Hsapiens.UCSC.hg38.knownGene::",
        "TxDb.Hsapiens.UCSC.hg38.knownGene"
      ),
      "mm10" = base::paste0(
        "TxDb.Mmusculus.UCSC.mm10.knownGene::",
        "TxDb.Mmusculus.UCSC.mm10.knownGene"
      ),
      .stop_no_txdb_default(genome)
    )
  }
  if (base::is.null(organism)) {
    organism_pkg <- base::switch(
      genome,
      "hg38" = "org.Hs.eg.db",
      "mm10" = "org.Mm.eg.db",
      .stop_no_orgdb_default(genome)
    )
  } else {
    validate_character_param(
      organism, "organism", "plot_quick_view"
    )
    organism_pkg <- organism
  }
  validate_genome_matches_txdb(
    genome, txdb, "plot_quick_view"
  )

  track_connection <- .build_quick_view_track_connection(
    bw_paths, colors
  )
  minimal_dB <- .build_minimal_dB(
    genome, txdb, organism_pkg
  )

  # Shared display params for delegation
  display_params <- base::list(
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
    quantile_cap = quantile_cap,
    scale_factor = scale_factor,
    export = export,
    width = width,
    height = height
  )

  if (!base::is.null(gene)) {
    # --- Gene mode: delegate to track_layer_gene ---
    base::do.call(
      plot_gene_tracks,
      base::c(base::list(
        gene_symbol = gene,
        database = minimal_dB,
        track_connection = track_connection,
        chromatin_states = NULL,
        padding = padding,
        show_bigwig = TRUE,
        show_chromatin = FALSE,
        show_annotations = FALSE,
        show_gene_model = TRUE,
        show_enhancer_highlight = FALSE
      ), display_params)
    )
  } else {
    # --- Region mode: build synthetic enhanceosome ---
    synth_obj <- .build_region_synth_obj(
      region, genome, txdb, organism_pkg, minimal_dB
    )
    base::do.call(
      plot_tracks,
      base::c(base::list(
        putative_enhanceosome = synth_obj,
        index = 1L,
        database = minimal_dB,
        track_connection = track_connection,
        chromatin_states = NULL,
        gene_lookup = TRUE,
        show_bigwig = TRUE,
        show_chromatin = FALSE,
        show_annotations = FALSE,
        show_gene_model = TRUE,
        show_enhancer_highlight = FALSE
      ), display_params)
    )
  }
}

## Internal: raise a clean error when genome has no default TxDb/org.db
## package wired in. Kept as a named helper so condition-signal calls
## pass a single-literal format to sprintf/stop (BiocCheck-friendly).
.stop_no_txdb_default <- function(genome) {
  fmt <- paste0(
    "plot_quick_view: no built-in TxDb default for genome '%s'. ",
    "Please supply 'txdb' explicitly as 'package::object' ",
    "(e.g. 'TxDb.Rnorvegicus.UCSC.rn6.refGene::",
    "TxDb.Rnorvegicus.UCSC.rn6.refGene')."
  )
  base::stop(base::sprintf(fmt, genome))
}

.stop_no_orgdb_default <- function(genome) {
  fmt <- paste0(
    "plot_quick_view: no built-in org.db default for genome '%s'. ",
    "Please supply 'organism' explicitly (e.g. 'org.Rn.eg.db')."
  )
  base::stop(base::sprintf(fmt, genome))
}
