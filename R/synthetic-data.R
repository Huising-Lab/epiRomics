# =============================================================================
# Exported example-data helpers for epiRomics man-page examples and vignettes
#
# These functions build fully self-contained synthetic epiRomicsS4 objects and
# BigWig files from in-memory GRanges — no network access, no external files.
# They are shipped with the package so that every @examples block in every
# exported function can run cleanly under R CMD check --as-cran --run-donttest.
#
# They also power the vignettes (via devtools::load_all() / library(epiRomics))
# and the testthat helpers in tests/testthat/helper-synthetic-builders.R.
#
# Author: Dr. Alex M. Mawla, PhD <ammawla@ucdavis.edu>
# =============================================================================

# ---------------------------------------------------------------------------
# Internal builder (mirrors make_synthetic_dB in helper-synthetic-builders.R)
# ---------------------------------------------------------------------------
.build_epiRomicsS4_from_lists <- function(marks_list,
                                           genome = "hg38",
                                           chip_list = NULL,
                                           functional_list = NULL) {

  build_entries <- function(entry_list, type_label) {
    if (base::is.null(entry_list) || base::length(entry_list) == 0L) {
      return(list(
        gr = GenomicRanges::GRanges(),
        meta = base::data.frame()
      ))
    }
    nms <- base::names(entry_list)
    typed <- base::lapply(nms, function(nm) {
      gr <- entry_list[[nm]]
      gr$type <- base::paste0(genome, "_custom_", nm)
      gr
    })
    combined_gr <- base::do.call(base::c, typed)
    meta_df <- base::data.frame(
      name = nms,
      type = base::rep(type_label, base::length(nms)),
      stringsAsFactors = FALSE
    )
    list(gr = combined_gr, meta = meta_df)
  }

  hist_e <- build_entries(marks_list,      "histone")
  chip_e <- build_entries(chip_list,       "chip")
  func_e <- build_entries(functional_list, "functional")

  all_gr <- base::do.call(base::c, base::list(hist_e$gr, chip_e$gr, func_e$gr))
  meta <- base::rbind(hist_e$meta, chip_e$meta, func_e$meta)

  methods::new(
    "epiRomicsS4",
    annotations = all_gr,
    meta        = meta,
    genome      = genome,
    txdb = paste0("TxDb.Hsapiens.UCSC.hg38.knownGene::",
                  "TxDb.Hsapiens.UCSC.hg38.knownGene"),
    organism    = "org.Hs.eg.db"
  )
}


# ---------------------------------------------------------------------------
# 1. make_example_database
# ---------------------------------------------------------------------------

#' Build a synthetic epiRomicsS4 database for use in examples
#'
#' Constructs a fully populated \code{epiRomicsS4} object from in-memory
#' \code{GRanges}, with no network access and no external files required.
#' The returned object contains five histone marks (h3k4me1, h3k27ac,
#' h3k27me3, h3k4me3, h3k36me3), two ChIP-seq TF tracks (TF1, TF2), and
#' one functional annotation track (fantom), all anchored on chr1 of hg38.
#'
#' This function is the canonical data source for all \code{@@examples} blocks
#' in the epiRomics package. It is also used in vignettes and testthat helpers.
#'
#' @param genome Character string naming the genome assembly (default
#'   \code{"hg38"}). Only the label is changed; the synthetic GRanges are
#'   always on chr1 regardless of this value, which is sufficient for
#'   example purposes.
#' @return An \code{epiRomicsS4} object with:
#'   \describe{
#'     \item{annotations}{GRanges containing all synthetic peak calls}
#'     \item{meta}{data.frame with columns \code{name} and \code{type}}
#'     \item{genome}{the value of \code{genome}}
#'     \item{txdb}{\code{"TxDb.Hsapiens.UCSC.hg38.knownGene::
#'       TxDb.Hsapiens.UCSC.hg38.knownGene"}}
#'     \item{organism}{\code{"org.Hs.eg.db"}}
#'   }
#' @family synthetic example data helpers
#' @seealso \code{\link{make_example_putative_enhancers}},
#'   \code{\link{make_example_enhanceosome}},
#'   \code{\link{make_example_bigwig}}
#' @export
#' @examples
#' db <- make_example_database()
#' genome(db)
#' nrow(meta(db))
make_example_database <- function(genome = "hg38") {
  marks <- base::list(
    h3k4me1 = GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges   = IRanges::IRanges(
        start = base::c(1000L, 5000L, 10000L, 20000L, 50000L),
        end   = base::c(2000L, 6000L, 11000L, 21000L, 51000L)
      )
    ),
    h3k27ac = GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges   = IRanges::IRanges(
        start = base::c(1000L, 5000L, 30000L, 50000L),
        end   = base::c(2000L, 6000L, 31000L, 51000L)
      )
    ),
    h3k27me3 = GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges   = IRanges::IRanges(
        start = base::c(10000L, 40000L),
        end   = base::c(11000L, 41000L)
      )
    ),
    h3k4me3 = GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges   = IRanges::IRanges(
        start = base::c(1000L, 50000L, 60000L),
        end   = base::c(2000L, 51000L, 61000L)
      )
    ),
    h3k36me3 = GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges   = IRanges::IRanges(
        start = base::c(70000L, 80000L),
        end   = base::c(71000L, 81000L)
      )
    )
  )

  chip <- base::list(
    TF1 = GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges   = IRanges::IRanges(
        start = base::c(1000L, 5000L, 50000L),
        end   = base::c(2000L, 6000L, 51000L)
      )
    ),
    TF2 = GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges   = IRanges::IRanges(
        start = base::c(1000L, 20000L, 50000L),
        end   = base::c(2000L, 21000L, 51000L)
      )
    )
  )

  functional <- base::list(
    fantom = GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges   = IRanges::IRanges(
        start = base::c(1000L, 50000L),
        end   = base::c(2000L, 51000L)
      )
    )
  )

  .build_epiRomicsS4_from_lists(marks, genome = genome,
                                  chip_list = chip,
                                  functional_list = functional)
}


# ---------------------------------------------------------------------------
# 2. make_example_putative_enhancers
# ---------------------------------------------------------------------------

#' Build a synthetic putative-enhancer result for use in examples
#'
#' Returns a \code{data.frame} matching the output format of
#' \code{\link{find_putative_enhancers}}. When \code{database} is \code{NULL}
#' (the default), a fresh synthetic database is created internally via
#' \code{\link{make_example_database}} and then
#' \code{\link{find_putative_enhancers}} is called on it to produce a genuine
#' result. This keeps the helper fast (< 2 s) while ensuring the output
#' structure is always consistent with the real function.
#'
#' @param database An \code{epiRomicsS4} object, or \code{NULL} (default). When
#'   \code{NULL}, a synthetic database is built automatically.
#' @return A \code{data.frame} with columns produced by
#'   \code{\link{find_putative_enhancers}}: \code{putative_id}, \code{chr},
#'   \code{start}, \code{end}, \code{width}, \code{source},
#'   \code{chromatin_state}, \code{chromatin_state_detail},
#'   \code{histone_marks}, \code{n_histone_marks}, \code{h2az},
#'   \code{tf_names}, \code{n_tfs}.
#' @family synthetic example data helpers
#' @seealso \code{\link{make_example_database}},
#'   \code{\link{find_putative_enhancers}}
#' @export
#' @examples
#' pe <- make_example_putative_enhancers()
#' nrow(pe)
#' head(pe[, c("chr", "start", "end", "chromatin_state")])
make_example_putative_enhancers <- function(database = NULL) {
  if (base::is.null(database)) {
    database <- make_example_database()
  }
  find_putative_enhancers(database)
}


# ---------------------------------------------------------------------------
# 3. make_example_enhanceosome
# ---------------------------------------------------------------------------

#' Build a synthetic enhanceosome result for use in examples
#'
#' Returns an \code{epiRomicsS4} object whose \code{annotations} slot holds a
#' synthetic enhanceosome \code{GRanges} with the same column structure as the
#' output of \code{\link{find_enhanceosomes}} (one integer count column per
#' ChIP TF in the \code{\link{meta}()} table of \code{database}, plus a
#' \code{ChIP_Hits} total column).
#'
#' The result is built directly from in-memory structures (no ChIPseeker
#' annotation, no TxDb lookup) so the example completes in under one second
#' under \command{R CMD check}. It is fit-for-purpose for exercising the
#' enhanceosome-consuming functions
#' (\code{\link{analyze_tf_cobinding}}, \code{\link{analyze_tf_overlap}}, etc.)
#' without the cost of the full enhancer-calling pipeline.
#'
#' @param database An \code{epiRomicsS4} object, or \code{NULL} (default).
#'   When \code{NULL}, \code{\link{make_example_database}} is called
#'   internally.
#' @return An \code{epiRomicsS4} object whose \code{annotations} slot contains
#'   the synthetic enhanceosome GRanges.
#' @family synthetic example data helpers
#' @seealso \code{\link{make_example_database}},
#'   \code{\link{find_enhanceosomes}}
#' @export
#' @examples
#' eso <- make_example_enhanceosome()
#' length(annotations(eso))
make_example_enhanceosome <- function(database = NULL) {
  if (base::is.null(database)) {
    database <- make_example_database()
  }

  ## Build the enhanceosome GRanges synthetically (no ChIPseeker annotation,
  ## no TxDb lookup) so the example runs in < 1 s under R CMD check. The
  ## ranges, ChIP_Hits column, and per-TF count columns mirror the output
  ## format of find_enhanceosomes() so downstream consumers (analyze_tf_*,
  ## filter_enhancers, etc.) work without modification.
  genome <- genome(database)
  meta   <- meta(database)
  tf_names <- meta[meta$type == "chip", "name"]

  ranges <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges   = IRanges::IRanges(
      start = base::c(1000L, 5000L, 50000L),
      end   = base::c(2000L, 6000L, 51000L)
    )
  )

  ## Deterministic synthetic ChIP hit counts matching find_enhanceosomes()
  ## output shape: one integer column per TF plus a ChIP_Hits total.
  if (base::length(tf_names) > 0L) {
    counts <- base::matrix(
      base::c(1L, 1L, 1L, 1L, 0L, 1L),
      nrow = base::length(ranges),
      ncol = base::length(tf_names),
      byrow = FALSE
    )[, base::seq_along(tf_names), drop = FALSE]
    base::colnames(counts) <- tf_names
    for (j in base::seq_along(tf_names)) {
      GenomicRanges::mcols(ranges)[[tf_names[j]]] <- counts[, j]
    }
    GenomicRanges::mcols(ranges)$ChIP_Hits <- base::rowSums(counts)
  } else {
    GenomicRanges::mcols(ranges)$ChIP_Hits <-
      base::integer(base::length(ranges))
  }
  base::names(ranges) <- base::as.character(base::seq_len(base::length(ranges)))

  methods::new(
    "epiRomicsS4",
    annotations = ranges,
    meta        = meta,
    genome      = genome,
    txdb        = txdb(database),
    organism    = organism(database)
  )
}


# ---------------------------------------------------------------------------
# 4. make_example_bigwig
# ---------------------------------------------------------------------------

#' Create a temporary synthetic BigWig file for use in examples
#'
#' Writes a BigWig file from synthetic \code{GRanges} and numeric scores.
#' The file is created in \code{base::tempdir()} by default.
#'
#' \strong{Caller responsibility}: The returned path points to a temporary file
#' that persists for the R session. Call \code{file.remove(path)} when the
#' file is no longer needed to avoid accumulating temporary files.
#'
#' @param regions_gr A \code{GRanges} object.  When \code{NULL} (default), a
#'   small five-region GRanges on chr1 is used.
#' @param scores A numeric vector the same length as \code{regions_gr}.  When
#'   \code{NULL} (default), scores 1--5 are used.
#' @param path Character file path for output.  When \code{NULL} (default), a
#'   path is generated via \code{base::tempfile(fileext = ".bw")}.
#' @return Character string giving the path to the created BigWig file.
#' @family synthetic example data helpers
#' @seealso \code{\link{make_example_database}},
#'   \code{rtracklayer::export}
#' @export
#' @examples
#' bw_path <- make_example_bigwig()
#' file.exists(bw_path)
#' file.remove(bw_path)
make_example_bigwig <- function(regions_gr = NULL,
                                 scores     = NULL,
                                 path       = NULL) {
  if (base::is.null(regions_gr)) {
    regions_gr <- GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges   = IRanges::IRanges(
        start = base::c(1000L, 5000L, 10000L, 20000L, 50000L),
        end   = base::c(2000L, 6000L, 11000L, 21000L, 51000L)
      )
    )
  }
  if (base::is.null(scores)) {
    scores <- base::seq_len(base::length(regions_gr))
  }
  if (base::is.null(path)) {
    path <- base::tempfile(fileext = ".bw")
  }

  regions_gr$score <- base::as.numeric(scores)

  si <- GenomeInfoDb::Seqinfo(
    seqnames   = base::as.character(
      base::unique(GenomicRanges::seqnames(regions_gr))
    ),
    seqlengths = base::rep(
      250000000L,
      base::length(base::unique(GenomicRanges::seqnames(regions_gr)))
    )
  )
  GenomeInfoDb::seqinfo(regions_gr) <- si

  rtracklayer::export(regions_gr, path, format = "BigWig")
  base::return(path)
}
