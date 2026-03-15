#' Resolve TxDb object from a package::object string
#'
#' Parses a TxDb specification string (e.g., "TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene")
#' and returns the actual TxDb object using getExportedValue.
#'
#' @param txdb_string character string of TxDb specification
#' @return TxDb object
#' @noRd
resolve_txdb <- function(txdb_string) {
  parts <- base::strsplit(txdb_string, "::")[[1]]
  if (base::length(parts) == 2) {
    base::getExportedValue(parts[1], parts[2])
  } else {
    base::stop(base::sprintf(
      "Invalid TxDb format: '%s'. Expected 'package::object' format (e.g., 'TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene')",
      txdb_string
    ))
  }
}

#' Calculate maximum coverage from a BigWig file for given genomic regions
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `maxCovBw()` is a legacy utility with no internal callers. Use
#' \code{.fast_bw_signal(bw, gr, type = "max")} for batch-efficient
#' BigWig max-score queries, or \code{\link{maxCovBwCached}} for
#' repeated multi-region access with optional RDS caching.
#'
#' Returns the maximum signal score from a BigWig file overlapping the given
#' genomic regions. Uses region-specific import via BigWigSelection for
#' efficiency.
#'
#' @param bw Character string path to BigWig file
#' @param gr GenomicRanges object containing regions to evaluate
#' @return Numeric value representing the maximum coverage across the regions
#' @export
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1000, 2000))
#' \donttest{
#' gr <- GenomicRanges::GRanges(
#'   seqnames = "chr1",
#'   ranges = IRanges::IRanges(start = 1000, end = 2000)
#' )
#' max_cov <- maxCovBw("path/to/file.bw", gr)
#' }
maxCovBw <- function(bw, gr) {
  .Deprecated("maxCovBwCached",
    msg = "maxCovBw() is deprecated. Use maxCovBwCached() or .fast_bw_signal() instead."
  )
  ovlp <- IRanges::subsetByOverlaps(
    rtracklayer::import(bw, format = "BigWig", selection = rtracklayer::BigWigSelection(gr)),
    gr
  )
  if (base::length(ovlp) > 0) {
    max_cov <- base::max(ovlp$score)
  } else {
    base::message("The selected genomic region has no coverage value in the BigWig. Coverage set to 0.")
    max_cov <- 0
  }
  return(max_cov)
}

#' Calculate maximum coverage from multiple BigWig files for given genomic regions
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `maxCovFiles()` is a legacy utility with no internal callers. Use
#' \code{\link{maxCovFilesCached}} with \code{fast = TRUE} for batch-efficient
#' multi-file BigWig max-score queries.
#'
#' Iterates over each genomic range and computes the maximum BigWig signal
#' across all provided files. Returns the GRanges object with max coverage
#' values as metadata.
#'
#' @param bws Character vector of paths to BigWig files
#' @param gr GenomicRanges object containing regions to evaluate
#' @param fast Logical. When \code{TRUE}, uses the BigWig R-tree index
#'   (\code{.fast_bw_signal}) for faster batch queries — approximately
#'   5-10x faster on multi-region inputs. When \code{FALSE} (default),
#'   uses the original per-region import for exact regression-baseline
#'   compatibility.
#' @return GenomicRanges object with coverage values added as metadata column X
#' @export
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1000, 2000))
#' \donttest{
#' gr <- GenomicRanges::GRanges(
#'   seqnames = "chr1",
#'   ranges = IRanges::IRanges(start = 1000, end = 2000)
#' )
#' bw_files <- c("path/to/file1.bw", "path/to/file2.bw")
#' result <- maxCovFiles(bw_files, gr)
#' # Fast mode for many regions:
#' result_fast <- maxCovFiles(bw_files, gr, fast = TRUE)
#' }
maxCovFiles <- function(bws, gr, fast = FALSE) {
  .Deprecated("maxCovFilesCached",
    msg = "maxCovFiles() is deprecated. Use maxCovFilesCached(fast = TRUE) instead."
  )
  n <- base::length(gr)
  if (fast) {
    # Batch query: use BigWig R-tree index for all regions at once per file,
    # then take element-wise max across files
    signal_mat <- base::vapply(bws, function(f) {
      .fast_bw_signal(f, gr, type = "max")
    }, base::numeric(n))
    if (base::is.matrix(signal_mat)) {
      max_cov <- base::apply(signal_mat, 1L, base::max)
    } else {
      # Single file: vapply returns a vector, not a matrix
      max_cov <- signal_mat
    }
    max_cov <- base::round(max_cov, 2)
  } else {
    # Original per-region loop (exact regression baseline values)
    max_cov <- base::numeric(n)
    for (i in base::seq_len(n)) {
      my_feat <- gr[i, ]
      max_cov[i] <- base::round(
        base::max(
          base::vapply(bws, maxCovBw, gr = my_feat, numeric(1))
        ),
        2
      )
    }
  }
  GenomicRanges::values(gr) <- max_cov
  return(gr)
}

#' Cache BigWig data to RDS for repeated multi-region access
#'
#' Imports the full BigWig file and saves as RDS for repeated access.
#' Only beneficial when querying many regions (e.g., heatmaps over thousands
#' of enhancers). For single-region queries, direct BigWig import with
#' BigWigSelection is faster.
#'
#' @param bw_path Character string path to BigWig file
#' @param cache_dir Directory to store cached data (default: tempdir())
#' @param force_refresh Logical, whether to force refresh the cache
#' @return Path to cached RDS file
#' @noRd
cache_bigwig <- function(bw_path, cache_dir = base::tempdir(), force_refresh = FALSE) {
  validate_file_paths(bw_path, "cache_bigwig")

  if (!base::dir.exists(cache_dir)) {
    base::dir.create(cache_dir, recursive = TRUE)
  }

  bw_hash <- digest::digest(bw_path, algo = "md5")
  cache_file <- base::file.path(cache_dir, base::paste0("bw_cache_", bw_hash, ".rds"))

  if (!force_refresh && base::file.exists(cache_file)) {
    cache_time <- base::file.info(cache_file)$mtime
    if (base::difftime(base::Sys.time(), cache_time, units = "days") < 7) {
      return(cache_file)
    }
  }

  base::message(base::sprintf("Caching BigWig data for: %s", base::basename(bw_path)))
  bw_data <- rtracklayer::import(bw_path, format = "BigWig")
  base::saveRDS(bw_data, cache_file)
  return(cache_file)
}

#' Load cached BigWig data from RDS
#'
#' @param cache_file Path to cached RDS file
#' @param gr GenomicRanges object for subsetting (optional)
#' @return GenomicRanges object with BigWig data
#' @noRd
load_cached_bigwig <- function(cache_file, gr = NULL) {
  validate_file_paths(cache_file, "load_cached_bigwig")

  bw_data <- base::readRDS(cache_file)
  if (!base::is.null(gr)) {
    bw_data <- IRanges::subsetByOverlaps(bw_data, gr)
  }
  return(bw_data)
}

#' Cached BigWig coverage calculation (opt-in fast path)
#'
#' Alternative to \code{\link{maxCovBw}} that supports RDS caching for
#' repeated multi-region queries. Uses \code{na.rm = TRUE} for robustness
#' with potentially missing data. Note: may produce slightly different ylim
#' values than the default \code{maxCovBw} if BigWig files contain NAs.
#'
#' @param bw_path Character string path to BigWig file
#' @param gr GenomicRanges object containing regions to evaluate
#' @param use_cache Logical, whether to use RDS caching (default: FALSE).
#'   When TRUE, imports full BigWig to RDS on first call, then reads from cache.
#'   When FALSE, uses region-specific BigWigSelection import (same speed as maxCovBw).
#' @param cache_dir Directory to store cached data (default: tempdir())
#' @return Numeric value representing the maximum coverage across the regions
#' @export
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1000, 2000))
#' \donttest{
#' gr <- GenomicRanges::GRanges(
#'   seqnames = "chr1",
#'   ranges = IRanges::IRanges(start = 1000, end = 2000)
#' )
#' max_cov <- maxCovBwCached("path/to/file.bw", gr)
#' }
maxCovBwCached <- function(bw_path, gr, use_cache = FALSE, cache_dir = base::tempdir()) {
  validate_file_paths(bw_path, "maxCovBwCached")
  validate_genomic_ranges(gr, "maxCovBwCached")

  if (use_cache) {
    cache_file <- cache_bigwig(bw_path, cache_dir)
    bw_data <- load_cached_bigwig(cache_file, gr)
  } else {
    bw_data <- rtracklayer::import(bw_path, format = "BigWig",
                                    selection = rtracklayer::BigWigSelection(gr))
    bw_data <- IRanges::subsetByOverlaps(bw_data, gr)
  }

  if (base::length(bw_data) > 0) {
    max_cov <- base::max(bw_data$score, na.rm = TRUE)
  } else {
    base::warning("The selected genomic region has no coverage value in the BigWig")
    max_cov <- 0
  }
  return(max_cov)
}

#' Cached multiple BigWig coverage calculation (opt-in fast path)
#'
#' Alternative to \code{\link{maxCovFiles}} that supports RDS caching and
#' parallel processing for multi-region batch queries (e.g., heatmaps).
#' Default behavior (use_cache=FALSE) uses region-specific BigWig import.
#'
#' @param bw_paths Character vector of paths to BigWig files
#' @param gr GenomicRanges object containing regions to evaluate
#' @param use_cache Logical, whether to use RDS caching (default: FALSE)
#' @param cache_dir Directory to store cached data (default: tempdir())
#' @param parallel Logical, whether to use parallel processing (default: FALSE)
#' @param fast Logical, whether to use batch BigWig R-tree query path for all
#'   regions at once per file (default: FALSE). When TRUE, uses
#'   \code{.fast_bw_signal()} for efficient multi-region queries.
#' @return GenomicRanges object with coverage values added as metadata column X
#' @export
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1000, 2000))
#' \donttest{
#' gr <- GenomicRanges::GRanges(
#'   seqnames = "chr1",
#'   ranges = IRanges::IRanges(start = 1000, end = 2000)
#' )
#' bw_files <- c("path/to/file1.bw", "path/to/file2.bw")
#' result <- maxCovFilesCached(bw_files, gr)
#' }
maxCovFilesCached <- function(bw_paths, gr, use_cache = FALSE, cache_dir = base::tempdir(),
                              parallel = FALSE, fast = FALSE) {
  validate_file_paths(bw_paths, "maxCovFilesCached")
  validate_genomic_ranges(gr, "maxCovFilesCached")

  n_regions <- base::length(gr)

  if (fast) {
    # Batch query: use BigWig R-tree index for all regions at once per file,
    # then take element-wise max across files (same pattern as maxCovFiles fast path)
    signal_mat <- base::vapply(bw_paths, function(f) {
      .fast_bw_signal(f, gr, type = "max")
    }, base::numeric(n_regions))
    if (base::is.matrix(signal_mat)) {
      max_cov <- base::apply(signal_mat, 1L, base::max)
    } else {
      max_cov <- signal_mat
    }
    max_cov <- base::round(max_cov, 2)
  } else if (parallel && base::requireNamespace("parallel", quietly = TRUE)) {
    n_cores <- base::max(1L, parallel::detectCores() - 1L)
    max_cov <- base::unlist(parallel::mclapply(base::seq_len(n_regions), function(i) {
      my_feat <- gr[i, ]
      base::round(
        base::max(base::vapply(bw_paths, function(bw) {
          maxCovBwCached(bw, my_feat, use_cache, cache_dir)
        }, numeric(1))),
        2
      )
    }, mc.cores = n_cores))
  } else {
    max_cov <- base::numeric(n_regions)
    for (i in base::seq_len(n_regions)) {
      my_feat <- gr[i, ]
      max_cov[i] <- base::round(
        base::max(base::vapply(bw_paths, function(bw) {
          maxCovBwCached(bw, my_feat, use_cache, cache_dir)
        }, numeric(1))),
        2
      )
    }
  }

  GenomicRanges::values(gr) <- max_cov
  return(gr)
}

# --- BigWig pre-import for fast visualization ---

#' Pre-import BigWig files and build Gviz DataTrack objects from in-memory data
#'
#' Imports each BigWig file once for a viewing region, computes the shared
#' y-axis limit (max coverage), and returns DataTrack objects constructed from
#' in-memory GRanges. This eliminates the double-read bottleneck where
#' maxCovFiles reads each file and then Gviz re-reads at plotTracks time.
#'
#' @param bw_paths character vector of BigWig file paths
#' @param region GRanges viewing region (single range)
#' @param sample_names character vector of sample names for track labels
#' @param colors character vector of colors for each track
#' @param genome character genome name (e.g., "hg38")
#' @param chr character chromosome name
#' @return list with components:
#'   \describe{
#'     \item{tracks}{list of Gviz DataTrack objects}
#'     \item{max_cov}{numeric, computed ylim value (ceiling-rounded)}
#'   }
#' @noRd
.build_data_tracks <- function(bw_paths, region, sample_names, colors, genome, chr) {
  if (!base::requireNamespace("Gviz", quietly = TRUE)) {
    base::stop("Package 'Gviz' is required for .build_data_tracks(). ",
               "Install it with: BiocManager::install('Gviz')", call. = FALSE)
  }
  n_files <- base::length(bw_paths)

  # Fast max score computation via R-tree index (no full import needed)
  max_scores <- .fast_bw_signal(bw_paths[1L], region, type = "max")
  if (n_files > 1L) {
    max_scores <- base::vapply(bw_paths, function(f) {
      .fast_bw_signal(f, region, type = "max")
    }, numeric(1L))
  } else {
    max_scores <- max_scores[1L]
  }

  # Compute ylim (ceiling((max + 10%) / 10) * 10)
  overall_max <- base::max(max_scores)
  max_cov <- base::ceiling((overall_max + (overall_max * 0.1)) / 10) * 10

  # Import raw GRanges data for Gviz DataTrack (needs position-level data)
  imported_data <- base::lapply(bw_paths, function(f) {
    bw_data <- rtracklayer::import(f, format = "BigWig",
      selection = rtracklayer::BigWigSelection(region))
    IRanges::subsetByOverlaps(bw_data, region)
  })

  # Build DataTrack objects from in-memory GRanges
  tracks <- base::lapply(base::seq_len(n_files), function(j) {
    Gviz::DataTrack(
      range = imported_data[[j]],
      size = 10,
      col = colors[j],
      fill = colors[j],
      col.histogram = colors[j],
      ylim = base::c(0, max_cov),
      genome = genome,
      type = "hist",
      chromosome = chr,
      name = sample_names[j]
    )
  })

  return(base::list(tracks = tracks, max_cov = max_cov))
}

# --- IdeogramTrack caching ---
# Package-level environment for caching IdeogramTrack objects
.epiRomics_cache <- base::new.env(parent = base::emptyenv())

#' Get a cached IdeogramTrack (avoids repeated UCSC network calls)
#'
#' IdeogramTrack downloads cytoband data from UCSC on every call (~3-4s).
#' This function caches the result per genome+chromosome combination.
#'
#' @param genome Character string genome name (e.g., "hg38", "mm10")
#' @param chromosome Character string chromosome name (e.g., "chr1")
#' @return Gviz IdeogramTrack object
#' @noRd
get_ideogram_track <- function(genome, chromosome) {
  if (!base::requireNamespace("Gviz", quietly = TRUE)) {
    base::stop("Package 'Gviz' is required for get_ideogram_track(). ",
               "Install it with: BiocManager::install('Gviz')", call. = FALSE)
  }
  cache_key <- base::paste0("ideogram_", genome, "_", chromosome)
  if (base::exists(cache_key, envir = .epiRomics_cache)) {
    return(base::get(cache_key, envir = .epiRomics_cache))
  }
  itrack <- Gviz::IdeogramTrack(genome = genome, chromosome = chromosome)
  base::assign(cache_key, itrack, envir = .epiRomics_cache)
  return(itrack)
}

#' Validate epiRomics database object
#'
#' @param epiRomics_dB epiRomics database object to validate
#' @param function_name Name of the calling function for error messages
#' @return TRUE if valid, stops with error if invalid
#' @noRd
validate_epiRomics_dB <- function(epiRomics_dB, function_name = "function") {
  if (!methods::is(epiRomics_dB, "epiRomicsS4")) {
    base::stop(base::sprintf("%s: epiRomics_dB must be an epiRomicsS4 object", function_name))
  }
  if (base::is.null(epiRomics_dB@annotations) || base::length(epiRomics_dB@annotations) == 0) {
    base::stop(base::sprintf("%s: epiRomics_dB@annotations is empty or NULL", function_name))
  }
  if (base::is.null(epiRomics_dB@meta) || base::nrow(epiRomics_dB@meta) == 0) {
    base::stop(base::sprintf("%s: epiRomics_dB@meta is empty or NULL", function_name))
  }
  if (base::is.null(epiRomics_dB@genome) || base::length(epiRomics_dB@genome) == 0) {
    base::stop(base::sprintf("%s: epiRomics_dB@genome is empty or NULL", function_name))
  }
  return(TRUE)
}

#' Validate genomic ranges object
#'
#' @param gr GenomicRanges object to validate
#' @param function_name Name of the calling function for error messages
#' @return TRUE if valid, stops with error if invalid
#' @noRd
validate_genomic_ranges <- function(gr, function_name = "function") {
  if (!methods::is(gr, "GRanges")) {
    base::stop(base::sprintf("%s: gr must be a GRanges object", function_name))
  }
  if (base::length(gr) == 0) {
    base::stop(base::sprintf("%s: gr is empty", function_name))
  }
  return(TRUE)
}

#' Validate file paths
#'
#' @param file_paths Character vector of file paths to validate
#' @param function_name Name of the calling function for error messages
#' @param check_exists Whether to check if files exist
#' @return TRUE if valid, stops with error if invalid
#' @noRd
validate_file_paths <- function(file_paths, function_name = "function", check_exists = TRUE) {
  if (!base::is.character(file_paths)) {
    base::stop(base::sprintf("%s: file_paths must be a character vector", function_name))
  }
  if (base::length(file_paths) == 0) {
    base::stop(base::sprintf("%s: file_paths is empty", function_name))
  }
  if (check_exists) {
    missing_files <- file_paths[!base::file.exists(file_paths)]
    if (base::length(missing_files) > 0) {
      base::stop(base::sprintf(
        "%s: The following files do not exist: %s",
        function_name, base::paste(missing_files, collapse = ", ")
      ))
    }
  }
  return(TRUE)
}

#' Validate character parameters
#'
#' @param param Character parameter to validate
#' @param param_name Name of the parameter for error messages
#' @param function_name Name of the calling function for error messages
#' @param allow_empty Whether to allow empty strings
#' @return TRUE if valid, stops with error if invalid
#' @noRd
validate_character_param <- function(param, param_name, function_name = "function", allow_empty = FALSE) {
  if (!base::is.character(param)) {
    base::stop(base::sprintf("%s: %s must be a character string", function_name, param_name))
  }
  if (base::length(param) != 1) {
    base::stop(base::sprintf("%s: %s must be a single character string", function_name, param_name))
  }
  if (!allow_empty && (base::is.na(param) || param == "")) {
    base::stop(base::sprintf("%s: %s cannot be empty or NA", function_name, param_name))
  }
  return(TRUE)
}

#' Validate numeric parameters
#'
#' @param param Numeric parameter to validate
#' @param param_name Name of the parameter for error messages
#' @param function_name Name of the calling function for error messages
#' @param min_val Minimum allowed value (optional)
#' @param max_val Maximum allowed value (optional)
#' @return TRUE if valid, stops with error if invalid
#' @noRd
validate_numeric_param <- function(param, param_name, function_name = "function", min_val = NULL, max_val = NULL) {
  if (!base::is.numeric(param)) {
    base::stop(base::sprintf("%s: %s must be numeric", function_name, param_name))
  }
  if (base::length(param) != 1) {
    base::stop(base::sprintf("%s: %s must be a single numeric value", function_name, param_name))
  }
  if (base::is.na(param)) {
    base::stop(base::sprintf("%s: %s cannot be NA", function_name, param_name))
  }
  if (!base::is.null(min_val) && param < min_val) {
    base::stop(base::sprintf("%s: %s must be >= %s", function_name, param_name, min_val))
  }
  if (!base::is.null(max_val) && param > max_val) {
    base::stop(base::sprintf("%s: %s must be <= %s", function_name, param_name, max_val))
  }
  return(TRUE)
}

#' Chromatin state color palette for visualization
#'
#' Returns a named character vector of colors for chromatin state classification.
#' Used by visualization functions and track layers.
#'
#' @return named character vector of hex color codes
#' @noRd
.chromatin_state_palette <- function() {
  c(
    # 6 simplified chromatin states
    active    = "#2ECC40",      # Green -- H3K27ac-containing (active)
    poised    = "#FF851B",      # Orange -- H3K4me1 + H3K27me3
    repressed = "#E74C3C",      # Red -- H3K27me3/H3K9me3 silencing
    bivalent  = "#B10DC9",      # Purple -- H3K4me3 + H3K27me3
    primed    = "#FFDC00",      # Yellow -- H3K4me1 only
    unmarked  = "#D3D3D3",      # Light gray -- no marks
    # Special tracks (super-enhancer overlay)
    super_enhancer     = "#008080",      # Dark teal
    constituent_peak   = "#CCCCCC"       # Medium gray
  )
}

#' Fast BigWig mean signal extraction using R-tree index
#'
#' Uses \code{rtracklayer::summary()} which queries the BigWig R-tree zoom
#' levels directly, computing mean signal in O(log n) per region. This is
#' dramatically faster than import() + findOverlaps() + tapply() which reads
#' all underlying intervals. Falls back to import-based aggregation if
#' summary() fails (e.g., malformed BigWig or missing zoom levels).
#'
#' @param bw_path character path to BigWig file
#' @param regions GRanges of target regions
#' @param type summary statistic: "mean" (default), "max", "min"
#' @return numeric vector of length \code{length(regions)} with signal;
#'   regions with no signal return 0
#' @noRd
.fast_bw_signal <- function(bw_path, regions, type = "mean") {
  n <- base::length(regions)
  if (n == 0L) return(base::numeric(0L))
  base::tryCatch({
    bwf <- rtracklayer::BigWigFile(bw_path)
    s <- rtracklayer::summary(bwf, which = regions, size = 1L, type = type)
    # Vectorized extraction: unlist avoids per-element R function dispatch
    vals <- base::as.numeric(base::unlist(s, use.names = FALSE))
    vals[base::is.na(vals)] <- 0
    vals
  }, error = function(e) {
    # Fallback: import + aggregate (slower but robust)
    bw_data <- rtracklayer::import(bw_path, format = "BigWig",
      selection = rtracklayer::BigWigSelection(regions))
    agg_fun <- base::switch(type,
      "mean" = base::mean, "max" = base::max, "min" = base::min,
      base::mean)
    .aggregate_signal_over_regions(bw_data, regions, agg_fun = agg_fun)
  })
}

#' Fast multi-replicate BigWig mean signal
#'
#' Computes mean signal across multiple BigWig replicates without for-loops.
#' Uses \code{.fast_bw_signal()} for each file via \code{vapply}, then
#' averages with \code{rowMeans}. Replaces for-loop + import + aggregate
#' patterns for replicate averaging.
#'
#' @param bw_paths character vector of BigWig file paths (replicates)
#' @param regions GRanges of target regions
#' @param type summary statistic: "mean" (default), "max"
#' @return numeric vector of mean signal across replicates
#' @noRd
.fast_bw_replicate_mean <- function(bw_paths, regions, type = "mean") {
  n_regions <- base::length(regions)
  n_files <- base::length(bw_paths)
  if (n_regions == 0L || n_files == 0L) return(base::numeric(n_regions))
  if (n_files == 1L) {
    return(.fast_bw_signal(bw_paths[1L], regions, type = type))
  }
  # Build matrix: rows=regions, cols=replicates -- no for-loop
  signal_mat <- base::vapply(bw_paths, function(f) {
    .fast_bw_signal(f, regions, type = type)
  }, base::numeric(n_regions))
  base::rowMeans(signal_mat, na.rm = TRUE)
}

#' Aggregate BigWig signal over genomic regions (vectorized)
#'
#' Replaces per-region subsetByOverlaps loops with a single findOverlaps +
#' tapply call for O(n) instead of O(n^2) performance. This is the core
#' speed-critical helper used by signal enrichment, super-enhancer
#' quantification, ABC scoring, and heatmap binning. Also serves as fallback
#' for \code{.fast_bw_signal()} when BigWig summary() is unavailable.
#'
#' @param bw_data GRanges with \code{score} metadata column (from BigWig import)
#' @param regions GRanges of target regions to aggregate over
#' @param agg_fun aggregation function applied per region (default: base::mean)
#' @return numeric vector of length \code{length(regions)} with aggregated
#'   signal; regions with no overlapping signal return 0
#' @noRd
.aggregate_signal_over_regions <- function(bw_data, regions,
                                            agg_fun = base::mean) {
  signal <- base::numeric(base::length(regions))
  if (base::length(bw_data) == 0L || base::length(regions) == 0L) {
    return(signal)
  }
  hits <- GenomicRanges::findOverlaps(regions, bw_data)
  if (base::length(hits) == 0L) {
    return(signal)
  }
  q_idx <- S4Vectors::queryHits(hits)
  s_idx <- S4Vectors::subjectHits(hits)
  # data.table group-by replaces tapply for speed on large overlap sets
  dt_sig <- data.table::data.table(q = q_idx, score = bw_data$score[s_idx])
  if (base::identical(agg_fun, base::mean)) {
    dt_agg <- dt_sig[, base::list(val = base::mean(score)), by = q]
  } else if (base::identical(agg_fun, base::sum)) {
    dt_agg <- dt_sig[, base::list(val = base::sum(score)), by = q]
  } else if (base::identical(agg_fun, base::max)) {
    dt_agg <- dt_sig[, base::list(val = base::max(score)), by = q]
  } else {
    dt_agg <- dt_sig[, base::list(val = agg_fun(score)), by = q]
  }
  signal[dt_agg$q] <- dt_agg$val
  signal
}

#' Fast BigWig weighted signal using R-tree index
#'
#' Computes total signal (mean * region_width) per region using the fast
#' \code{rtracklayer::summary()} path. For BigWig format, mean * width gives
#' the exact weighted sum since mean = sum(score_i * interval_width_i) /
#' region_width. Falls back to import-based weighted aggregation on error.
#'
#' @param bw_path character path to BigWig file
#' @param regions GRanges of target regions
#' @return numeric vector of weighted signal sums per region
#' @noRd
.fast_bw_weighted_signal <- function(bw_path, regions) {
  n <- base::length(regions)
  if (n == 0L) return(base::numeric(0L))
  base::tryCatch({
    mean_signal <- .fast_bw_signal(bw_path, regions, type = "mean")
    widths <- BiocGenerics::width(regions)
    mean_signal * base::as.numeric(widths)
  }, error = function(e) {
    bw_data <- rtracklayer::import(bw_path, format = "BigWig",
      selection = rtracklayer::BigWigSelection(regions))
    .aggregate_weighted_signal(bw_data, regions)
  })
}

#' Aggregate weighted BigWig signal (score * width) over genomic regions
#'
#' Variant of \code{.aggregate_signal_over_regions} that computes the weighted
#' sum of \code{score * width} per region, as required by the ROSE
#' super-enhancer algorithm for total signal quantification.
#' Also serves as fallback for \code{.fast_bw_weighted_signal()}.
#'
#' @param bw_data GRanges with \code{score} metadata column
#' @param regions GRanges of target regions
#' @return numeric vector of weighted signal sums per region
#' @noRd
.aggregate_weighted_signal <- function(bw_data, regions) {
  signal <- base::numeric(base::length(regions))
  if (base::length(bw_data) == 0L || base::length(regions) == 0L) {
    return(signal)
  }
  hits <- GenomicRanges::findOverlaps(regions, bw_data)
  if (base::length(hits) == 0L) {
    return(signal)
  }
  q_idx <- S4Vectors::queryHits(hits)
  s_idx <- S4Vectors::subjectHits(hits)
  weighted_scores <- bw_data$score[s_idx] * BiocGenerics::width(bw_data)[s_idx]
  dt_ws <- data.table::data.table(q = q_idx, ws = weighted_scores)
  dt_agg <- dt_ws[, base::list(val = base::sum(ws)), by = q]
  signal[dt_agg$q] <- dt_agg$val
  signal
}

#' Validate logical parameters
#'
#' @param param Logical parameter to validate
#' @param param_name Name of the parameter for error messages
#' @param function_name Name of the calling function for error messages
#' @return TRUE if valid, stops with error if invalid
#' @noRd
validate_logical_param <- function(param, param_name, function_name = "function") {
  if (!base::is.logical(param)) {
    base::stop(base::sprintf("%s: %s must be logical", function_name, param_name))
  }
  if (base::length(param) != 1) {
    base::stop(base::sprintf("%s: %s must be a single logical value", function_name, param_name))
  }
  if (base::is.na(param)) {
    base::stop(base::sprintf("%s: %s cannot be NA", function_name, param_name))
  }
  return(TRUE)
}
