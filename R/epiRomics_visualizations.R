#' Compute per-region signal z-scores from a BigWig file
#'
#' Calculates the mean ATAC/ChIP signal per region from a BigWig file and
#' returns z-scores relative to the distribution across all input regions.
#' The null distribution is estimated from the signal across all regions,
#' assuming most of the genome has low/no signal. Regions with z-scores above
#' the threshold are called as accessible.
#'
#' @param bw_path character, path to a BigWig file.
#' @param regions GRanges object of regions to classify.
#' @param z_threshold numeric, z-score cutoff for calling accessibility
#'   (default: 1.0). Regions with z >= threshold are called accessible.
#'   Set to NULL to return z-scores without thresholding. Ignored when
#'   \code{auto_threshold = TRUE}.
#' @param auto_threshold logical. If TRUE, automatically detects the
#'   optimal threshold by finding the valley between the noise and signal
#'   peaks in the log-transformed signal distribution using kernel density
#'   estimation. This is the recommended approach -- let the data define
#'   background vs. signal rather than using a fixed z-score (default: FALSE).
#' @param return_scores logical. If TRUE, returns a list with z_scores,
#'   raw signal, and boolean calls. If FALSE (default), returns a logical
#'   vector for backward compatibility.
#' @return If \code{return_scores = FALSE} (default): a logical vector the
#'   same length as \code{regions}, TRUE for accessible regions.
#'   If \code{return_scores = TRUE}: a list with:
#'   \describe{
#'     \item{z_scores}{numeric vector of z-scores}
#'     \item{signal}{numeric vector of raw mean signal per region}
#'     \item{accessible}{logical vector (z >= threshold)}
#'     \item{mu}{mean signal (background estimate)}
#'     \item{sigma}{SD of signal}
#'   }
#' @export
#' @examples
#' tryCatch(call_accessible_regions("nonexistent.bw", GenomicRanges::GRanges()),
#'   error = function(e) message(e$message))
#' \donttest{
#' # Simple boolean call with fixed threshold
#' accessible <- call_accessible_regions("sample.bw", regions, z_threshold = 1)
#'
#' # Auto-detect threshold from signal distribution (recommended)
#' result <- call_accessible_regions("sample.bw", regions, auto_threshold = TRUE,
#'   return_scores = TRUE)
#' cat("Auto z-threshold:", result$auto_z, "\n")
#'
#' # Get full z-score information for manual threshold selection
#' result <- call_accessible_regions("sample.bw", regions, return_scores = TRUE)
#' hist(result$signal, breaks = 100)  # inspect distribution
#' }
call_accessible_regions <- function(bw_path, regions, z_threshold = 1.0,
                                     auto_threshold = FALSE,
                                     return_scores = FALSE) {
  if (!base::file.exists(bw_path)) {
    base::stop(base::sprintf("BigWig file not found: %s", bw_path))
  }
  if (!methods::is(regions, "GRanges") || base::length(regions) == 0) {
    base::stop("regions must be a non-empty GRanges object")
  }

  n_regions <- base::length(regions)

  # Fast BigWig signal extraction using R-tree index
  region_signal <- .fast_bw_signal(bw_path, regions, type = "mean")

  # Check if any signal was found
  if (base::all(region_signal == 0)) {
    if (return_scores) {
      return(base::list(
        z_scores = base::rep(0, n_regions),
        signal = base::rep(0, n_regions),
        accessible = base::rep(FALSE, n_regions),
        mu = 0, sigma = 0))
    }
    return(base::rep(FALSE, n_regions))
  }

  # Z-score: (x - mean) / sd across all regions
  mu <- base::mean(region_signal)
  sigma <- stats::sd(region_signal)
  if (sigma == 0 || base::is.na(sigma)) {
    if (return_scores) {
      return(base::list(
        z_scores = base::rep(0, n_regions),
        signal = region_signal,
        accessible = base::rep(FALSE, n_regions),
        mu = mu, sigma = 0))
    }
    return(base::rep(FALSE, n_regions))
  }

  z_scores <- (region_signal - mu) / sigma

  # Auto-detect threshold from signal distribution valley
  auto_z <- NA_real_
  if (auto_threshold) {
    log_sig <- base::log10(region_signal + 1)
    dens <- .safe_density(log_sig)
    valley_x <- .find_signal_valley(dens)

    valley_raw <- if (!base::is.na(valley_x)) {
      10^valley_x - 1
    } else {
      10^0.1 - 1  # Fallback: log10 = 0.1 (raw ~ 0.26)
    }
    auto_z <- (valley_raw - mu) / sigma
    z_threshold <- auto_z
  }

  accessible <- if (!base::is.null(z_threshold)) {
    z_scores >= z_threshold
  } else {
    base::rep(TRUE, n_regions)
  }

  if (return_scores) {
    return(base::list(
      z_scores = z_scores,
      signal = region_signal,
      accessible = accessible,
      mu = mu, sigma = sigma,
      auto_z = if (auto_threshold) auto_z else NA_real_))
  }
  return(accessible)
}


#' Plot signal distribution histogram for accessibility threshold selection
#'
#' Visualizes the distribution of BigWig signal across genomic regions to
#' help identify the noise floor and select an appropriate z-score or raw
#' signal threshold for calling accessible regions. The distribution is
#' typically bimodal: a large peak at zero/low signal (background noise)
#' and a smaller peak at higher signal (accessible regions).
#'
#' @param bw_paths named character vector of BigWig file paths.
#' @param regions GRanges object of regions to examine.
#' @param n_bins integer, number of histogram bins (default: 100).
#' @param show_z_lines numeric vector, z-score thresholds to show as
#'   vertical lines (default: c(1.0, 1.5, 2.0)).
#' @param log_scale logical, use log10(signal + 1) for x-axis (default: TRUE).
#' @param cex_label Numeric. Font size for threshold and annotation labels
#'   (default 0.7).
#' @param cex_title Numeric. Font size for histogram titles (default 1.0).
#' @return Invisible list of per-sample signal statistics, including
#'   mean, sd, quantiles, and suggested thresholds.
#' @export
#' @examples
#' tryCatch(plot_signal_histogram(c(Alpha = "nonexistent.bw"),
#'   GenomicRanges::GRanges()), error = function(e) message(e$message))
#' \donttest{
#' stats <- plot_signal_histogram(
#'   bw_paths = c(Alpha = "alpha.bw", Beta = "beta.bw"),
#'   regions = enhanceosome_regions
#' )
#' # Use stats to pick threshold:
#' stats$Alpha$suggested_threshold
#' }
plot_signal_histogram <- function(bw_paths, regions, n_bins = 100,
                                   show_z_lines = c(1.0, 1.5, 2.0),
                                   log_scale = TRUE,
                                   cex_label = 0.7,
                                   cex_title = 1.0) {
  if (!base::is.character(bw_paths) || base::length(bw_paths) == 0) {
    base::stop("bw_paths must be a character vector of BigWig file paths")
  }
  if (!methods::is(regions, "GRanges") || base::length(regions) == 0) {
    base::stop("regions must be a non-empty GRanges object")
  }

  ct_names <- if (!base::is.null(base::names(bw_paths))) {
    base::names(bw_paths)
  } else {
    base::paste0("Sample_", base::seq_along(bw_paths))
  }

  n_samples <- base::length(bw_paths)
  stats_list <- base::list()

  # Set up multi-panel plot
  old_par <- graphics::par(
    mfrow = base::c(base::ceiling(n_samples / 2), base::min(n_samples, 2)),
    mar = base::c(4, 4, 3, 1))
  base::on.exit(graphics::par(old_par), add = TRUE)

  for (i in base::seq_along(bw_paths)) {
    ct <- ct_names[i]
    result <- call_accessible_regions(bw_paths[i], regions,
      z_threshold = NULL, return_scores = TRUE)

    sig <- result$signal
    mu <- result$mu
    sigma <- result$sigma

    # Log transform for visualization if requested
    plot_vals <- if (log_scale) base::log10(sig + 1) else sig
    x_lab <- if (log_scale) "log10(mean signal + 1)" else "Mean signal"

    # Histogram — cap y-axis so signal is visible (noise peak dominates)
    h <- graphics::hist(plot_vals, breaks = n_bins, plot = FALSE)
    sorted_counts <- base::sort(h$counts, decreasing = TRUE)
    y_cap <- if (base::length(sorted_counts) >= 3) {
      sorted_counts[base::max(2, base::ceiling(base::length(sorted_counts) * 0.05))] * 1.3
    } else {
      base::max(h$counts)
    }
    graphics::hist(plot_vals, breaks = n_bins, main = ct,
      xlab = x_lab, ylab = "# Regions", col = "#CCDDEE", border = "#88AACC",
      freq = TRUE, ylim = base::c(0, y_cap))

    # Indicate if noise peak was truncated
    if (base::max(h$counts) > y_cap) {
      noise_peak_x <- h$mids[base::which.max(h$counts)]
      graphics::text(noise_peak_x, y_cap * 0.95,
        labels = base::paste0("(truncated: ", base::max(h$counts), ")"),
        cex = cex_label * 0.85, col = "#888888")
    }

    # Add density overlay
    dens <- .safe_density(plot_vals)
    dens_scale <- y_cap / base::max(dens$y) * 0.8
    graphics::lines(dens$x, dens$y * dens_scale, col = "#2C3E50", lwd = 2)

    # Add z-score threshold lines
    z_colors <- base::c("#2ECC40", "#FF851B", "#E74C3C")
    for (j in base::seq_along(show_z_lines)) {
      z_val <- show_z_lines[j]
      raw_threshold <- mu + z_val * sigma
      plot_threshold <- if (log_scale) base::log10(raw_threshold + 1) else raw_threshold
      line_col <- if (j <= base::length(z_colors)) z_colors[j] else "#333333"
      graphics::abline(v = plot_threshold, col = line_col, lwd = 2, lty = 2)
      n_above <- base::sum(sig >= raw_threshold)
      graphics::text(plot_threshold, graphics::par("usr")[4] * 0.95,
        labels = base::paste0("z=", z_val, "\n(n=", n_above, ")"),
        pos = 4, cex = cex_label, col = line_col)
    }

    # Suggested threshold via shared valley detection
    min_after_peak <- .find_signal_valley(dens)

    suggested_raw <- if (!base::is.na(min_after_peak)) {
      if (log_scale) 10^min_after_peak - 1 else min_after_peak
    } else {
      10^0.1 - 1
    }

    suggested_z <- if (sigma > 0) (suggested_raw - mu) / sigma else 1.0

    stats_list[[ct]] <- base::list(
      mu = mu, sigma = sigma,
      n_regions = base::length(sig),
      n_nonzero = base::sum(sig > 0),
      quantiles = stats::quantile(sig, probs = base::c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1)),
      suggested_threshold = suggested_raw,
      suggested_z = base::round(suggested_z, 2),
      n_accessible_z1 = base::sum(result$z_scores >= 1.0),
      n_accessible_z1.5 = base::sum(result$z_scores >= 1.5),
      n_accessible_z2 = base::sum(result$z_scores >= 2.0)
    )

    # Add suggested threshold line
    if (!base::is.na(min_after_peak)) {
      graphics::abline(v = min_after_peak, col = "#B10DC9", lwd = 2, lty = 1)
      graphics::text(min_after_peak, graphics::par("usr")[4] * 0.85,
        labels = base::paste0("suggested\nz=", base::round(suggested_z, 1)),
        pos = 4, cex = cex_label, col = "#B10DC9")
    }
  }

  base::invisible(stats_list)
}


#' Classify regions by cell-type-specific binary accessibility
#'
#' For each genomic region, determines whether it is accessible in each
#' cell type using per-sample z-score thresholding on BigWig ATAC-seq signal.
#' Unlike fold-change approaches (which compare relative enrichment between
#' samples), this method makes an independent binary open/closed call per
#' sample, then classifies regions based on the combination of binary calls.
#'
#' A region with high signal in BOTH cell types is correctly classified as
#' \code{"Shared"} regardless of the magnitude difference between them.
#' This avoids the fold-change pitfall where a region with signal 50 in one
#' cell type and 100 in the other is called "enriched" even though both have
#' accessible chromatin.
#'
#' @section Methodology:
#' For each BigWig file independently:
#' \enumerate{
#'   \item Import mean ATAC signal per region
#'   \item Compute genome-wide null distribution (mean and SD of all region
#'     signals)
#'   \item Apply z-score threshold: z = (signal - mean) / sd
#'   \item Regions with z >= threshold are called accessible (open)
#' }
#' The z-score approach normalizes within each sample, making the open/closed
#' call robust to differences in sequencing depth or library complexity
#' between samples (Buenrostro et al. 2013, \emph{Nat Methods};
#' Corces et al. 2018, \emph{Science}).
#'
#' @param bw_paths named character vector of BigWig file paths. Names are
#'   cell-type labels (e.g., \code{c(Alpha = "alpha.atac.bw",
#'   Beta = "beta.atac.bw")}).
#' @param regions GRanges object of regions to classify.
#' @param z_threshold numeric, z-score cutoff for calling a region as
#'   accessible (default: 1.0). Higher values are more stringent. Use
#'   \code{\link{plot_signal_histogram}} to visualize the signal distribution
#'   and identify an appropriate threshold for your data.
#' @param auto_threshold logical. If TRUE, automatically detect the optimal
#'   threshold from the signal distribution using kernel density valley
#'   detection, overriding \code{z_threshold} (default: FALSE).
#' @return A data.frame with columns:
#'   \describe{
#'     \item{region_idx}{integer, index into input regions}
#'     \item{cell_type}{character, cell-type label (single name if
#'       accessible in only one cell type), \code{"Shared"} if accessible
#'       in all cell types, or \code{"Closed"} if not accessible in any}
#'     \item{accessible_in}{character, comma-separated list of cell types
#'       where the region is accessible (empty string if closed)}
#'   }
#'   The binary accessibility matrix (logical, regions x cell types) is
#'   attached as the \code{"accessibility_matrix"} attribute.
#' @export
#' @examples
#' tryCatch(classify_celltype_accessibility(c(Alpha = "nonexistent.bw"),
#'   GenomicRanges::GRanges()), error = function(e) message(e$message))
#' \donttest{
#' ct <- classify_celltype_accessibility(
#'   bw_paths = c(Alpha = "alpha.atac.bw", Beta = "beta.atac.bw"),
#'   regions = enhanceosome_regions
#' )
#' table(ct$cell_type)  # Shared, Alpha, Beta, Closed
#' }
classify_celltype_accessibility <- function(bw_paths, regions, z_threshold = 1.0,
                                            auto_threshold = FALSE) {
  if (!base::is.character(bw_paths) || base::is.null(base::names(bw_paths))) {
    base::stop("bw_paths must be a named character vector (names = cell-type labels)")
  }
  if (base::any(base::names(bw_paths) == "")) {
    base::stop("All elements of bw_paths must be named")
  }
  if (!methods::is(regions, "GRanges") || base::length(regions) == 0) {
    base::stop("regions must be a non-empty GRanges object")
  }

  ct_names <- base::names(bw_paths)
  n_ct <- base::length(ct_names)
  n_regions <- base::length(regions)

  # Build accessibility matrix: regions x cell_types (logical)
  acc_matrix <- base::matrix(FALSE, nrow = n_regions, ncol = n_ct)
  base::colnames(acc_matrix) <- ct_names

  for (i in base::seq_along(ct_names)) {
    acc_matrix[, i] <- call_accessible_regions(bw_paths[i], regions, z_threshold,
      auto_threshold = auto_threshold)
  }

  # Classify each region
  n_accessible <- base::rowSums(acc_matrix)

  # Build accessible_in string column-wise (iterates over cell types, not regions)
  accessible_in <- base::rep("", n_regions)
  for (j in base::seq_along(ct_names)) {
    has_ct <- acc_matrix[, j]
    accessible_in[has_ct] <- base::ifelse(
      accessible_in[has_ct] == "",
      ct_names[j],
      base::paste0(accessible_in[has_ct], ",", ct_names[j])
    )
  }

  # Cell type assignment logic
  cell_type <- base::character(n_regions)
  cell_type[n_accessible == 0] <- "Closed"
  cell_type[n_accessible == n_ct] <- "Shared"

  # Regions accessible in exactly one cell type (vapply for type safety)
  single_ct <- n_accessible == 1
  if (base::any(single_ct)) {
    cell_type[single_ct] <- base::vapply(
      base::which(single_ct),
      function(idx) ct_names[base::which(acc_matrix[idx, ])[1]],
      character(1)
    )
  }

  # Regions accessible in >1 but <all cell types (partial overlap)
  mixed <- n_accessible > 1 & n_accessible < n_ct
  if (base::any(mixed)) {
    cell_type[mixed] <- "Shared"
  }

  result <- base::data.frame(
    region_idx = base::seq_len(n_regions),
    cell_type = cell_type,
    accessible_in = accessible_in,
    stringsAsFactors = FALSE
  )
  base::attr(result, "accessibility_matrix") <- acc_matrix
  return(result)
}


# ---- Internal helpers (not exported) ----

#' Compute kernel density with SJ bandwidth, fallback to nrd0
#' @param x numeric vector of values
#' @param n integer, number of density points (default: 2048)
#' @return density object
#' @noRd
.safe_density <- function(x, n = 2048L) {
  base::tryCatch(
    stats::density(x, n = n, bw = "SJ"),
    error = function(e) stats::density(x, n = n, bw = "nrd0")
  )
}


#' Find valley between noise and signal peaks in density
#'
#' Given a density object (from stats::density), finds the first local minimum
#' after the global maximum (noise peak) where the density has dropped by at
#' least 5% of the peak height. Used for automatic threshold selection in
#' signal distributions. Values are capped at 0.5 (assumes log10 scale).
#'
#' @param dens density object from stats::density
#' @param peak_fraction_drop numeric, minimum drop from peak height to qualify
#'   as a valley (default: 0.05)
#' @param cap numeric, maximum allowed valley position (default: 0.5)
#' @param fallback numeric, default value when cap is exceeded (default: 0.1)
#' @return numeric scalar, x-coordinate of the valley, or NA_real_ if none found
#' @noRd
.find_signal_valley <- function(dens, peak_fraction_drop = 0.05,
                                 cap = 0.5, fallback = 0.1) {
  d_diff <- base::diff(dens$y)
  peak_idx <- base::which.max(dens$y)
  peak_height <- dens$y[peak_idx]
  min_drop <- peak_height * peak_fraction_drop
  valley_x <- NA_real_

  if (peak_idx < base::length(d_diff)) {
    current_min <- peak_height
    for (k in (peak_idx + 1):base::length(d_diff)) {
      current_min <- base::min(current_min, dens$y[k])
      if (d_diff[k] > 0 && d_diff[k - 1] <= 0 &&
          (peak_height - current_min) >= min_drop) {
        valley_x <- dens$x[k]
        break
      }
    }
  }

  # Cap: threshold should be near noise/signal boundary
  if (!base::is.na(valley_x) && valley_x > cap) {
    valley_x <- fallback
  }

  valley_x
}
