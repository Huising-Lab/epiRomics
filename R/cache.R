#' Download or retrieve the epiRomics data archive from BiocFileCache
#'
#' Queries BiocFileCache for an existing cached archive. If not found,
#' downloads it. If force_update is TRUE, removes the old entry and
#' re-downloads. Returns the local path to the archive file.
#'
#' @param cache A BiocFileCache object.
#' @param data_url Character URL of the data archive.
#' @param force_update Logical; if TRUE, re-download even if cached.
#' @return Character string path to the local archive file.
#' @noRd
.download_or_retrieve_archive <- function(cache, data_url, force_update) {
    query_result <- BiocFileCache::bfcquery(
        cache, "epiromics_cache_data", "rname", exact = TRUE
    )

    if (base::nrow(query_result) == 0L) {
        ## First time: download and add to cache
        base::message(
            "Downloading epiRomics example data (~1.3 GB).\n",
            "This only needs to happen once..."
        )
        archive_path <- BiocFileCache::bfcadd(
            cache,
            rname   = "epiromics_cache_data",
            fpath   = data_url,
            action  = "copy"
        )
    } else if (force_update) {
        ## Force re-download
        base::message("Re-downloading epiRomics example data...")
        rid <- query_result$rid[1L]
        BiocFileCache::bfcremove(cache, rids = rid)
        archive_path <- BiocFileCache::bfcadd(
            cache,
            rname   = "epiromics_cache_data",
            fpath   = data_url,
            action  = "copy"
        )
    } else {
        ## Already cached
        archive_path <- BiocFileCache::bfcrpath(
            cache, rids = query_result$rid[1L]
        )
    }

    base::return(archive_path)
}


#' Download and Cache epiRomics Example Data
#'
#' Downloads the epiRomics example dataset (histone marks, ChIP-seq peaks,
#' BED annotations, BigWig signal files, and differential analysis results)
#' from a remote archive and caches it locally using
#' \pkg{BiocFileCache}. Subsequent calls return the cached path without
#' re-downloading.
#'
#' @param force_update Logical; if \code{TRUE}, re-download the archive
#'   even when a cached copy exists.
#' @param ask Logical; passed to \code{BiocFileCache::BiocFileCache()}.
#'   Set to \code{FALSE} (default) for non-interactive use.
#'
#' @return A character string giving the path to the local directory
#'   containing the example data files (BigWigs, BED annotations,
#'   ChIP peaks, histones, CSV files, and template sheets).
#'
#' @details
#' The example dataset is approximately 1.3 GB compressed and includes:
#' \describe{
#'   \item{Histone/}{H3K27ac, H3K4me1, H3K27me3, H3K9me3, H3K4me3,
#'     H3K36me3, and H2A.Z peak calls (BED format)}
#'   \item{ChIP/}{Transcription factor peak calls for FOXA2, MAFB,
#'     NKX2.2, NKX6.1, and PDX1 (BED format)}
#'   \item{BED_Annotation/}{FANTOM5 enhancers, UCNEs, and Human Islet
#'     Regulome active/super enhancers (BED format)}
#'   \item{BigWigs/}{ATAC-seq and RNA-seq signal tracks for human
#'     pancreatic islet alpha and beta cells (bigWig format)}
#'   \item{CSV files}{DiffBind differential accessibility results and
#'     RNA-seq differential expression results}
#' }
#'
#' Data are downloaded once and stored in the BiocFileCache directory
#' (typically \code{~/.cache/R/BiocFileCache} on Linux/macOS).
#'
#' @section Internet requirement:
#' The first call requires internet access to download the data archive.
#' Subsequent calls work offline using the local cache. Use
#' \code{\link{has_cache}} to test whether the data is already
#' available before attempting to build vignettes or run examples.
#'
#' @seealso \code{\link{has_cache}} to check data availability
#'   without triggering a download.
#'
#' @export
#'
#' @examples
#' ## Check whether cached data is already available (no network).
#' has_cache()
#'
#' ## A lightweight toy subset is always bundled with the package and
#' ## does not require cache_data(). Use it for quick demos:
#' toy_dir <- system.file("extdata", "toy", package = "epiRomics")
#' list.files(toy_dir)
#'
#' \donttest{
#' ## Network download (~1.3 GB). Wrapped in \donttest{} per Bioconductor
#' ## guidance. Only run interactively:
#' if (interactive()) {
#'   cache_dir <- cache_data()
#' }
#' }
cache_data <- function(force_update = FALSE, ask = FALSE) {

    ## ---- check BiocFileCache availability ----
    if (!base::requireNamespace("BiocFileCache", quietly = TRUE)) {
        base::stop(
            "Package 'BiocFileCache' is required to download example data.\n",
            "Install with: BiocManager::install('BiocFileCache')",
            call. = FALSE
        )
    }

    ## ---- data URL (Zenodo stable archive) ----
    data_url <- paste0(
        "https://zenodo.org/records/19189987/files/",
        "epiromics_cache_data.tar.gz?download=1"
    )

    ## ---- set up BiocFileCache ----
    cache <- BiocFileCache::BiocFileCache(ask = ask)

    ## ---- download or retrieve cached archive via helper ----
    archive_path <- .download_or_retrieve_archive(cache, data_url, force_update)

    ## ---- extract archive ----
    cache_base   <- base::dirname(archive_path)
    extract_dir  <- base::file.path(cache_base, "epiRomics_extdata")
    sentinel     <- base::file.path(extract_dir, ".epiRomics_extracted")

    if (!base::file.exists(sentinel) || force_update) {
        base::message("Extracting epiRomics data archive...")
        if (base::dir.exists(extract_dir)) {
            base::unlink(extract_dir, recursive = TRUE)
        }
        base::dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
        utils::untar(archive_path, exdir = extract_dir)
        base::writeLines(base::as.character(base::Sys.time()), sentinel)
        base::message(
            "epiRomics data cached successfully at:\n  ",
            extract_dir
        )
    }

    ## ---- locate data root ----
    .find_data_root(extract_dir)
}


#' Check Whether epiRomics Example Data is Cached
#'
#' Tests whether the example data archive has been previously downloaded
#' and extracted using \code{\link{cache_data}}. This function
#' never triggers a download.
#'
#' @return Logical; \code{TRUE} if the data is cached and extracted,
#'   \code{FALSE} otherwise.
#'
#' @seealso \code{\link{cache_data}} to download the data.
#'
#' @export
#'
#' @examples
#' ## Check data availability (does NOT download)
#' if (has_cache()) {
#'   message("Example data is available locally.")
#' } else {
#'   message("Run cache_data() to download example data.")
#' }
has_cache <- function() {
    base::tryCatch({
        if (!base::requireNamespace("BiocFileCache", quietly = TRUE)) {
            return(FALSE)
        }
        cache <- BiocFileCache::BiocFileCache(ask = FALSE)
        query_result <- BiocFileCache::bfcquery(
            cache, "epiromics_cache_data", "rname", exact = TRUE
        )
        if (base::nrow(query_result) == 0L) return(FALSE)

        archive_path <- BiocFileCache::bfcrpath(
            cache, rids = query_result$rid[1L]
        )
        extract_dir <- base::file.path(
            base::dirname(archive_path), "epiRomics_extdata"
        )
        sentinel <- base::file.path(
            extract_dir, ".epiRomics_extracted"
        )
        base::file.exists(sentinel)
    }, error = function(e) FALSE)
}


#' Get Path to Cached epiRomics Data (No Download)
#'
#' Returns the path to previously cached example data, or \code{NULL}
#' if the data has not been downloaded yet.
#'
#' @return Character string path, or \code{NULL} if data is not cached.
#'
#' @seealso \code{\link{cache_data}},
#'   \code{\link{has_cache}}
#'
#' @export
#'
#' @examples
#' cache_path <- get_cache_path()
#' if (!is.null(cache_path)) {
#'   list.files(cache_path)
#' }
get_cache_path <- function() {
    if (!has_cache()) return(NULL)

    cache <- BiocFileCache::BiocFileCache(ask = FALSE)
    query_result <- BiocFileCache::bfcquery(
        cache, "epiromics_cache_data", "rname", exact = TRUE
    )
    archive_path <- BiocFileCache::bfcrpath(
        cache, rids = query_result$rid[1L]
    )
    extract_dir <- base::file.path(
        base::dirname(archive_path), "epiRomics_extdata"
    )

    .find_data_root(extract_dir)
}


# ---- internal helpers ------------------------------------------------

#' Locate the data root inside the extraction directory
#'
#' If the archive contains a single top-level directory, that
#' directory is the data root. Otherwise the extraction directory
#' itself is used.
#'
#' @param extract_dir Path to the extraction directory.
#' @return Normalized path to the data root.
#' @noRd
.find_data_root <- function(extract_dir) {
    contents <- base::list.files(extract_dir, full.names = TRUE)
    ## Exclude the sentinel file
    contents <- contents[
        !base::grepl(
            "^\\.epiRomics_extracted$",
            base::basename(contents)
        )
    ]

    ## Single subdirectory -> use that as root
    if (base::length(contents) == 1L &&
        base::dir.exists(contents)) {
        return(
            base::normalizePath(contents, mustWork = TRUE)
        )
    }

    base::normalizePath(extract_dir, mustWork = TRUE)
}
