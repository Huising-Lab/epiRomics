#' @importFrom GenomicRanges GRanges reduce countOverlaps findOverlaps
#'   trim resize mcols intersect values setdiff
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom BiocGenerics start end width strand
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom GenomeInfoDb seqnames genome
#' @importFrom GenomicFeatures genes exonsBy
#' @importFrom rtracklayer import BigWigSelection BigWigFile
#' @importFrom ChIPseeker annotatePeak
#' @importFrom annotatr builtin_annotations read_annotations
#'   build_annotations annotatr_cache
#' @importFrom AnnotationDbi mapIds select
#'
#' @importFrom methods is new setValidity setMethod setClass setGeneric
#' @importFrom data.table %like% fcase data.table fread
#' @importFrom digest digest
#'
# Note: parallel::mclapply() / parallel::detectCores() are called in R/,
# but `parallel` lives in Suggests (it is a recommended base-R package,
# not required at load time). All parallel:: call sites are guarded by
# base::requireNamespace("parallel", quietly = TRUE), so we intentionally
# do NOT declare @importFrom parallel here — the fully-qualified calls
# resolve at runtime only when the package is present.
#
#' @importFrom graphics text segments par rect plot polygon mtext abline
#'   hist arrows points lines layout axis
#' @importFrom grDevices adjustcolor tiff postscript png pdf dev.off
#' @importFrom stats setNames density sd quantile p.adjust hclust
#'   fisher.test as.dist
#' @importFrom utils untar head read.csv packageVersion
#'   globalVariables
#' @importFrom tools file_ext
"_PACKAGE"

# Annotation packages (org.Hs.eg.db, org.Mm.eg.db, TxDb.Hsapiens.UCSC.hg38.knownGene,
# TxDb.Mmusculus.UCSC.mm10.knownGene) live in DESCRIPTION: Suggests.
# They are never called as pkg::symbol in runtime code; they appear only
# as string literals in S4 slot defaults and are resolved by the caller's
# own requireNamespace check. No @importFrom tag is declared for them.
NULL

# The following block is used by usethis to automatically
# manage roxygen namespace tags. Modify with
# care!  usethis namespace: start usethis namespace: end
NULL
