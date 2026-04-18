## make-toy-data.R
##
## Provenance / reproducibility script for the toy dataset that ships
## under inst/extdata/toy/ and backs the "Getting Started with
## epiRomics" vignette.
##
## The toy dataset is a narrow subset of the full 1.3 GB epiRomics
## example archive hosted on Zenodo.
##

##
## The full archive is used verbatim by the companion vignette
## "Full Analysis with epiRomics". This script converts that archive
## into the small bundle the getting-started vignette needs.
##
## WHAT THIS SCRIPT DOES
##   1. Calls epiRomics::cache_data() to download and extract the full
##      1.3 GB Zenodo archive via BiocFileCache (runs once; subsequent
##      runs reuse the cache).
##   2. Selects:
##        - BIGWIG_WINDOWS: a single 400 kb window centred on the hg38
##          INS locus (chr11:1,900,000-2,300,000). Keeps the BigWigs
##          small — they dominate bundle size.
##        - BED_WINDOWS:   the BigWig window UNIONED with three
##          additional islet loci (GCG on chr2, PDX1 on chr13, MAFB
##          on chr20). BEDs are peak-call text files — adding extra
##          loci costs kilobytes and gives the TF co-binding analysis
##          in the vignette enough observations to be informative.
##   3. Subsets every track (BigWigs, histone BEDs, TF ChIP BEDs,
##      enhancer annotation BEDs) to its respective window and writes
##      the results to inst/extdata/toy/.
##   4. Regenerates the toy manifests (example_epiRomics_Db_sheet.csv
##      and example_epiRomics_BW_sheet.csv) and the toy README.
##
## EXECUTION POLICY
##   This script is a reproducibility aid. It is NOT executed by
##   R CMD build / check — `inst/scripts/` is shipped verbatim but
##   never sourced automatically. To regenerate the toy dataset
##   manually, run it in an interactive R session:
##
##       source("inst/scripts/make-toy-data.R")
##
##   The guard below refuses to run non-interactively so the script
##   can never be picked up by an automated pipeline (e.g. R CMD
##   check --run-donttest).
##
## ====================================================================

if (!interactive()) {
    message(
        "make-toy-data.R is a reproducibility aid. Run it interactively:\n",
        "  R> source(\"inst/scripts/make-toy-data.R\")"
    )
    return(invisible(NULL))
}

suppressPackageStartupMessages({
    library(rtracklayer)
    library(GenomicRanges)
    library(IRanges)
    library(epiRomics)
})

## ---- 1. Download + extract full Zenodo archive ---------------------
## cache_data() is idempotent — it reuses the BiocFileCache entry on
## every subsequent run. Force a re-pull by passing force_update = TRUE.
message("Resolving full Zenodo cache via epiRomics::cache_data() ...")
epiRomics::cache_data()
src <- epiRomics::get_cache_path()
if (is.null(src) || !dir.exists(src)) {
    stop("Full epiRomics cache not found after cache_data() — aborting.",
         call. = FALSE)
}
message("Full cache resolved at: ", src)

## ---- 2. Configure toy windows --------------------------------------
## BigWigs: keep tiny. Single 400 kb window around hg38 INS locus.
BIGWIG_WINDOWS <- GenomicRanges::GRanges(
    seqnames = "chr11",
    ranges   = IRanges::IRanges(start = 1900000L, end = 2300000L)
)

## BEDs: union of BigWig window + additional islet loci so the TF
## co-binding Fisher test has enough observations to be informative.
## These extra windows cost only kilobytes because BED files are sparse
## peak calls rather than continuous signal.
BED_WINDOWS <- c(
    BIGWIG_WINDOWS,
    ## GCG locus (glucagon) — alpha-cell master TF binding
    GenomicRanges::GRanges("chr2",  IRanges::IRanges(162800000L, 163200000L)),
    ## PDX1 locus — beta-cell master TF binding
    GenomicRanges::GRanges("chr13", IRanges::IRanges(27700000L, 28300000L)),
    ## MAFB locus — alpha-cell enriched TF binding
    GenomicRanges::GRanges("chr20", IRanges::IRanges(40600000L, 41000000L))
)

## ---- 3. Prepare output tree ----------------------------------------
pkg_root <- normalizePath(getwd())
toy_root <- file.path(pkg_root, "inst", "extdata", "toy")
if (!dir.exists(toy_root)) dir.create(toy_root, recursive = TRUE)
for (sub in c("Histone", "ChIP", "BED_Annotation", "BigWigs")) {
    d <- file.path(toy_root, sub)
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

## ---- 4. Subsetting helpers -----------------------------------------
subset_bed <- function(src_bed, out_bed, w) {
    gr  <- rtracklayer::import(src_bed)
    sub <- IRanges::subsetByOverlaps(gr, w)
    rtracklayer::export(sub, out_bed, format = "BED")
    invisible(length(sub))
}

subset_bw <- function(src_bw, out_bw, w) {
    sig <- rtracklayer::import(src_bw, which = w)
    GenomeInfoDb::seqlevels(sig, pruning.mode = "coarse") <-
        as.character(unique(GenomicRanges::seqnames(w)))
    rtracklayer::export.bw(sig, out_bw)
    invisible(length(sig))
}

## ---- 5. Subset BED files to the WIDE window ------------------------
tracks <- list(
    list(sub = "Histone",
         files = list.files(file.path(src, "Histone"),
                            pattern = "\\.bed$", full.names = TRUE)),
    list(sub = "ChIP",
         files = list.files(file.path(src, "ChIP"),
                            pattern = "\\.bed$", full.names = TRUE)),
    list(sub = "BED_Annotation",
         files = list.files(file.path(src, "BED_Annotation"),
                            pattern = "\\.bed$", full.names = TRUE))
)
for (t in tracks) {
    for (f in t$files) {
        out <- file.path(toy_root, t$sub, basename(f))
        n <- subset_bed(f, out, BED_WINDOWS)
        message(sprintf("  [%s] %-55s %7d features", t$sub, basename(f), n))
    }
}

## ---- 6. Subset BigWigs to the NARROW window ------------------------
bw_src <- list.files(
    file.path(src, "BigWigs"),
    pattern = "multi_chr\\.bigwig$",
    full.names = TRUE
)
for (bw in bw_src) {
    base <- basename(bw)
    cell <- if (grepl("Alpha", base)) "Alpha" else "Beta"
    kind <- if (grepl("\\.atac\\.", base)) "ATAC" else "RNA"
    out_name <- sprintf("%s_%s.toy.bw", cell, kind)
    out <- file.path(toy_root, "BigWigs", out_name)
    n <- subset_bw(bw, out, BIGWIG_WINDOWS)
    message(sprintf("  [BigWigs] %-55s %7d intervals", out_name, n))
}

## ---- 7. Manifests --------------------------------------------------
db <- data.frame(
    name   = c("h3k27ac", "h3k4me1", "h3k27me3", "h3k9me3", "h3k4me3",
               "h3k36me3", "h2az",
               "foxa2", "mafb", "nkx2_2", "nkx6_1", "pdx1",
               "fantom", "regulome_active", "regulome_super"),
    path   = c("Histone/H3k27ac_hg38.bed", "Histone/H3K4me1_hg38.bed",
               "Histone/H3K27me3_hg38.bed", "Histone/H3K9me3_hg38.bed",
               "Histone/H3K4me3_hg38.bed",  "Histone/H3K36me3_hg38.bed",
               "Histone/H2AZ_hg38.bed",
               "ChIP/FOXA2_hg38.bed", "ChIP/MAFB_hg38.bed",
               "ChIP/NKX2_2_hg38.bed", "ChIP/NKX6_1_hg38.bed",
               "ChIP/PDX1_hg38.bed",
               "BED_Annotation/Fantom_5.hg38.enhancers.bed",
               "BED_Annotation/Human_Regulome_hg38_Active_Enhancers.bed",
               "BED_Annotation/Human_Regulome_hg38_Super_Enhancers.bed"),
    genome = "hg38",
    format = "bed",
    type   = c(rep("histone", 7), rep("chip", 5), rep("functional", 3)),
    stringsAsFactors = FALSE
)
write.csv(db, file.path(toy_root, "example_epiRomics_Db_sheet.csv"),
          row.names = FALSE)

bw_manifest <- data.frame(
    path  = c("BigWigs/Alpha_ATAC.toy.bw", "BigWigs/Beta_ATAC.toy.bw",
              "BigWigs/Beta_RNA.toy.bw",   "BigWigs/Alpha_RNA.toy.bw"),
    name  = c("Alpha_ATAC", "Beta_ATAC", "Beta_RNA", "Alpha_RNA"),
    color = c("#e13c64",    "#008b00",   "#008b00",  "#e13c64"),
    type  = c("atac",       "atac",      "rna",      "rna"),
    stringsAsFactors = FALSE
)
write.csv(bw_manifest,
          file.path(toy_root, "example_epiRomics_BW_sheet.csv"),
          row.names = FALSE)

## ---- 8. README -----------------------------------------------------
writeLines(c(
    "# epiRomics toy dataset",
    "",
    "Small subset of the curated epiRomics example archive hosted on",
    "Zenodo (https://zenodo.org/records/19189987). The underlying sequencing",
    "data were not generated by the curator; they were assembled from",
    "public resources and reprocessed through a uniform pipeline. The",
    "full curation methodology is documented in the package README.",
    "",
    "## Windows",
    "- BigWigs (signal):  hg38 chr11:1,900,000-2,300,000 (INS locus)",
    "- BEDs (peak calls): INS locus UNION three additional islet loci",
    "  (GCG on chr2, PDX1 on chr13, MAFB on chr20) so the TF",
    "  co-binding demo has enough observations to be informative.",
    "",
    "## Contents",
    "- BigWigs/        (4 tracks: alpha/beta x ATAC/RNA)",
    "- Histone/        (7 histone mark peak calls)",
    "- ChIP/           (5 TF peak calls)",
    "- BED_Annotation/ (FANTOM5, Regulome active/super enhancers)",
    "- example_epiRomics_Db_sheet.csv",
    "- example_epiRomics_BW_sheet.csv",
    "",
    "## Regeneration",
    "The script inst/scripts/make-toy-data.R reproduces this bundle",
    "from the full Zenodo archive. It is not executed by R CMD check;",
    "source it in an interactive R session to rebuild the toy data."
), file.path(toy_root, "README.md"))

message("\nToy dataset generated at: ", toy_root)
message("Total size: ",
    format(sum(file.info(
        list.files(toy_root, recursive = TRUE, full.names = TRUE)
    )$size) / (1024^2), digits = 3), " MB")
