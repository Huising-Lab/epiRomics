# epiRomics dependency installer
# This script ensures all dependencies are available.

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  utils::install.packages("BiocManager")
}
BiocManager::install(update = FALSE, ask = FALSE)

# CRAN packages
to_install_cran <- c(
  "data.table",
  "digest",
  "igraph"
)
for (i in to_install_cran) {
  message(paste("Checking for", i))
  if (!requireNamespace(i, quietly = TRUE)) {
    message(paste("  Installing", i))
    utils::install.packages(i, repos = "https://cran.r-project.org")
  }
}

# Bioconductor packages
to_install_bc <- c(
  "AnnotationDbi",
  "annotatr",
  "BiocGenerics",
  "ChIPseeker",
  "GenomeInfoDb",
  "GenomicFeatures",
  "GenomicRanges",
  "Gviz",
  "IRanges",
  "rtracklayer",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "TxDb.Mmusculus.UCSC.mm10.knownGene"
)

for (i in to_install_bc) {
  message(paste("Checking for", i))
  if (!requireNamespace(i, quietly = TRUE)) {
    message(paste("  Installing", i))
    BiocManager::install(i, ask = FALSE)
  }
}
