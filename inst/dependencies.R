# No Remotes ----
# Attachments ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  utils::install.packages("BiocManager")
BiocManager::install()

to_install_cran <-
  c("data.table",
    "party",
    "plyr")
for (i in to_install_cran) {
  base::packageStartupMessage(paste("looking for ", i))
  if (!requireNamespace(i)) {
    base::packageStartupMessage(paste("     installing", i))
    utils::install.packages(i)
  }
}

to_install_bc <-
  c(
    "AnnotationDbi",
    "annotatr",
    "BiocGenerics",
    "GenomicFeatures",
    "GenomicRanges",
    "Gviz",
    "IRanges",
    "rtracklayer",
    "org.Hs.eg.db",
    "TxDb.Hsapiens.UCSC.hg38.knownGene"
  )

for (i in to_install_bc) {
  base::packageStartupMessage(paste("looking for ", i))
  if (!requireNamespace(i)) {
    base::packageStartupMessage(paste("     installing", i))
    BiocManager::install(i, type = "source")
  }
}