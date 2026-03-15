# epiRomics: a multi-omics R package to identify and visualize enhancers

[![R-CMD-check](https://github.com/Huising-Lab/epiRomics/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Huising-Lab/epiRomics/actions/workflows/R-CMD-check.yaml) [![Codecov](https://codecov.io/gh/Huising-Lab/epiRomics/branch/main/graph/badge.svg)](https://codecov.io/gh/Huising-Lab/epiRomics) [![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic--2.0-blue.svg)](https://opensource.org/licenses/Artistic-2.0) [![BiocCheck](https://github.com/Huising-Lab/epiRomics/actions/workflows/BiocCheck.yaml/badge.svg)](https://github.com/Huising-Lab/epiRomics/actions/workflows/BiocCheck.yaml) [![Altmetric](https://img.shields.io/badge/Altmetric-10.1101%2F2021.08.19.456732-blue.svg)](https://www.altmetric.com/details/112068712)
  
An R package designed to integrate and visualize various levels of epigenomic information, including but not limited to: ChIP, Histone, ATAC, and RNA sequencing. epiRomics is also designed for regulatory network analysis in order to identify enhancer and enhanceosome regions from these data. 

Please contact <ammawla@ucdavis.edu> for suggestions, feedback, or bug reporting.

**Paper:** [Chromatin accessibility differences between alpha, beta, and delta cells identifies common and cell type-specific enhancers](https://link.springer.com/article/10.1186/s12864-023-09293-6) — BMC Genomics (2023)

**Package:** [epiRomics: a multi-omics R package to identify and visualize enhancers](https://www.biorxiv.org/content/10.1101/2021.08.19.456732v2) — bioRxiv (2021)

## Shiny epiRomics

Explore epiRomics results interactively through our companion Shiny web application:

**[Shiny epiRomics](https://huisinglab.com/epiRomics_2021/index.html)** — Interactive browser for pancreatic islet enhancer and enhanceosome analysis results.

## Features

- **Multi-omics Integration**: Seamlessly integrate ChIP-seq, ATAC-seq, RNA-seq, and histone modification data
- **Enhancer Identification**: Identify putative enhancer regions using histone mark combinations
- **Enhanceosome Analysis**: Detect enhanceosome regions through co-transcription factor analysis
- **Memory-Efficient Processing**: Built-in BigWig caching for large-scale data analysis
- **Comprehensive Visualization**: Generate publication-ready genomic tracks and plots
- **Robust Error Handling**: Comprehensive parameter validation and error messages

## Installation

### Prerequisites

epiRomics requires R version 4.4.0 or higher and several Bioconductor packages. Make sure you have the latest version of R and Bioconductor installed.

### Install from GitHub

```r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install epiRomics
devtools::install_github("Huising-Lab/epiRomics")

# Load the package
library(epiRomics)

# Download example data (~1.3 GB, one-time)
epiRomics_cache_data()
```

### Install from Bioconductor

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("epiRomics")

# Download example data (~1.3 GB, one-time)
epiRomics::epiRomics_cache_data()
```

## Quick Start

```r
# Load the package
library(epiRomics)

# Build epiRomics database
epiRomics_dB <- epiRomics_build_dB(
  epiRomics_db_file = "path/to/your/data.csv",
  txdb_organism = "TxDb.Hsapiens.UCSC.hg38.knownGene",
  epiRomics_genome = "hg38",
  epiRomics_organism = "org.Hs.eg.db"
)

# Identify putative enhancers
enhancers <- epiRomics_enhancers_co_marks(
  epiRomics_dB = epiRomics_dB,
  epiRomics_histone_mark_1 = "h3k4me1",
  epiRomics_histone_mark_2 = "h3k27ac"
)

# Identify enhanceosomes
enhanceosomes <- epiRomics_enhanceosome(
  epiRomics_putative_enhancers = enhancers,
  epiRomics_dB = epiRomics_dB
)

# Visualize results
epiRomics_track_layer(
  epiRomics_putative_enhanceosome = enhanceosomes,
  epiRomics_index = 1,
  epiRomics_dB = epiRomics_dB,
  epiRomics_track_connection = your_track_data
)
```

## Documentation

- **Vignette (HTML)**: [Getting Started with epiRomics](https://github.com/Huising-Lab/epiRomics/blob/main/doc/getting-started-with-epiRomics.html)
- **Vignette (PDF)**: [Getting Started with epiRomics](https://github.com/Huising-Lab/epiRomics/blob/main/doc/getting-started-with-epiRomics.pdf)
- **Reference Manual**: [Package Documentation](https://github.com/Huising-Lab/epiRomics/blob/main/doc/epiRomics-manual.pdf)
- **pkgdown site**: [https://Huising-Lab.github.io/epiRomics/](https://Huising-Lab.github.io/epiRomics/)
- **Help**: `help(package = 'epiRomics', help_type = 'html')`

## Example Data

Example data is downloaded on demand via `epiRomics_cache_data()` (~1.3 GB, cached locally after first download). The dataset focuses on delineating putative human pancreatic islet enhancers between alpha and beta cells and includes:

- Human pancreatic islet alpha and beta ATAC-seq data from GEO accession [GSE76268](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76268)
- ChIP-seq data for transcription factors Foxa2, MafB, Nkx2.2, Nkx6.1, and Pdx1
- Histone modification data (H3K27ac, H3K4me1, H2A.Z)
- FANTOM5 human enhancer database
- Human ultra-conserved non-coding elements (UCNEs)

## Memory Optimization

For large-scale analyses, epiRomics includes memory-efficient BigWig caching:

```r
# Use cached BigWig processing for memory efficiency
result <- maxCovFilesCached(
  bw_paths = your_bigwig_files,
  gr = your_genomic_regions,
  use_cache = TRUE,
  parallel = TRUE  # Enable parallel processing
)
```

## Citation

If you use epiRomics, please cite:

Alex M. Mawla, Talitha van der Meulen, Mark O. Huising. (2023). Chromatin accessibility differences between alpha, beta, and delta cells identifies common and cell type-specific enhancers. *BMC Genomics*. [10.1186/s12864-023-09293-6](https://link.springer.com/article/10.1186/s12864-023-09293-6)

Alex M. Mawla & Mark O. Huising. (2021). epiRomics: a multi-omics R package to identify and visualize enhancers. *bioRxiv* 2021.08.19.456732. [10.1101/2021.08.19.456732](https://www.biorxiv.org/content/10.1101/2021.08.19.456732v2)

## Contributing

We welcome contributions! Please feel free to submit issues, feature requests, or pull requests on our [GitHub repository](https://github.com/Huising-Lab/epiRomics).

## License

This package is licensed under the [Artistic License 2.0](https://opensource.org/licenses/Artistic-2.0). See the [LICENSE](LICENSE) file for details.

## Contact

For questions, suggestions, or bug reports, please contact:
- **Maintainer**: Alex M. Mawla <ammawla@ucdavis.edu>
- **GitHub Issues**: [https://github.com/Huising-Lab/epiRomics/issues](https://github.com/Huising-Lab/epiRomics/issues)


**Resources**:

[Getting Started with epiRomics (PDF)](https://github.com/Huising-Lab/epiRomics/blob/main/doc/getting-started-with-epiRomics.pdf)

[epiRomics Reference Manual](https://github.com/Huising-Lab/epiRomics/blob/main/doc/epiRomics-manual.pdf)

[Methods for sample data in vignette](https://www.biorxiv.org/content/biorxiv/early/2021/12/04/2021.08.19.456732/DC1/embed/media-1.docx?download=true)

[Accession table for sample data](https://www.biorxiv.org/content/biorxiv/early/2021/12/04/2021.08.19.456732/DC2/embed/media-2.xlsx?download=true)

# Vignette Data Notes

This package includes some example data to get you started, focusing on delineating putative human pancreatic islet enhancers between alpha and beta cells. 

Human pancreatic islet alpha and beta ATAC- and companion RNA- Seq data were retrieved from GEO accession [GSE76268](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76268) (Ackermann, et al., 2016). ATAC samples were processed using the [ENCODE-DCC ATAC sequencing pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline), aligning to the hg38 (Harrow, et al., 2012) build of the human genome (Consortium, 2012; Davis, et al., 2018). Peak calls generated through the pipeline using [MACS2](https://github.com/macs3-project/MACS) (Zhang, et al., 2008) were analyzed downstream through the Bioconductor package [DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html) (Ross-Innes, et al., 2012) in order to identify differentially enriched chromatin regions between the two cell types. 

RNA samples were quality controlled using the tool [Fastp](https://github.com/OpenGene/fastp) (Chen, et al., 2018), and aligned using [STAR](https://github.com/alexdobin/STAR) (Dobin, et al., 2013) to the hg38 build of the human genome. Wiggle files produced by the [STAR](https://github.com/alexdobin/STAR) aligner were then merged by cell type using [UCSC command line tools](https://github.com/ENCODE-DCC/kentUtils) (Kent, et al., 2002).

All bigwigs merged by cell type were subsetted to chromosome 1 using [UCSC command line tools](https://github.com/ENCODE-DCC/kentUtils). 

ChIP-sequencing peak calls generated using [MACS2](https://github.com/macs3-project/MACS) for human pancreatic islet transcription factors Foxa2, MafB, Nkx2.2, Nkx6.1, and Pdx1 were retrieved from the EMBL-EBI repository database [E-MTAB-1919](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1919/) (Pasquali, et al., 2014). All peak calls were lifted over to the hg38 genome build using the [UCSC genome browser liftOver tool](http://genome.ucsc.edu/cgi-bin/hgLiftOver) (Kent, et al., 2002). 

Histone-sequencing peak calls generated using [MACS2](https://github.com/macs3-project/MACS) for histones H3k27ac and H3k4me1 were retrieved from GEO accession [GSE16256](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16256) (Bernstein, et al., 2010), and for histone H2A.Z from the EMBL-EBI repository database [E-MTAB-1919](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1919/) (Pasquali, et al., 2014). All peak calls were lifted over to the hg38 genome build using the [UCSC genome browser liftOver tool](http://genome.ucsc.edu/cgi-bin/hgLiftOver).


The [FANTOM5 human enhancer database](https://fantom.gsc.riken.jp/5/) (Lizio, et al., 2015) was retrieved, and all regions were lifted over to the hg38 genome build using the UCSC genome browser liftOver tool. 

Human ultra-conserved non-coding elements (UCNEs) were retrieved from the [UCNE database](https://epd.expasy.org/ucnebase/) (Dimitrieva and Bucher, 2012), and all regions were lifted over to the hg38 genome build using the [UCSC genome browser liftOver tool](http://genome.ucsc.edu/cgi-bin/hgLiftOver).

The [human islet regulome database](http://pasqualilab.upf.edu/app/isletregulome) was retrieved (Miguel-Escalada, et al., 2019) and all regions were lifted over to the hg38 genome build using the  [UCSC genome browser liftOver tool](http://genome.ucsc.edu/cgi-bin/hgLiftOver).



# Citation Notes
If you use epiRomics, please cite:

Alex M. Mawla, Talitha van der Meulen, Mark O. Huising. (2023). [Chromatin accessibility differences between alpha, beta, and delta cells identifies common and cell type-specific enhancers.](https://link.springer.com/article/10.1186/s12864-023-09293-6) BMC Genomics. 10.1186/s12864-023-09293-6

Alex M. Mawla & Mark O. Huising. (2021). [epiRomics: a multi-omics R package to identify and visualize enhancers.](https://www.biorxiv.org/content/10.1101/2021.08.19.456732v2) bioRxiv 2021.08.19.456732

## Troubleshooting

### Installation Issues

**Problem**: Package dependencies fail to install
**Solution**: Make sure you have the latest version of R (≥ 4.4.0) and Bioconductor:

```r
# Update R packages
update.packages(ask = FALSE)

# Update Bioconductor
BiocManager::install(version = BiocManager::version())

# If installation fails, try installing core dependencies manually:
BiocManager::install(c("AnnotationDbi", "annotatr", "BiocGenerics",
                       "GenomeInfoDb", "GenomicFeatures", "GenomicRanges",
                       "IRanges", "org.Hs.eg.db", "rtracklayer",
                       "TxDb.Hsapiens.UCSC.hg38.knownGene", "data.table"))
```

**Problem**: Compilation errors on Linux
**Solution**: Install required system dependencies:

```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install -y libxml2-dev libcurl4-openssl-dev libssl-dev

# CentOS/RHEL  
sudo yum install libxml2-devel openssl-devel libcurl-devel
```

**Problem**: Memory issues with large BigWig files
**Solution**: Use the built-in caching system:

```r
# Enable caching for memory efficiency
result <- maxCovFilesCached(bw_files, gr, use_cache = TRUE)
```

### Graphics Issues

**Problem**: "Too many stacks to draw" error in visualization
**Solution**: Reduce the plotting region or increase device size:

```r
# Option 1: Focus on smaller genomic regions
smaller_gr <- resize(your_granges, width = 50000, fix = "center")

# Option 2: Increase device size
pdf("output.pdf", width = 12, height = 8)
your_plotting_function()
dev.off()
```

### Getting Help

- **Issues**: Report bugs at [GitHub Issues](https://github.com/Huising-Lab/epiRomics/issues)
- **Questions**: Email ammawla@ucdavis.edu
- **Documentation**: See the [package vignette](https://Huising-Lab.github.io/epiRomics/articles/getting-started-with-epiRomics.html)
