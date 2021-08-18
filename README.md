# epiRomics
An R package designed to integrate and visualize various levels of epigenomic information, including but not limited to: ChIP, Histone, ATAC, and RNA sequencing. epiRomics is also designed for regulatory network analysis in order to identify enhancer and enhanceosome regions from these data. 

This package is currently in development. Please contact <ammawla@ucdavis.edu> for suggestions, feedback, or bug reporting.

[**Getting Started with EpiRomics Walkthrough**](https://github.com/Huising-Lab/epiRomics/blob/main/vignettes/Getting%20Started%20with%20EpiRomics.pdf)

[**epiRomics Reference Manual**](https://github.com/Huising-Lab/epiRomics/blob/main/doc/epiRomics_0.1.2%20Reference%20Manual.pdf)

# Package installation in R

# Use if you do not have devtools currently installed
install.packages("devtools")

devtools::install_github("hadley/devtools")

# Load devtools library
library(devtools)

# Install 

install_github(repo="Huising-Lab/epiRomics")


# Load library

library(epiRomics)

# Accessing the documentation

[Getting Started with EpiRomics Walkthrough](https://github.com/Huising-Lab/epiRomics/blob/main/vignettes/Getting%20Started%20with%20EpiRomics.pdf)

[**epiRomics Reference Manual**](https://github.com/Huising-Lab/epiRomics/blob/main/doc/epiRomics_0.1.2%20Reference%20Manual.pdf)

or

help(package = 'epiRomics', help_type = 'html')


# Vignette Data Notes

This package includes some example data to get you started, delineating human pancreatic islet enhancers between alpha and beta cells. Human pancreatic islet alpha and beta ATAC- and companion RNA- Seq data were retrieved from GEO accession GSE76268 (Ackermann, et al., 2016). ATAC samples were processed using the ENCODE-DCC ATAC sequencing pipeline, aligning to the hg38 (Harrow, et al., 2012) build of the human genome (Consortium, 2012; Davis, et al., 2018). Peak calls generated through the pipeline using MACS2 (Zhang, et al., 2008) were analyzed downstream through the BioConductor package DiffBind (Ross-Innes, et al., 2012) in order to identify differentially enriched chromatin regions between the two cell types. 

RNA samples were quality controlled using the tool fastp (Chen, et al., 2018), and aligned using STAR (Dobin, et al., 2013) to the hg38 build of the human genome. Wiggle files produced by the STAR aligner were then merged by cell type using UCSC command line tools. 

All bigwigs were merged by cell type were subsetted to chromosome 1 using UCSC command line tools (Kent, et al., 2010). 

ChIP-sequencing peak calls generated using MACS2 for human pancreatic islet transcription factors Foxa2, MafB, Nkx2.2, Nkx6.1, and Pdx1 were retrieved from the EMBL-EBI repository database E-MTAB-1919 (Pasquali, et al., 2014). All peak calls were lifted over to the hg38 genome build using the UCSC genome browser liftOver tool (Kent, et al., 2002). 

Histone-sequencing peak calls generated using MACS2 for histones H3k27ac and H3k4me1 were retrieved from GEO accession GSE16256 (Bernstein, et al., 2010), and for histone H2A.Z from the EMBL-EBI repository database E-MTAB-1919 (Pasquali, et al., 2014). All peak calls were lifted over to the hg38 genome build using the UCSC genome browser liftOver tool. 


The FANTOM5 human enhancer database (Lizio, et al., 2015) was retrieved, and all regions were lifted over to the hg38 genome build using the UCSC genome browser liftOver tool. 

Human ultra-conserved non-coding elements (UCNEs) were retrieved form the UCNE database (Dimitrieva and Bucher, 2012), and all regions were lifted over to the hg38 genome build using the UCSC genome browser liftOver tool.

The human islet regulome database was retrieved (Miguel-Escalada, et al., 2019) and all regions were lifted over to the hg38 genome build using the UCSC genome browser liftOver tool.



# Citation Notes
If you use epiRomics, please cite: 

<in process of biorXiv submission as well as peer-reviewed publication. TBD. >
  
  

