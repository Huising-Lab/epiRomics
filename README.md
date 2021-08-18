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

This package includes some example data to get you started, focusing on delineating putative human pancreatic islet enhancers between alpha and beta cells. 

Human pancreatic islet alpha and beta ATAC- and companion RNA- Seq data were retrieved from GEO accession [GSE76268](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76268) (Ackermann, et al., 2016). ATAC samples were processed using the [ENCODE-DCC ATAC sequencing pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline), aligning to the hg38 (Harrow, et al., 2012) build of the human genome (Consortium, 2012; Davis, et al., 2018). Peak calls generated through the pipeline using [MACS2](https://github.com/macs3-project/MACS) (Zhang, et al., 2008) were analyzed downstream through the BioConductor package [DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html) (Ross-Innes, et al., 2012) in order to identify differentially enriched chromatin regions between the two cell types. 

RNA samples were quality controlled using the tool [Fastp](https://github.com/OpenGene/fastp) (Chen, et al., 2018), and aligned using [STAR](https://github.com/alexdobin/STAR) (Dobin, et al., 2013) to the hg38 build of the human genome. Wiggle files produced by the STAR aligner were then merged by cell type using [UCSC command line tools](https://github.com/ENCODE-DCC/kentUtils) (Kent, et al., 2002)

All bigwigs were merged by cell type were subsetted to chromosome 1 using [UCSC command line tools](https://github.com/ENCODE-DCC/kentUtils). 

ChIP-sequencing peak calls generated using [MACS2](https://github.com/macs3-project/MACS) for human pancreatic islet transcription factors Foxa2, MafB, Nkx2.2, Nkx6.1, and Pdx1 were retrieved from the EMBL-EBI repository database [E-MTAB-1919](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1919/) (Pasquali, et al., 2014). All peak calls were lifted over to the hg38 genome build using the [UCSC genome browser liftOver tool](http://genome.ucsc.edu/cgi-bin/hgLiftOver) (Kent, et al., 2002). 

Histone-sequencing peak calls generated using [MACS2](https://github.com/macs3-project/MACS) for histones H3k27ac and H3k4me1 were retrieved from GEO accession [GSE16256](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16256) (Bernstein, et al., 2010), and for histone H2A.Z from the EMBL-EBI repository database E-MTAB-1919](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1919/) (Pasquali, et al., 2014). All peak calls were lifted over to the hg38 genome build using the [UCSC genome browser liftOver tool](http://genome.ucsc.edu/cgi-bin/hgLiftOver).


The [FANTOM5 human enhancer database](https://fantom.gsc.riken.jp/5/) (Lizio, et al., 2015) was retrieved, and all regions were lifted over to the hg38 genome build using the UCSC genome browser liftOver tool. 

Human ultra-conserved non-coding elements (UCNEs) were retrieved form the [UCNE database](https://ccg.epfl.ch/UCNEbase/) (Dimitrieva and Bucher, 2012), and all regions were lifted over to the hg38 genome build using the [UCSC genome browser liftOver tool](http://genome.ucsc.edu/cgi-bin/hgLiftOver).

The [human islet regulome database](http://pasqualilab.upf.edu/app/isletregulome) was retrieved (Miguel-Escalada, et al., 2019) and all regions were lifted over to the hg38 genome build using the  [UCSC genome browser liftOver tool](http://genome.ucsc.edu/cgi-bin/hgLiftOver).



# Citation Notes
If you use epiRomics, please cite: 

Alex M. Mawla & Mark O. Huising. (2021). epiRomics: a multi-omics R package to identify and visualize enhancers. *in process of biorXiv submission as well as peer-reviewed publication. TBD.*
  
  

