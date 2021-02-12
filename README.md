# epiRomics
An R package designed to integrate and visualize various levels of epigenomic information, including but not limited to: ChIP, Histone, ATAC, and RNA sequencing. epiRomics is also designed for regulatory network analysis in order to identify enhancer and enhanceosome regions from these data. 

This package is currently in development. Please contact <ammawla@ucdavis.edu> for suggestions, feedback, or bug reporting.

# Package installation in R

# Use if you do not have devtools currently installed
install.packages("devtools")

devtools::install_github("hadley/devtools")

# Load devtools library
library(devtools)

# Install 
#Note this repository is private. Please do not share the line below or it's authentication token with anyone outside the lab.

install_github(repo="Huising-Lab/epiRomics", auth_token="5edce3483dbcbb83dc24dfbcef9f0685f1b9fa30")

# Load library

library(epiRomics)

# Accessing the documentation

help(package = 'epiRomics', help_type = 'html')

# Citation Notes
If you use epiRomics, please cite: 

<Methods paper citation>
  
For some background on the implementation of epiRomics, see:

<Alpha, Beta, Delta ATAC manuscript citation>

<ObOb/Lean Beta ATAC manuscript citation>
