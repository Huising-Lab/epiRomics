% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epiRomics_track_layer_human.R
\name{epiRomics_track_layer_human}
\alias{epiRomics_track_layer_human}
\title{Visualizes data from epiRomics results}
\usage{
epiRomics_track_layer_human(
  epiRomics_putative_enhanceosome,
  epiRomics_index,
  epiRomics_dB,
  epiRomics_track_connection,
  epiRomics_keep_epitracks = TRUE
)
}
\arguments{
\item{epiRomics_putative_enhanceosome}{epiRomics class database containing putative enhanceosome calls}

\item{epiRomics_index}{numeric of row value from epiRomics_putative_enhanceosome to visualize}

\item{epiRomics_dB}{epiRomics class database containing all data initially loaded}

\item{epiRomics_track_connection}{data frame containing bigwig track locations and their names}

\item{epiRomics_keep_epitracks}{logical indicating whether to show enhancer and chip tracks, default is TRUE}
}
\value{
GViz plot
}
\description{
Visualizes data from epiRomics results
}
