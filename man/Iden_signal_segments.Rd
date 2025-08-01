% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scdna_subclones.R
\name{Iden_signal_segments}
\alias{Iden_signal_segments}
\title{Identify informative segments in scDNA-seq matrix using mixture model}
\usage{
Iden_signal_segments(scdna_matrix, output_dir = ".")
}
\arguments{
\item{scdna_matrix}{A matrix of scDNA-seq CNV values (segments x cells)}

\item{output_dir}{Directory to save output files, default is current working directory}
}
\value{
A vector of indices for informative segments
}
\description{
This function selects informative segments based on the standard deviation of CNVs across segments,
using a mixture model to distinguish signal from background.
}
