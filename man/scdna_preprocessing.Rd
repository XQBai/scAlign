% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scdna_preprocessing.R
\name{scdna_preprocessing}
\alias{scdna_preprocessing}
\title{Preprocessing function for scDNA-seq data}
\usage{
scdna_preprocessing(
  genome_reference_file,
  scdna_tsv,
  barcode_txt,
  bin_size = 50,
  dims_reduce = 50,
  output_dir = "."
)
}
\arguments{
\item{genome_reference_file}{Path to genome reference bed file}

\item{scdna_tsv}{Path to scDNA-seq matrix tsv file}

\item{barcode_txt}{Path to barcode txt file}

\item{bin_size}{Number of bins to merge, default 50}

\item{dims_reduce}{Number of dimensions for reduction, default 50}

\item{output_dir}{Directory to save output files, default is current working directory}
}
\value{
list(scdna_matrix_merge,scdna_matrix_locs)
}
\description{
Preprocessing function for scDNA-seq data
}
