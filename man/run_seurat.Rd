% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scdna_subclones.R
\name{run_seurat}
\alias{run_seurat}
\title{Run Seurat clustering on scDNA-seq matrix}
\usage{
run_seurat(seurat_obj, dims.reduce = 50, resolution = 0.5, output_dir = ".")
}
\arguments{
\item{seurat_obj}{A Seurat object}

\item{dims.reduce}{Number of dimensions for reduction (default 50)}

\item{resolution}{Clustering resolution parameter (default: 0.1 if <=600 cells, else 0.5)}

\item{output_dir}{Directory to save output files, default is current working directory}
}
\value{
A Seurat object with clustering results
}
\description{
This function performs feature selection, dimensionality reduction, and clustering on a scDNA-seq matrix.
}
