% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scAlign_integration.R
\name{assign_scrna_subclones}
\alias{assign_scrna_subclones}
\title{Assign scRNA cells to scDNA subclones using distance or correlation}
\usage{
assign_scrna_subclones(
  scdna_gene_subclone,
  scrna_proj_matrix,
  method = "euclidean"
)
}
\arguments{
\item{scdna_gene_subclone}{scDNA subclone matrix (subclones x genes)}

\item{scrna_proj_matrix}{Projected scRNA matrix in scDNA space (cells x genes)}

\item{method}{Distance method: "euclidean" or "correlation" (default "euclidean")}
}
\value{
Named vector of subclone assignments for each scRNA cell
}
\description{
Assign scRNA cells to scDNA subclones using distance or correlation
}
