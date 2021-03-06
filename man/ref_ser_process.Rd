% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ref_ser_process.R
\name{ref_ser_process}
\alias{ref_ser_process}
\title{Processes counts into a Seurat object prepared for alignment based on reference object.}
\usage{
ref_ser_process(ser, res = 0.8, ref_ser = NULL)
}
\arguments{
\item{ser}{Seurat object to process.}

\item{res}{Resolution for clustering}

\item{ref_ser}{A processed reference Seurat object used to obtain metadata
for scaling and cell selection.}
}
\value{
Outputs a processed Seurat outputs (PCA, UMAP, Phate)
}
\description{
Reads in a blank ser object (usually from align_sers.R) and processes
with a traditional Seurat pipeline. All scaling will be performed using
metadata from the reference Seurat object
}
