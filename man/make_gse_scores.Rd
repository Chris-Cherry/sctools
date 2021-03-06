% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_gse_scores.R
\name{make_gse_scores}
\alias{make_gse_scores}
\title{Makes gene set enrichment scores by cell}
\usage{
make_gse_scores(
  ser,
  directory = NULL,
  from_gene = "HGNC",
  to_gene = "MGI",
  gene_sets = NULL,
  csv_dir = NULL,
  type = "Real"
)
}
\arguments{
\item{ser}{Seurat object to process}

\item{directory}{Directory to get gmt file}

\item{from_gene}{"ENSG" or "ENSMUSG" or "HGNC" or "MGI"}

\item{to_gene}{"MGI" or "HGNC"}

\item{gene_sets}{User defined subset of gene set (Optional)}

\item{csv_dir}{Directory to get csv file (Optional)}

\item{type}{Output score in Absolute or Real values}
}
\value{
Outputs seurat object
}
\description{
This function will read in the gene sets from an existing list, gmt file, 
or csv file and score cells based on scaled gene expression. The results 
will be appended to the set_scores assay in the provided Seurat object or 
the assay will be created.
}
