% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rank_gse.R
\name{rank_gse}
\alias{rank_gse}
\title{This function will perform fast preranked gene set encrichment analysis}
\usage{
rank_gse(
  directory,
  rank_data,
  out_dir = "3_gse",
  from_gene,
  to_gene,
  cluster_name,
  rank_by = "logFC"
)
}
\arguments{
\item{directory}{Directory to get gmt file}

\item{rank_data}{This can be a Directory to get RDS file or directly passing in dataframe}

\item{out_dir}{Output directory}

\item{from_gene}{"ENSG" or "ENSMUSG" or "HGNC" or "MGI"}

\item{to_gene}{"MGI" or HGNC}

\item{cluster_name}{Cluster name that user want to do the geneset analysis on}

\item{rank_by}{User defined data ranking by logFC, Pvalue, or sign_Pvalue (Default setting at logFC)}
}
\value{
Outputs dataframe from fgsea function
}
\description{
This function will perform fast preranked gene set encrichment analysis
}
