% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s1_normalize_raw_counts_function.r
\name{s1_normalize_raw_counts}
\alias{s1_normalize_raw_counts}
\title{s1 normalize raw counts function following Gale lab protocol}
\usage{
s1_normalize_raw_counts(
  countfile,
  targetfile,
  gene_conversion_file = FALSE,
  target_column = 10,
  batch_column = FALSE,
  blocking_column = FALSE,
  visualize_data = TRUE,
  filter_genes_below_counts = 0,
  filter_method = "sum",
  norm_method = "none",
  results_folder = "s1_norm_raw_counts_results",
  figres = 150
)
}
\arguments{
\item{countfile}{raw counts table (generally output by htseq).}

\item{targetfile}{target file.}

\item{target_column}{columns from the target file to build design matrix for future DE analysis}

\item{blocking_column}{column to account sampling from the same animal multiple times}

\item{filter_genes_below_counts}{filter out genes with counts below a certain value, (default set to 0).}

\item{filter_method}{method by which to filter genes (sum or mean, sum is defualt)}

\item{norm_method}{normalization method to use when normalizing LogCPM}

\item{results_folder}{folder to store results in defualt is s1_norm_raw_counts_results}

\item{figres}{resolution at which to output figures (default is 300).}

\item{vizualize_data}{whether or not to generate figures (default set to true).}
}
\description{
This function allows the normalization of raw counts following a general protocol developed by previous members of the gale lab
}
\examples{
normalize_raw_counts(countfile="./p1_modified_count_matrix_results/count_file.txt", targetfile="./p1_modified_count_matrix_results/target_file.csv", gene_conversion_file="rhesus2human.csv", target_column=10, blocking_column=2, vizualize_data=TRUE, filter_genes_below_counts=0, figres=100)
}
\keyword{expression}
\keyword{gene}
\keyword{normalization}
