% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p1_modify_count_matrix_function.r
\name{p1_modify_count_matrix}
\alias{p1_modify_count_matrix}
\title{p1 modifying count matrix modifies count matrix}
\usage{
p1_modify_count_matrix(
  countfile,
  targetfile,
  samples_to_remove_count_matrix = FALSE
)
}
\arguments{
\item{countfile}{raw counts table (generally output by htseq).}

\item{targetfile}{target file.}

\item{samples_to_remove_count_matrix}{text file (no header) or r list of samples to remove (defualt is FALSE)}
}
\description{
This function allows remove samples from count matrix if they
 are not to be included in downstream gene expression analysis
}
\examples{
p1_modify_count_matrix(count_file.txt, target_file.csv, samples_to_remove_count_matrix=c(sample1,sample2))
}
\keyword{remove}
\keyword{samples}
