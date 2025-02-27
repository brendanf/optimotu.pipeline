% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{fasta_deline}
\alias{fasta_deline}
\title{Convert a multiline FASTA file into a single-line FASTA file}
\usage{
fasta_deline(infile, outfile, compress = FALSE)
}
\arguments{
\item{infile}{(\code{character(1)}) FASTA file, can be compressed}

\item{outfile}{(\code{character(1)}) name of output FASTA file}

\item{compress}{(\code{logical(1)}) compress output file}
}
\value{
\code{character(1)} output file name
}
\description{
Convert a multiline FASTA file into a single-line FASTA file
}
