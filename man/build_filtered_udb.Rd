% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/usearch.R
\name{build_filtered_udb}
\alias{build_filtered_udb}
\title{Build a usearch database (UDB) file using USEARCH, removing blacklisted
sequences}
\usage{
build_filtered_udb(
  infile,
  outfile,
  type = c("usearch", "sintax", "ublast"),
  blacklist,
  usearch = find_usearch()
)
}
\arguments{
\item{infile}{(\code{character} string) input file name}

\item{outfile}{(\code{character} string) output file name}

\item{type}{(\code{character} string) type of database to build}

\item{blacklist}{(\code{character} vector) names  of sequences to remove from the
database}

\item{usearch}{(\code{character} string) path to the usearch executable}
}
\value{
(\code{character} string) the name of the output file
}
\description{
Build a usearch database (UDB) file using USEARCH, removing blacklisted
sequences
}
