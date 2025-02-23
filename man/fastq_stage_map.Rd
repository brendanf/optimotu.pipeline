% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{fastq_stage_map}
\alias{fastq_stage_map}
\title{Map reads from a FASTQ file to processing stages}
\usage{
fastq_stage_map(raw, stages)
}
\arguments{
\item{raw}{(\code{character}) FASTQ file path, optionally gzipped. The raw reads.}

\item{stages}{(\code{character}) FASTQ file paths, optionally gzipped. The
processed reads after multiple ordered stages.}
}
\value{
a \code{data.frame} with columns \code{read_idx} (\code{integer}) and \code{flags}
(\code{raw}) indicating which processing stages each read was found in.
}
\description{
Map reads from a FASTQ file to processing stages
}
