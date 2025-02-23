% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{fastq_stage_flags}
\alias{fastq_stage_flags}
\alias{fastq_stage_flag2}
\title{Flag the fate of reads from a FASTQ file through processing stages}
\usage{
fastq_stage_flags(raw, stages)

fastq_stage_flag2(raw, stage1, stage2)
}
\arguments{
\item{raw}{(\code{character}) FASTQ file path, optionally gzipped. The raw reads.}

\item{stages}{(\code{character}) FASTQ file paths, optionally gzipped. The
processed reads after multiple ordered stages.}

\item{stage1}{(\code{character}) FASTQ file path, optionally gzipped, the first
processed stage}

\item{stage2}{(\code{character}) FASTQ file path, optionally gzipped, the second
processed stage}
}
\value{
a \code{data.frame} with columns \code{read_idx} (\code{integer}) and \code{flags}
(\code{raw}) indicating which processing stages each read was found in.
}
\description{
Flag the fate of reads from a FASTQ file through processing stages
}
