% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dada2_map.R
\name{seq_map}
\alias{seq_map}
\title{Map the fate of individual reads through merging to find unique reads}
\usage{
seq_map(
  sample,
  fq_raw,
  fq_trim,
  fq_filt,
  dadaF,
  derepF,
  dadaR,
  derepR,
  merged,
  seq_all,
  rc = FALSE
)
}
\arguments{
\item{sample}{(\code{character}) name of the sample}

\item{fq_raw}{(\code{character}) name of the raw fastq R1 file}

\item{fq_trim}{(\code{character}) name of the trimmed fastq R1 file}

\item{fq_filt}{(\code{character}) name of the filtered fastq R1 file}

\item{dadaF}{(\code{\link[dada2:dada-class]{dada2::dada}}) denoised R1}

\item{derepF}{(\code{\link[dada2:derep-class]{dada2::derep}}) dereplicated R1}

\item{dadaR}{(\code{\link[dada2:dada-class]{dada2::dada}}) denoised R2}

\item{derepR}{(\code{\link[dada2:dada-class]{dada2::dada}}) dereplicated R2}

\item{merged}{(\code{data.frame} as returned by \code{dada2::mergePairs()}) result of
merging \code{dadaF} and \code{dadaR}}

\item{seq_all}{(\code{character}) unique ASV sequences}

\item{rc}{(\code{logical}) if \code{TRUE}, sequences in \code{merged} are reverse-complemented
relative to \code{seq_all}.}
}
\value{
\code{data.frame} with columns:
-\code{sample} (character) the sample name
-\code{raw_idx} (integer) the index of the sequence in the raw file; see \code{seq_map()}
-\code{seq_idx} (integer) the index of the sequence in seq_all
-\code{flags} (raw) bitset indicating the presence of the sequence at different stages:
0x01 = trimmed
0x02 = filtered
0x04 = denoised & merged
}
\description{
Map the fate of individual reads through merging to find unique reads
}
