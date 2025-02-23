% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vsearch.R
\name{vsearch_usearch_global_closed_ref}
\alias{vsearch_usearch_global_closed_ref}
\title{Perform closed-reference clustering using vsearch}
\usage{
vsearch_usearch_global_closed_ref(query, ref, threshold, ...)
}
\arguments{
\item{query}{(\code{data.frame}, \code{\link[Biostrings:XStringSet-class]{Biostrings::DNAStringSet}}, \code{character} vector,
or file name) query sequences}

\item{ref}{(\code{data.frame}, \code{Biostrings::DNAStringSet}\code{\link[Biostrings:XStringSet-class]{Biostrings::DNAStringSet}}, \code{character} vector, or
file name) reference sequences}

\item{threshold}{(\code{numeric} scalar) identity threshold, in range 0.0-1.0}

\item{...}{additional arguments to pass to \code{vsearch_usearch_global()}}
}
\value{
\code{tibble::tibble} with columns \code{seq_id} and \code{cluster}, where \code{seq_id}
is the name of a sequence from \code{query}, and \code{cluster} is the closest match
to that sequence in \code{ref}
}
\description{
Perform closed-reference clustering using vsearch
}
