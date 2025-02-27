% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{unnest_yaml_list}
\alias{unnest_yaml_list}
\title{Flatten a nested list of length-1 lists into a single list}
\usage{
unnest_yaml_list(x)
}
\arguments{
\item{x}{(\code{list}) the list to flatten}
}
\value{
a \code{list} with all nested lists unnested
}
\description{
This function is useful for unnesting lists originally imported from YAML files
formatted like:

\if{html}{\out{<div class="sourceCode">}}\preformatted{list:
  - key1: value1
  - key2: value2
}\if{html}{\out{</div>}}

which are imported as \code{list(list(key1 = value1), list(key2 = value2))}.
}
