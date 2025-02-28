# optimotu.pipeline (development version)
* Fixed implementation of `fastx_split()` and `fastx_combine()` to work with
  files which contain whitespace in the header.
* Added `fasta_deline()` which converts a fasta file with sequences split into
  multiple lines into a fasta file with each sequence on a single line.

# optimotu.pipeline 0.5.1
* `ensure_directory()` now returns the name of the input file invisibly,
so that it can be used as a wrapper around the name of an output file in a
function call.

# optimotu.pipeline 0.5.0
* `calc_taxon_thresholds()` and `calc_subtaxon_thresholds()` have been moved to
the `optimotu` package (version >= 0.9)

# optimotu.pipeline 0.4.0
* Change arguments of `calc_taxon_thresholds()` and `calc_subtaxon_thresholds()`
to accept optional arguments "conf_level" and "metric", which are used to filter
the "optima" table (which is renamed from "fmeasure_optima" because a different
optimization target may be used instead)

# optimotu.pipeline 0.3.0
* Move a large number of functions from `optimotu_targets` to `optimotu.pipeline`
* New C++ implementations of `fastq_names()`, `fastx_split()` and
`fastx_combine()`

# optimotu.pipeline 0.2.1

* Fix `fastq_qual_bins()` incorrectly detecting carriage return as quality score
of -20 on Windows.

# optimotu.pipeline 0.2.0

* Add `fastq_qual_bins()` to detect which quality scores are actually used in a
fastq file.
