# optimotu.pipeline (development version)
* Incorporated bug fixes from `optimotu_targets`.
* Added more convenience functions for accessing taxonomy:
  - `unknown_ranks()` returns the complement of `known_ranks()`
  - `tax_rank_vars()`, `known_rank_vars()`, `unknown_rank_vars()`,
    `subrank_vars()`, and `superrank_vars()` return the same ranks as their
    equivalent `*_ranks()` functions, but as a list of symbols instead of a
    character vector.
  - `define_taxonomy()` sets up the taxonomy options for use in a pipeline.
* Moved parsing of `pipeline_options.yaml` into package, using top-level
  function `parse_pipeline_options()`.
* Added functionality to run taxonomic classifiers SINTAX (`sintax()` function;
  uses external VSEARCH); BayesANT(`bayesant()` function; uses package BayesANT);
  and EPA-ng/Gappa (`epa_ng()` and `gappa_assign()` functions; uses external
  epa-ng and Gappa executables); as well as to parse pipeline options to
  configure these.
* Add `empirical_transition_matrix()` function to calculate the empirical
  transition matrix for a set of reads which have been mapped to known
  true sequences (e.g. from a mock community). This can be used to calibrate
  the DADA2 error model.
* Add functions `fastq(_pair)?_sample_(fraction|number)(_multiple)?`
  to repeatably sample reads from a fast file or pair of fastq files. The
  `_fraction_` variants take a numerator and denominator; the `_number_`
  variants take a target number of reads. The `_multiple` variants take multiple
  values for `numerator` or `number`, and produce multiple output files; for
  these the larger subsamples are guaranteed to include the same reads as the
  smaller subsamples.
* Add `tar_substitute()` function to aid in editing `targets` pipelines
  programmatically.
* Add `read_sample_table()` and `infer_sample_table()` functions to read in a
  custom sample table from a file, or infer it from the names of input read
  files.
* Add top-level functions `sample_table()`, which reads/infers the sample table
  the first time it is called in a session, and returns a cached version
  thereafter; and `sample_table_hash()`, which can be used to track changes in
  the sample table.
* Add `*_path()` functions which return the paths used to read/store
  various input, intermediate, and output files.
* Add `tar_merge()` to merge plans (typically products of `tar_map()`)
  element-wise by name.
* `tar_map_bind_rows()`, `tar_map_c()` and `tar_map_list()` now work when
  `tar_map()` was applied to `values` with 0 rows.
* `fastq_seq_map()` now correctly detects when read names are not _entirely_
  hexadecimal numbers.
* Add parsing functionality for basic raw-read rarefaction.
  
# optimotu.pipeline 0.5.2
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
