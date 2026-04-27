# optimotu.pipeline (development version)
* Add `parse_reference_taxonomy()` to parse common taxonomy formats from
  fasta/fastq headers used by common reference databases as well as TSV files.
* Broadened IQ-TREE model parsing for EPA-ng to cover time-reversible models
  with fitted parameterization (rates, frequencies, and heterogeneity),
  with explicit failure when parsing is not possible.
* Fixed parsing of options for BayesANT taxonomic classifier.
* Added options for optimizing cluster thresholds in `pipeline_options.yaml`.
* Exported threshold optimization option accessors
  (`do_optimize_thresholds*()` and `optimize_thresholds_file()`), and aligned
  parsing so optimization is enabled when a thresholds file is configured.
* Added new formats "qs", "qs2", "qdata" for `write_and_return_file.default()`.
  These, as well as "rds", will be autodetected from the file extension if
  possible. This might break some old usages where a file was saved in RDS
  format without an ".rds" extension (not case-sensitive).
* Fixed a bug in `fasta_rename()` and `fastq_rename()` which caused them to
  only rename the first sequence in the file. This affected pipeline outputs
  for `otu_(plausible|reliable).fasta.gz`.
* Switch from command-line `FastqIndEx` to R-native `fastqindexr` for extracting
  sequence subsets.  This should increase speed, reduce memory overhead, and
  prevent failures due to too many open file handles. The old wrappers
  `fastq_gz_extract()` and `fastq_gz_random_access_extract()` have been modified
  to use the new interface internally, so existing workflows should not break,
  but these are likely to be deprecated in the future.
* Improved robustness against 0-length and gzipped input files
* New helper function `read_long_sequence_table()` reads an OptimOTU-style long
  sequence table from a tsv file.
# optimotu.pipeline 0.6.2
* Add option `force_denovo` in `pipeline_options.yaml` to force de novo
  clustering at some or all taxonomic ranks, i.e., to ignore taxonomic
  constraints at those ranks.
* LULU implementation now works with multiple input files, allowing use in
  OptimOTU pipeline with model-aligned amplicons.
* LULU implementation now correctly maps its internal OTU indices to
  provided names/indices when the input is not already sorted.

# optimotu.pipeline 0.6.1
* Improve handling of small datasets where some objects may end up empty.
* (Re)implement option for merging multiple samples with the same name.
* Fix for parsing distance configuration from `pipeline_options.yaml`.

# optimotu.pipeline 0.6.0
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
* Add functions `fastq(_pair)?_sample_(fraction|number)?(_multiple)?`
  to repeatably sample reads from a fast file or pair of fastq files. The
  `_fraction_` variants take a numerator and denominator; the `_number_`
  variants take a target number of reads. The `_multiple` variants take multiple
  values for `numerator` or `number`, and produce multiple output files; for
  these the larger subsamples are guaranteed to include the same reads as the
  smaller subsamples. The version with neither `fraction` nor `number` takes
  an externally supplied shuffle of the sequences and a target number of reads.
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
* Add functionality pipeline rarefaction.
* Taxonomic helper functions now have a `tax_ranks` argument for use in
  an environment other than the one where pipeline options were parsed, e.g. on
  crew workers.
* `cutadapt_filter_trim()` and `cutadapt_paired_filter_trim()` now correctly
  handle multiple values in the `cut`, `cut_R1`, or `cut_R2` options, to
  unconditionally trim both ends of the read.
* Add functions `(full/large/small)_preclosed_taxon_table()` to prepare the
  pre-closed-reference-clustering table, and `do_closed_ref_cluster()` as a
  wrapper around `optimotu::closed_ref_cluster()`. The large/small variants
  select only those taxa which will be either a large or small clustering job,
  so that these can be treated differently with respect to batching and
  parallelization.
* Add functions `lulu_map()` and `lulu_table()` to run the LULU curation
  algorithm on a long OTU table.
* Add functions to parse LULU settings from pipeline options, and
 `lulu_distmx()` to calculate the distance matrix required by LULU in the
 context of the pipeline.

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
