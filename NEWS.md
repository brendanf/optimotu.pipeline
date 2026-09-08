# optimotu.pipeline (development version)
* Add `fasta_header_seq_ids()` to extract `seq_id` from SINTAX/BOLD/UNITE/
  BayesANT FASTA headers for threshold training with
  `optimotu::optimize_thresholds(seq_names = ...)`, avoiding a duplicate
  annotated reference FASTA.
* Add `compression_level` to `cutadapt_options()` and
  `cutadapt_paired_options()`, with default value 6.
* Add `write_fastqindexr_index()` and treat `.qs2` paths like `.fqi` in
  `fastx_gz_extract()` / `fastx_gz_hash()`, seq-batch helpers, and
  `lulu_distmx()`, so plan scripts can store indexes as qs2 files and pass
  literal paths without targets tracking. Reading a large index from qs2
  takes tens of milliseconds instead of tens of seconds, which dominated
  per-batch runtime when every batch reloaded the `.fqi` index.
  `bayesant()` gains `...` for dependency-tracking arguments such as
  `hash = seqbatch_hash`.
* Promote read-pair merging to a top-level `merging:` section
  (`min_overlap`, `max_mismatch`) with method-specific defaults (dada2:
  10/1; unoise: 16/5). Remove `denoising.unoise.merge` /
  `unoise_merge_*()`; configs that still use the old location get a clear
  error. Fractional `max_mismatch` is rejected for dada2.
* Add top-level `dist_config:` inherited by `lulu` and `clustering`
  (section keys win; same-method inheritance only). Omitting
  `dist_config` in a section no longer errors; default is `usearch`.
* Add top-level `executables:` map consulted first by `find_executable()`,
  so every tool lookup can be overridden from `pipeline_options.yaml`.
* Move `added_reference` under `taxonomy.protax.added_reference`. Top-level
  `added_reference` remains accepted with a deprecation warning when
  populated; empty stubs are ignored; a populated top-level block with a
  non-Protax classifier now errors.
* Breaking: rename per-read fate map API to be denoiser-neutral.
  `seq_map()` → `dada2_read_map()`, `unoise_seq_map()` → `unoise_read_map()`,
  `merge_seq_maps()` → `merge_read_maps()`, `add_lulu_to_seq_map()` →
  `add_lulu_to_read_map()`, `add_uncross_to_seq_map()` →
  `add_uncross_to_read_map()`, `with_seqmap_annotate()` →
  `with_read_map_annotate()`. LULU's pre-parent column is renamed from
  `denoise_idx` to `prelulu_idx`.
* Add `make_denoise_map()` / `denoise_map_to_seqtable()` so matching against
  `seq_all` happens once per sample chunk; both community tables and read
  maps consume that shared map. `dada2_read_map()` and `unoise_read_map()`
  now take `denoise_map` instead of `seq_all`/`rc`.
  `make_mapped_sequence_table()` is a thin wrapper over these helpers.
* `seq_map()` (now named `dada2_read_map()`) and `unoise_seq_map()` (now
  named `unoise_read_map()`) are now vectorized across samples; these and
  `make_mapped_sequence_table.list()` (per-chunk lookup now in
  `make_denoise_map()`) avoid repeated work between samples, leading to
  large speedups when many samples are chunked.
* Split clustering job sizing: `min_ops` (default `1e6`) is the large/small
  parallel-efficiency cutoff, while `max_ops` (default `1e10`) packs both
  large and small taxa into execution groups. These are configured in
  `pipeline_options.yaml` as `clustering.min_parallel_ops` and
  `clustering.max_batch_ops`, with suitable defaults for the selected
  distance algorithm. This change also fixes batching for small clustering jobs.
* Fix `lulu_map_lowmem()` on crew remote workers: resolve and load upstream
  targets from the worker subpipeline instead of calling
  `targets::tar_meta()` / `tar_read()` during the pipeline.
* Add new helpers (and fix old helpers) for read-fate maps so that they work
  correctly with LULU enabled.
* Fix `do_denovo_cluster()` with USEARCH distances: call the exported
  `seq_cluster_usearch()` generic on a `DNAStringSet` so clustering does not
  look up the unexported `seq_cluster_usearch.DNAStringSet` method in the
  caller.
* Add vsearch UNOISE as an alternative ASV denoiser: `vsearch_fastq_merge_pairs()`
  (empty FASTQ inputs are skipped; `shards` is clamped to the number of files;
  gzipped output is written via vsearch stdout `-`, not a file named `--`),
  `merged_filter_options()`, `vsearch_cluster_unoise2()`, and
  `unoise_seq_map()` (now named `unoise_read_map()`),
  plus `make_mapped_sequence_table()` methods for `uc_cluster` objects.
* Add optional `denoising:` section in `pipeline_options.yaml` (`method`, `pool`,
  `unoise` sub-options) with accessors `denoising_method()`, `do_dada2()`,
  `do_unoise()`, `unoise_alpha()`, `unoise_minsize()`. Pair-merge settings now
  live in top-level `merging:` (see above). Default remains DADA2 when the
  section is omitted.
* Extend `filtering:` to cover both paired-read (DADA2 `maxEE_R1`/`maxEE_R2`)
  and merged-read (UNOISE `maxEE`, `maxEE_rate`, `maxNs`, `maxLen`, `minLen`)
  keys, with warnings when keys do not apply to the selected denoiser.
* Add `output` section support in `pipeline_options.yaml` via
  `parse_output_options()`: configurable tabular output formats (`rds`, `tsv`,
  `csv`, `xlsx`, `fst`, `feather`, `parquet`, `qs2`, `qdata`, `rdata`), plus
  `wide_table`
  (with legacy `dense_table` YAML alias). New accessors: `output_formats()`,
  `output_table_formats()`, `do_output_rdata()`. Default per-file formats are
  `rds` and `tsv`.
* Add `write_tabular_outputs()` to write a table to multiple formats from one
  stem path.
* Extend `write_and_return_file()` for tabular formats; add
  `write_and_return_file.list()` for bundled `.RData` export when `type` is
  `rdata` (named lists, symbol lists, or string name lists). Add
  `write_and_return_file.matrix()` with tibble coercion for tabular formats.
* Fix `parse_otu_table_options()` / `do_wide_otu_table()` so either
  `wide_table: yes` or legacy `dense_table: yes` alone enables wide output.
* Fix `hmmalign()` bugs when running with a single input sequence.
* Require `fastqindexr` (>= 0.1.0).
* Wrapper functions for many external tools (`hmmalign()`, `hmmsearch()`,
  `nhmmer()`, `run_protax()`, `run_protax_animal()`, `run_protax_besthit()`,
  `detect_numts()`, `sintax()`, `vsearch_uchime_ref()`, `bayesant()`) now use
  a more unified API which allows them to extract a subset of sequences from
  input files via `fastqindexr`. `fastq_gz_index_extract()`,
  `fastq_gz_random_access_extract()`, and `fastx_gz_hash()` have been
  refactored to also use `fastqindexr` internally for better performance.
* `ensure_directory()` is now vectorized over multiple inputs.
* Helper functions for clustering now support taxonomic sorting maps using
  character IDs as an alternative to numeric IDs.
* Sequence helper functions now use intermediate `BStringSet` instead of
  `DNAStringSet` where possible, in order to prevent unwanted alphabet
  conversion when the sequences contain non-DNA characters (including
  lower-case letters).
* Fix error when supplying non-default values for `min_taxa` or `min_refseq`
  for threshold optimization in `pipeline_options.yaml`.
* Reduce memory usage of internal C++ structures for `lulu_map()` (but memory
  usage is typically dominated by R input data).
* Add a `lulu_map_lowmem()`, which further reduces peak memory usage by
  `lulu_map()` in the context of a `targets` pipeline, by internally managing
  loading of the OTU table and match list from disk.

# optimotu.pipeline 0.6.3
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
* New helper function `tax_table_wide_to_long()`.
* Add option group `supplemental_asv` in `pipeline_options.yaml`, which allow
  ASVs from an external source (i.e., another study or reference data) to be
  included in the clustering and (optionally) taxonomy stages of the pipeline.

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
