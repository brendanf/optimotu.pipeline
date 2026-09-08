# optimotu.pipeline: Project Notes

This file is loaded automatically as project context. It summarizes
package-specific responsibilities, hotspots, and conventions for
`optimotu.pipeline`.

## 1) What `optimotu.pipeline` is responsible for

`optimotu.pipeline` is the workflow glue layer between project pipelines and
algorithmic backends. It is responsible for:

- parsing and validating `pipeline_options.yaml` semantics
- sample table inference/loading and path conventions
- wrappers around external tools (cutadapt, usearch, vsearch, hmmer, infernal)
- taxonomy adapter flows (Protax, SINTAX, BayesANT, EPA)
- taxon-aware clustering orchestration helpers built on `optimotu`
- `targets` metaprogramming helpers (`tar_map_*`, `tar_substitute`, `tar_merge`)

It is more domain/workflow-centric than `optimotu`.

## 2) Key module layout

Configuration and orchestration:

- `R/pipeline_options.R`
- `R/sample_table.R`
- `R/targets.R`
- `R/cluster.R`
- `R/taxonomy.R`

External tool wrappers:

- `R/cutadapt.R`, `R/usearch.R`, `R/vsearch.R`
- `R/hmmer.R`, `R/infernal.R`, `R/epa.R`, `R/protax.R`, `R/bayesant.R`
- `R/dada2_wrappers.R`, `R/dada2_map.R`, `R/dada2_chimera.R` (DADA2
  filtering/denoising/dereplication/merge options and per-read fate mapping,
  plus per-sample *de novo* chimera detection). `dada2_read_map()` takes a
  shared `denoise_map` from `make_denoise_map()` rather than matching
  `seq_all` itself. `seq_idx` is the current community-table id; after LULU
  that is the parent.
- `R/vsearch.R` — vsearch wrappers including paired-read merging
  (`vsearch_fastq_merge_pairs()`) and UNOISE clustering
  (`vsearch_cluster_unoise2()`, which pipes `--fastx_uniques` into
  `--cluster_unoise`)
- `R/unoise.R` — UNOISE per-read fate mapping (`unoise_read_map()`); also
  takes a shared `denoise_map`. Sequence tables from `uc_cluster` objects
  use `make_denoise_map()` / `denoise_map_to_seqtable()`.
- `R/denoise_map.R` — shared denoiser-local-id → `seq_idx` map
  (`make_denoise_map()`, `denoise_map_to_seqtable()`). One `seq_all` match
  per chunk feeds both community tables and read maps.
- `R/read_map.R` — denoiser-neutral read-map helpers (`seq_idx_lookup()`,
  `match_to_seq_all()`, `merge_read_maps()`, `with_read_map_annotate()`)

Denoising / read-quality post-processing:

- `R/numt.R` — detection of Nuclear Mitochondrial Paralogs (NUMTs) from
  `hmmalign()`/A2M alignment output
- `R/uncross.R` — removal of suspected tag-jump/cross-talk from sequence
  tables; `add_uncross_to_read_map()` sets fate-map bit `0x08` (survived
  UNCROSS). Join after LULU remap so `seq_idx` is the parent. Bits
  `0x10`--`0x80` are reserved for `asv_map$result`, not read flags.

Secondary clustering / OTU curation:

- `R/lulu_long.R` — LULU secondary denoising (Frøslev et al. 2017), merges
  putative artifact OTUs into parent OTUs post-clustering.
  `add_lulu_to_read_map()` rewrites per-read `seq_idx` to the parent and
  stores `prelulu_idx` (no extra flag bit). `lulu_map_lowmem()` must run
  with `retrieval = "none"`; on crew workers it reads deps from the
  target subpipeline rather than `tar_meta()` / `tar_read()`.
- `src/lulu.cpp` — native backend for LULU's pairwise comparisons
- configured via `parse_lulu_options()` in `R/pipeline_options.R`
  (`do_lulu`, `lulu_dist_type`, `lulu_max_dist`, etc.)

Reporting helpers:

- `R/krona.R` — generates KronaTools XML for taxonomic composition
- `R/target_taxa.R` — finds OTUs with any probability of containing a
  target taxon

Sequence/file helpers:

- `R/sequence_helpers.R`, `R/sequence_table.R`, `R/fastx_split_combine.R`,
  `R/fastq_index.R` (`write_fastqindexr_index()` for qs2 index files;
  `.qs2` / `.fqi` paths accepted by extract helpers), `R/seq_batch_input.R`
  (shared internals for indexed / sequential sequence batching used by
  `hmmalign()` and future tool wrappers), `R/write_and_return_file.R`

Small utilities:

- `R/sample_key.R` — derive a sample key from a file name
- `R/external.R` — locate external executables on the system
- `R/distribute_tasks.R` — greedy-algorithm task binning
- `R/parse_taxonomy.R` — parse taxonomy from TSV/FASTA headers;
  `fasta_header_seq_ids()` extracts bare IDs from annotated headers for
  `optimotu::optimize_thresholds(seq_names = ...)`
- `R/util.R` — misc helpers (SLURM/Snakemake/local execution detection, CPU
  count)

Native helpers:

- `src/` includes fastq/fastx, LULU, and other helper routines accessed via
  Rcpp.

## 3) Public API shape

This package exports many helpers used directly in targets plans. High-impact
families include:

- pipeline option accessors/parsers (`parse_pipeline_options()`, `do_*()`)
- sample and path accessors (`sample_table()`, `*_path()`)
- taxonomy utilities and rank helpers (`tax_ranks()`, `known_ranks()`,
  `superranks()`, `subranks()`, `build_taxonomy()`)
- clustering orchestration helpers (`do_closed_ref_cluster()`,
  `do_denovo_cluster()`, preclosed/predenovo table constructors).
  `min_ops` is the large/small parallel cutoff; `max_ops` packs both
  into execution groups. YAML keys are `clustering.min_parallel_ops` and
  `clustering.max_batch_ops`; omitted values default from `dist_config`
  method and are messaged at parse time.
- external command wrappers and sequence IO helpers
- Tool paths from `find_executable()` / `find_*()` must be resolved at plan
  definition time and passed into helpers (defaults remain for interactive
  use). Target scripts should inject with `!!`.
- DADA2 wrapper/option classes and per-read fate mapping (`dada2_wrappers.R`,
  `dada2_map.R`, `dada2_chimera.R`)
- UNOISE (vsearch) merge/cluster wrappers and per-read mapping (`vsearch.R`,
  `unoise.R`), configured through `parse_denoising_options()` /
  `parse_merge_options()` / `parse_filter_options()` (`do_unoise()`,
  `do_dada2()`, `denoising_method()`, `merge_min_overlap()`,
  `merge_max_mismatch()`)
- Shared denoise maps and read maps (`make_denoise_map()`,
  `denoise_map_to_seqtable()`, `dada2_read_map()`, `unoise_read_map()`,
  `merge_read_maps()`, `with_read_map_annotate()`)
- Top-level `dist_config:` / `executables:` / `merging:` options; Protax
  `added_reference` lives under `taxonomy.protax`
- LULU secondary-clustering entry points (`lulu_long.R`), configured through
  `parse_lulu_options()`, including `add_lulu_to_read_map()` for fate maps
  and `with_read_map_annotate()` to wrap those steps in target commands

Changes to exported option helpers are high-risk because downstream
`optimotu_targets` scripts often quote/unquote these calls inside target
commands.

## 4) Typical edit zones

- change pipeline option behavior/validation:
  - `R/pipeline_options.R`
- change sample detection and file naming:
  - `R/sample_table.R`
- change taxonomic rank or classifier handling:
  - `R/taxonomy.R`, `R/protax.R`, `R/sintax`-related helpers
- change cluster orchestration around `optimotu`:
  - `R/cluster.R`
- change external tool invocation:
  - wrapper files in `R/` listed above
- change DADA2 denoising/chimera behavior:
  - `R/dada2_wrappers.R`, `R/dada2_map.R`, `R/dada2_chimera.R`
- change UNOISE denoising or merged-read filtering:
  - `R/vsearch.R`, `R/unoise.R`, `parse_denoising_options()` /
    `parse_filter_options()` in `R/pipeline_options.R`
- change LULU secondary clustering behavior or its options:
  - `R/lulu_long.R`, `src/lulu.cpp`, `parse_lulu_options()` in
    `R/pipeline_options.R`
- change NUMT or tag-jump/cross-talk filtering:
  - `R/numt.R`, `R/uncross.R`
- change per-read fate maps (FASTA/`rc` matching, LULU remap, UNCROSS bits):
  - `R/denoise_map.R`, `R/read_map.R`, `R/dada2_map.R`, `R/unoise.R`,
    `add_lulu_to_read_map()`, `add_uncross_to_read_map()`,
    `with_read_map_annotate()`

## 5) Testing workflow

- Tests live in `tests/testthat/`.
- Validate parser and wrapper behavior when changing options or tool interfaces.
- LULU is tested extensively in `tests/testthat/test_lulu.R`.
- Run R and external tools (vsearch, cutadapt, etc.) in the development
  Apptainer image
  `../optimotu_targets/apptainer/OptimOTU_v7_dev.sif`, which has the current
  working-tree `optimotu.pipeline` plus those tools. Example:
  `apptainer exec ../optimotu_targets/apptainer/OptimOTU_v7_dev.sif R ...`
  Set `TMPDIR` and `APPTAINER_TMPDIR` to an exec-capable directory such as
  `$HOME/tmp` (host `/tmp` is `noexec`).
- If the development image is not available or does not contain required
  development package versions, prompt the user to rebuild it.
- For behavior consumed by `optimotu_targets`, also verify integration in that
  repo when possible.

## 6) Documentation and style conventions

- Use roxygen2 markdown docs.
- Do not hand-edit `NAMESPACE` or `.Rd`; regenerate with roxygen2.
- Keep comments concise and intent-focused, especially around option parsing
  and metaprogramming-heavy target factory helpers.
- Preserve `targets` dependency semantics (`!!`, `!!!`, quoted calls) when
  editing helpers used in target commands.
- UNOISE uses merge-then-denoise (vsearch `--fastq_mergepairs`, then
  `--cluster_unoise`). Do not reorder those steps to match DADA2's
  denoise-then-merge workflow.

## 7) Relationship to sister repos

- `optimotu` provides lower-level algorithmic engines used here.
- `optimotu_targets` consumes this package for concrete project targets plans.
- API/semantic changes here can ripple into both repos quickly.

## 8) Keep this file current

Before doing substantial work in this repository — option parser changes,
exported helper/API changes, taxonomy/cluster orchestration changes, or
external tool wrapper changes — review this file.

If your changes alter architecture, API semantics, or conventions documented
here, update this file in the same task. Common triggers:

- option key/semantics changes in parser/accessor functions
- exported helper additions/removals/renames
- changes to external tool assumptions
- taxonomy rank/classifier behavior changes
- shifts in boundaries between `optimotu.pipeline` and sister repos

## 9) Suggested reading order

1. `DESCRIPTION`
2. `R/pipeline_options.R`
3. `R/sample_table.R`
4. `R/targets.R`
5. `R/taxonomy.R`
6. `R/cluster.R`
7. wrapper/helper files relevant to your task
8. matching `tests/testthat/` files
