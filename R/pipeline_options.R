#' Flatten a nested list of length-1 lists into a single list
#'
#' This function is useful for unnesting lists originally imported from YAML
#' files formatted like:
#' ```
#' list:
#'   - key1: value1
#'   - key2: value2
#' ```
#' which are imported as `list(list(key1 = value1), list(key2 = value2))`.
#'
#' @param x (`list`) the list to flatten
#' @return a `list` with all nested lists unnested
#' @keywords internal
unnest_yaml_list <- function(x) {
  checkmate::assert_list(x)
  if (
    is.null(names(x)) &&
      checkmate::check_list(x, types = "list") &&
      all(vapply(x, length, 1L) == 1L)
  ) {
    do.call(c, x)
  } else {
    x
  }
}

#### project_name ####

#' Parse pipeline options from a YAML file
#'
#' These functions parse the pipeline options from a YAML file and set
#' package options accordingly.
#'
#' @rdname parse_pipeline_options
#' @param pipeline_options (`list`) the pipeline options to parse
#' @keywords internal
parse_project_name <- function(pipeline_options) {
  if (
    !("project_name" %in% names(pipeline_options)) ||
      length(pipeline_options$project_name) == 0
  ) {
    warning(
      "Missing project name in 'pipeline_options.yaml'.\n",
      "Using default name 'metabarcoding_project'."
    )
    return("metabarcoding_project")
  } else if (length(pipeline_options$project_name) > 1) {
    stop("Project can only have one name (file: pipeline_options.yaml)")
  } else if (pipeline_options$project_name == "metabarcoding_project") {
    message(
      "Option 'project_name' is the default value 'metabarcoding_project'.\n",
      "You can change it by editing the file 'pipeline_options.yaml'"
    )
  } else if (!grepl("^[[:alnum:]_-]+$", pipeline_options$project_name)) {
    stop(
      "Project name should consist of alphanumeric characters, '_', and",
      " '-'. (file:pipeline_options.yaml)"
    )
  }
  options(optimotu.pipeline.project_name = pipeline_options$project_name)
}

#' Get the currently defined project name
#' @return (`character` string) the project name
#' @export
project_name <- function() {
  getOption("optimotu.pipeline.project_name", "metabarcoding_project")
}

#### file_extension ####
#' @rdname parse_pipeline_options
#' @keywords internal
parse_file_extension <- function(pipeline_options) {
  checkmate::assert_string(
    pipeline_options$file_extension,
    null.ok = TRUE
  )
  if (!is.null(pipeline_options$file_extension)) {
    options(
      optimotu.pipeline.read_file_extension = pipeline_options$file_extension
    )
  }
}

#' Functions to access pipeline-wide options
#'
#' When used in a targets plan, these should always be
#' pre-evaluated with `!!` or `!!!` to ensure proper dependency tracking.
#'
#' @name pipeline_options
NULL

#' @rdname pipeline_options
#' @export
read_file_extension <- function() {
  getOption("optimotu.pipeline.read_file_extension", "f(?:ast)?q(?:[.]gz)?")
}

#### read orientation ####
#' @rdname parse_pipeline_options
#' @keywords internal
parse_orient <- function(pipeline_options) {
  checkmate::assert_string(
    pipeline_options$orient,
    null.ok = TRUE,
    pattern = "^(fwd|rev|mixed|custom)$"
  )
  if (is.null(pipeline_options$orient)) {
    message(
      "No read orientation specified in 'pipeline_options.yaml'.\n",
      "Using default: 'forward'"
    )
    options(optimotu.pipeline.read_orientation = "fwd")
  } else {
    options(optimotu.pipeline.read_orientation = pipeline_options$orient)
  }
}

#' Get the read orientation setting
#' @return (`character` string) the read orientation setting, one of
#' `"fwd"`, `"rev"`, `"mixed"`, or `"custom"`
#' @export
read_orientation <- function() {
  getOption("optimotu.pipeline.read_orientation", "fwd")
}

#### duplicate sample policy ####
#' @rdname parse_pipeline_options
#' @keywords internal
parse_duplicate_policy <- function(pipeline_options) {
  checkmate::assert_string(
    pipeline_options$duplicate_samples,
    null.ok = TRUE
  )
  if (!is.null(pipeline_options$duplicate_samples)) {
    checkmate::assert_choice(
      pipeline_options$duplicate_samples,
      choices = c("merge", "error")
    )
    options(
      optimotu.pipeline.duplicate_samples = pipeline_options$duplicate_samples
    )
  }
}

#' @rdname pipeline_options
#' @export
duplicate_samples <- function() {
  getOption("optimotu.pipeline.duplicate_samples", "error")
}

#### custom sample table ####
#' @rdname parse_pipeline_options
#' @keywords internal
parse_custom_sample_table <- function(pipeline_options) {
  checkmate::assert(
    checkmate::check_null(pipeline_options$custom_sample_table),
    checkmate::check_false(pipeline_options$custom_sample_table),
    checkmate::check_file_exists(pipeline_options$custom_sample_table)
  )
  if (is.character(pipeline_options$custom_sample_table)) {
    options(
      optimotu.pipeline.custom_sample_table = pipeline_options$custom_sample_table
    )
  }
}

#' @rdname pipeline_options
#' @export
custom_sample_table <- function() {
  getOption("optimotu.pipeline.custom_sample_table", NULL)
}

#' @rdname pipeline_options
do_custom_sample_table <- function() {
  is.character(custom_sample_table())
}

#### supplemental_asv ####
#' @rdname parse_pipeline_options
#' @param set_name (`character` scalar) name of the supplemental ASV set
#' @keywords internal
assert_supplemental_header_ids <- function(path, set_name) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(con), add = TRUE)
  repeat {
    chunk <- readLines(con, n = 10000L, warn = FALSE)
    if (length(chunk) == 0L) {
      break
    }
    header_idx <- startsWith(chunk, ">")
    if (!any(header_idx)) {
      next
    }
    headers <- sub("^>", "", chunk[header_idx])
    bad <- grepl("[[:space:];]", headers)
    if (any(bad)) {
      bad_header <- headers[which(bad)[1]]
      stop(
        "Invalid FASTA header in supplemental ASV set '",
        set_name,
        "'. ",
        "Sequence IDs must not contain whitespace or semicolons. ",
        "First invalid header: '",
        bad_header,
        "'."
      )
    }
  }
  invisible(NULL)
}

#' @rdname parse_pipeline_options
#' @keywords internal
parse_supplemental_asv_options <- function(pipeline_options) {
  options(
    optimotu.pipeline.do_supplemental_asv = FALSE,
    optimotu.pipeline.supplemental_asv = list()
  )
  if (is.null(pipeline_options$supplemental_asv)) {
    return(invisible(NULL))
  }
  checkmate::assert_list(pipeline_options$supplemental_asv)
  supplemental_opts <- unnest_yaml_list(pipeline_options$supplemental_asv)
  if (
    is.null(names(supplemental_opts)) || any(names(supplemental_opts) == "")
  ) {
    stop("'supplemental_asv' must be a named list of sets.")
  }
  duplicate_set_names <- unique(names(supplemental_opts)[duplicated(names(
    supplemental_opts
  ))])
  if (length(duplicate_set_names) > 0L) {
    stop(
      "Duplicate supplemental ASV set name(s): ",
      paste(sprintf("'%s'", duplicate_set_names), collapse = ", "),
      ". Set names must be unique."
    )
  }

  parsed <- lapply(names(supplemental_opts), function(set_name) {
    set_opts <- supplemental_opts[[set_name]]
    checkmate::assert_list(set_opts, names = "named")
    checkmate::assert_names(
      names(set_opts),
      subset.of = c("sequences", "sample_table", "taxonomy", "enabled")
    )
    set_opts <- unnest_yaml_list(set_opts)

    checkmate::assert_string(set_opts$sequences)
    checkmate::assert_file_exists(set_opts$sequences, access = "r")

    if (!is.null(set_opts$sample_table) && !isFALSE(set_opts$sample_table)) {
      checkmate::assert_string(set_opts$sample_table)
      checkmate::assert_file_exists(set_opts$sample_table, access = "r")
    }

    checkmate::assert(
      checkmate::check_null(set_opts$taxonomy),
      checkmate::check_false(set_opts$taxonomy),
      checkmate::check_flag(set_opts$taxonomy),
      checkmate::check_string(set_opts$taxonomy)
    )

    taxonomy_mode <- "none"
    taxonomy_file <- NULL
    if (isTRUE(set_opts$taxonomy)) {
      taxonomy_mode <- "header"
    } else if (is.character(set_opts$taxonomy)) {
      taxonomy_mode <- "file"
      taxonomy_file <- set_opts$taxonomy
      checkmate::assert_file_exists(taxonomy_file, access = "r")
    }
    if (!identical(taxonomy_mode, "header")) {
      assert_supplemental_header_ids(set_opts$sequences, set_name)
    }

    enabled <- if (is.null(set_opts$enabled)) TRUE else set_opts$enabled
    checkmate::assert_flag(enabled)

    list(
      name = set_name,
      sequences = set_opts$sequences,
      sample_table = if (is.character(set_opts$sample_table)) {
        set_opts$sample_table
      } else {
        NULL
      },
      taxonomy_mode = taxonomy_mode,
      taxonomy_file = taxonomy_file,
      enabled = enabled
    )
  })
  names(parsed) <- names(supplemental_opts)

  active <- parsed[vapply(parsed, `[[`, logical(1), "enabled")]

  options(
    optimotu.pipeline.do_supplemental_asv = length(active) > 0L,
    optimotu.pipeline.supplemental_asv = parsed
  )
}

#' @rdname pipeline_options
#' @export
do_supp_asv <- function() {
  getOption("optimotu.pipeline.do_supplemental_asv", FALSE)
}

#' @rdname pipeline_options
#' @param enabled_only (`logical` flag) if TRUE, only return enabled sets
#' @export
supp_asv_sets <- function(enabled_only = TRUE) {
  checkmate::assert_flag(enabled_only)
  sets <- getOption("optimotu.pipeline.supplemental_asv", list())
  if (isTRUE(enabled_only)) {
    sets <- sets[vapply(sets, `[[`, logical(1), "enabled")]
  }
  sets
}

#' @rdname pipeline_options
#' @export
supp_asv_set_names <- function(enabled_only = TRUE) {
  names(supp_asv_sets(enabled_only = enabled_only))
}

#' @rdname pipeline_options
#' @param set_name (`character` vector) name(s) of the supplemental ASV set(s)
#'   to access
#' @export
supp_asv_set <- function(set_name) {
  checkmate::assert_character(set_name, min.len = 1, any.missing = FALSE)
  sets <- supp_asv_sets(enabled_only = FALSE)
  unknown <- setdiff(unique(set_name), names(sets))
  if (length(unknown) > 0) {
    stop(
      "Unknown supplemental ASV set(s): ",
      paste(sprintf("'%s'", unknown), collapse = ", "),
      "."
    )
  }
  out <- sets[set_name]
  if (length(set_name) == 1L) {
    out[[1]]
  } else {
    out
  }
}

#' @rdname pipeline_options
#' @export
supp_asv_sequences <- function(set_name) {
  out <- supp_asv_set(set_name)
  if (length(set_name) == 1L) {
    out$sequences
  } else {
    vapply(out, `[[`, character(1), "sequences")
  }
}

#' @rdname pipeline_options
#' @export
supp_asv_sample_table <- function(set_name) {
  out <- supp_asv_set(set_name)
  if (length(set_name) == 1L) {
    out$sample_table
  } else {
    vapply(
      out,
      function(x) {
        if (is.null(x$sample_table)) NA_character_ else x$sample_table
      },
      character(1)
    )
  }
}

#' @rdname pipeline_options
#' @export
supp_asv_taxonomy_mode <- function(set_name) {
  out <- supp_asv_set(set_name)
  if (length(set_name) == 1L) {
    out$taxonomy_mode
  } else {
    vapply(out, `[[`, character(1), "taxonomy_mode")
  }
}

#' @rdname pipeline_options
#' @export
supp_asv_taxonomy_file <- function(set_name) {
  out <- supp_asv_set(set_name)
  if (length(set_name) == 1L) {
    out$taxonomy_file
  } else {
    vapply(
      out,
      function(x) {
        if (is.null(x$taxonomy_file)) NA_character_ else x$taxonomy_file
      },
      character(1)
    )
  }
}

#### added_reference ####
#' Apply added_reference fasta/table paths to options
#' @param added_reference (`list`) with optional `fasta` and `table` entries
#' @noRd
apply_added_reference <- function(added_reference) {
  checkmate::assert_list(added_reference)
  added_reference <- unnest_yaml_list(added_reference)
  checkmate::assert_string(added_reference$fasta, null.ok = TRUE)
  checkmate::assert_string(added_reference$table, null.ok = TRUE)

  if (xor(is.null(added_reference$fasta), is.null(added_reference$table))) {
    stop(
      "If one of 'added_reference: fasta' and 'added_reference: table' is ",
      "given in 'pipeline_options.yaml', then both must be given.",
      call. = FALSE
    )
  }
  if (!is.null(added_reference$fasta)) {
    checkmate::assert_file_exists(added_reference$fasta, access = "r")
    checkmate::assert_file_exists(added_reference$table, access = "r")
    options(
      optimotu.pipeline.do_added_reference = TRUE,
      optimotu.pipeline.added_reference_fasta = added_reference$fasta,
      optimotu.pipeline.added_reference_table = added_reference$table
    )
  }
}

#' Whether an added_reference block is populated (not just an empty stub)
#' @param added_reference (`list` or `NULL`)
#' @return (`logical`)
#' @noRd
added_reference_is_populated <- function(added_reference) {
  if (is.null(added_reference) || !is.list(added_reference)) {
    return(FALSE)
  }
  added_reference <- unnest_yaml_list(added_reference)
  !is.null(added_reference$fasta) || !is.null(added_reference$table)
}

#' @rdname parse_pipeline_options
#' @keywords internal
parse_added_reference <- function(pipeline_options) {
  if (is.null(pipeline_options$added_reference)) {
    return(invisible())
  }
  checkmate::assert_list(pipeline_options$added_reference)
  if (!added_reference_is_populated(pipeline_options$added_reference)) {
    return(invisible())
  }
  if (do_added_reference()) {
    stop(
      "added_reference is specified both at the top level and under ",
      "taxonomy.protax. Prefer taxonomy.protax.added_reference.",
      call. = FALSE
    )
  }
  if (!do_protax()) {
    stop(
      "Top-level 'added_reference' is only valid with the Protax classifier ",
      "(taxonomy.protax). Move it under taxonomy.protax.added_reference, or ",
      "remove it.",
      call. = FALSE
    )
  }
  warning(
    "Top-level 'added_reference' is deprecated; move it under ",
    "taxonomy.protax.added_reference in pipeline_options.yaml.",
    call. = FALSE
  )
  apply_added_reference(pipeline_options$added_reference)
}

#' @rdname pipeline_options
#' @export
do_added_reference <- function() {
  getOption("optimotu.pipeline.do_added_reference", FALSE)
}

#' @rdname pipeline_options
#' @export
added_reference_fasta <- function() {
  getOption("optimotu.pipeline.added_reference_fasta", NULL)
}

#' @rdname pipeline_options
#' @export
added_reference_table <- function() {
  getOption("optimotu.pipeline.added_reference_table", NULL)
}

#### executables ####
#' @rdname parse_pipeline_options
#' @keywords internal
parse_executable_options <- function(pipeline_options) {
  executables <- pipeline_options$executables
  if (is.null(executables)) {
    return(invisible())
  }
  checkmate::assert_list(executables)
  executables <- unnest_yaml_list(executables)
  if (length(executables) == 0L) {
    return(invisible())
  }
  checkmate::assert_names(names(executables), type = "unique")
  for (nm in names(executables)) {
    checkmate::assert_string(executables[[nm]], min.chars = 1)
  }
  options(
    optimotu.pipeline.executables = stats::setNames(
      as.character(unlist(executables, use.names = FALSE)),
      names(executables)
    )
  )
}

#' @rdname pipeline_options
#' @export
configured_executables <- function() {
  getOption("optimotu.pipeline.executables", character())
}

#### parallelism ####
#' @rdname parse_pipeline_options
#' @keywords internal
parse_parallel_options <- function(pipeline_options) {
  checkmate::assert_count(
    pipeline_options$local_threads,
    positive = TRUE,
    null.ok = TRUE
  )
  if (!is.null(pipeline_options$local_threads)) {
    options(optimotu_num_threads = pipeline_options$local_threads)
  }

  checkmate::assert_count(
    pipeline_options$max_batchsize,
    na.ok = TRUE,
    null.ok = TRUE
  )
  max_batchsize <- NULL
  if (checkmate::test_count(pipeline_options$max_batchsize, positive = TRUE)) {
    max_batchsize <- pipeline_options$max_batchsize
  }

  workers_per_seqrun <- 1L
  checkmate::assert_count(
    pipeline_options$workers_per_seqrun,
    positive = TRUE,
    null.ok = TRUE
  )
  checkmate::assert_count(
    pipeline_options$jobs_per_seqrun,
    positive = TRUE,
    null.ok = TRUE
  )
  if (!is.null(pipeline_options$workers_per_seqrun)) {
    if (!is.null(pipeline_options$jobs_per_seqrun)) {
      warning(
        "both 'workers_per_seqrun' and 'jobs_per_seqrun' (deprecated) were ",
        "given in 'pipeline_options.yaml'. Using 'workers_per_seqrun' (=",
        pipeline_options$workers_per_seqrun,
        ")"
      )
    }
    workers_per_seqrun <- pipeline_options$workers_per_seqrun
  } else if (!is.null(pipeline_options$jobs_per_seqrun)) {
    message(
      "Option 'jobs_per_seqrun' is deprecated in 'pipeline_options.yaml'.",
      "Please use 'workers_per_seqrun' instead."
    )
    workers_per_seqrun <- pipeline_options$jobs_per_seqrun
  }

  min_workers <- 1L
  checkmate::assert_int(
    pipeline_options$min_workers,
    lower = 1L,
    null.ok = TRUE
  )
  if (!is.null(pipeline_options$min_workers)) {
    min_workers <- pipeline_options$min_workers
  }

  max_workers <- Inf
  checkmate::assert_int(
    pipeline_options$max_workers,
    lower = min_workers,
    null.ok = TRUE
  )
  if (!is.null(pipeline_options$max_workers)) {
    max_workers <- pipeline_options$max_workers
  }

  options(
    optimotu.pipeline.max_workers = max_workers,
    optimotu.pipeline.min_workers = min_workers,
    optimotu.pipeline.workers_per_seqrun = workers_per_seqrun,
    optimotu.pipeline.max_batchsize = max_batchsize
  )
}

#' @rdname pipeline_options
#' @export
max_workers <- function() {
  getOption("optimotu.pipeline.max_workers", Inf)
}

#' @rdname pipeline_options
#' @export
min_workers <- function() {
  getOption("optimotu.pipeline.min_workers", 1L)
}

#' @rdname pipeline_options
#' @export
workers_per_seqrun <- function() {
  getOption("optimotu.pipeline.workers_per_seqrun", 1L)
}

#' @rdname pipeline_options
#' @export
max_batchsize <- function() {
  getOption("optimotu.pipeline.max_batchsize", NULL)
}

#' @rdname pipeline_options
#' @export
n_seqrun <- function() {
  dplyr::n_distinct(sample_table()$seqrun)
}

#' @rdname pipeline_options
#' @export
n_orient_seqrun <- function() {
  sample_table <- sample_table()
  dplyr::n_distinct(sample_table$seqrun, sample_table$orient)
}

#' @rdname pipeline_options
#' @export
n_workers <- function() {
  max(
    min_workers(),
    min(
      max_workers(),
      n_orient_seqrun() * workers_per_seqrun()
    )
  )
}

#### primers ####

#' @rdname parse_pipeline_options
#' @keywords internal
parse_forward_primer <- function(pipeline_options) {
  checkmate::assert_string(
    pipeline_options$forward_primer,
    null.ok = TRUE,
    min.chars = 10,
    pattern = "[ACGTSWRYMKBDHVIN]+",
    ignore.case = TRUE
  )
  if (is.null(pipeline_options$forward_primer)) {
    message(
      "Forward primer string missing (file: pipeline_options.yaml)\n",
      "Using default: GCATCGATGAAGAACGCAGC"
    )
  } else {
    options(optimotu.pipeline.forward_primer = pipeline_options$forward_primer)
  }
}

#' @rdname pipeline_options
#' @return (`character` string) the forward primer
#' @export
forward_primer <- function() {
  getOption("optimotu.pipeline.forward_primer", "GCATCGATGAAGAACGCAGC")
}

#' @rdname parse_pipeline_options
#' @export
parse_reverse_primer <- function(pipeline_options) {
  checkmate::assert_string(
    pipeline_options$reverse_primer,
    null.ok = TRUE,
    min.chars = 10,
    pattern = "[ACGTSWRYMKBDHVIN]+",
    ignore.case = TRUE
  )
  if (is.null(pipeline_options$reverse_primer)) {
    message(
      "Reverse primer string missing (file: pipeline_options.yaml)\n",
      "Using default: TCCTCCGCTTATTGATATGC"
    )
  } else {
    options(optimotu.pipeline.reverse_primer = pipeline_options$reverse_primer)
  }
}

#' @rdname pipeline_options
#' @return (`character` string) the reverse primer
#' @export
reverse_primer <- function() {
  getOption("optimotu.pipeline.reverse_primer", "TCCTCCGCTTATTGATATGC")
}

#' Get the combined primer string for trimming reads
#' @return (`character`) string with the primer setting
#' @export
#' @rdname trim_primer
trim_primer_R1 <- function() {
  sprintf("%s...%s;optional", forward_primer(), dada2::rc(reverse_primer()))
}

#' @rdname trim_primer
#' @export
trim_primer_R2 <- function() {
  sprintf("%s...%s;optional", reverse_primer(), dada2::rc(forward_primer()))
}

#' @rdname trim_primer
#' @export
trim_primer_merged <- function() {
  sprintf("%s...%s", forward_primer(), dada2::rc(reverse_primer()))
}

#### trimming settings ####
#' @rdname parse_pipeline_options
#' @export
#' @param pipeline_options (`list`) the pipeline options to parse
parse_trim_options <- function(pipeline_options) {
  checkmate::assert_list(pipeline_options$trimming, null.ok = TRUE)
  if (is.null(pipeline_options$trimming)) {
    message(
      "No 'trimming' options given in 'pipeline_options.yaml'\n",
      "Using defaults."
    )
  } else {
    options(
      optimotu.pipeline.trim_options = do.call(
        cutadapt_paired_options,
        unnest_yaml_list(pipeline_options$trimming)
      )
    )
  }
}

#' @rdname pipeline_options
#' @export
trim_options <- function() {
  getOption("optimotu.pipeline.trim_options", cutadapt_paired_options())
}

#### denoising settings ####
paired_filter_option_names <- c("maxEE_R1", "maxEE_R2")

#' @rdname parse_pipeline_options
#' @keywords internal
parse_unoise_options <- function(unoise_options) {
  if (is.null(unoise_options)) {
    unoise_options <- list()
  } else {
    unoise_options <- unnest_yaml_list(unoise_options)
  }
  checkmate::assert_list(unoise_options)
  if ("merge" %in% names(unoise_options)) {
    stop(
      "'denoising.unoise.merge' has moved to the top-level 'merging:' ",
      "section of pipeline_options.yaml. Use merging: min_overlap / ",
      "max_mismatch instead of denoising.unoise.merge: min_overlap / ",
      "max_diffs.",
      call. = FALSE
    )
  }
  if (length(unoise_options) > 0L) {
    checkmate::assert_names(
      names(unoise_options),
      subset.of = c("alpha", "minsize")
    )
  }
  if ("alpha" %in% names(unoise_options)) {
    checkmate::assert_number(unoise_options$alpha, lower = 0, finite = TRUE)
    options(optimotu.pipeline.unoise_alpha = unoise_options$alpha)
  }
  if ("minsize" %in% names(unoise_options)) {
    checkmate::assert_count(unoise_options$minsize, positive = TRUE)
    options(
      optimotu.pipeline.unoise_minsize = as.integer(unoise_options$minsize)
    )
  }
}

#' @rdname parse_pipeline_options
#' @keywords internal
parse_denoising_options <- function(pipeline_options) {
  denoising <- pipeline_options$denoising
  if (is.null(denoising)) {
    options(
      optimotu.pipeline.denoising_method = "dada2",
      optimotu.pipeline.denoising_pool = "sample"
    )
    return(invisible())
  }
  if (checkmate::test_string(denoising)) {
    denoising <- list(method = denoising)
  } else {
    checkmate::assert_list(denoising)
    denoising <- unnest_yaml_list(denoising)
  }
  checkmate::assert_names(
    names(denoising),
    subset.of = c("method", "pool", "unoise", "poanoise")
  )
  method <- denoising$method
  if (is.null(method)) {
    method <- "dada2"
  }
  checkmate::assert_choice(method, c("dada2", "unoise", "poanoise"))
  if (identical(method, "poanoise")) {
    stop(
      "denoising method 'poanoise' is reserved for a future release and is ",
      "not implemented yet.",
      call. = FALSE
    )
  }
  pool <- denoising$pool
  if (is.null(pool)) {
    pool <- "sample"
  }
  checkmate::assert_choice(pool, c("sample", "asv", "derep", "project"))
  if (!identical(pool, "sample")) {
    stop(
      "denoising pool '",
      pool,
      "' is reserved for a future release; only 'sample' is implemented.",
      call. = FALSE
    )
  }
  options(
    optimotu.pipeline.denoising_method = method,
    optimotu.pipeline.denoising_pool = pool
  )
  if (identical(method, "unoise")) {
    parse_unoise_options(denoising$unoise)
  } else if (!is.null(denoising$unoise)) {
    warning(
      "Ignoring 'denoising.unoise' options because method is '",
      method,
      "'.",
      call. = FALSE
    )
  }
}

#' @rdname pipeline_options
#' @export
denoising_method <- function() {
  getOption("optimotu.pipeline.denoising_method", "dada2")
}

#' @rdname pipeline_options
#' @export
do_dada2 <- function() {
  identical(denoising_method(), "dada2")
}

#' @rdname pipeline_options
#' @export
do_unoise <- function() {
  identical(denoising_method(), "unoise")
}

#' @rdname pipeline_options
#' @export
denoising_pool <- function() {
  getOption("optimotu.pipeline.denoising_pool", "sample")
}

#' @rdname pipeline_options
#' @export
unoise_alpha <- function() {
  getOption("optimotu.pipeline.unoise_alpha", 2)
}

#' @rdname pipeline_options
#' @export
unoise_minsize <- function() {
  getOption("optimotu.pipeline.unoise_minsize", 8L)
}

#### merging settings ####
#' @rdname parse_pipeline_options
#' @keywords internal
parse_merge_options <- function(pipeline_options) {
  merging <- pipeline_options$merging
  if (is.null(merging)) {
    return(invisible())
  }
  checkmate::assert_list(merging)
  merging <- unnest_yaml_list(merging)
  checkmate::assert_names(
    names(merging),
    subset.of = c("min_overlap", "max_mismatch")
  )
  if ("min_overlap" %in% names(merging)) {
    checkmate::assert_integerish(merging$min_overlap, lower = 5, len = 1)
    options(
      optimotu.pipeline.merge_min_overlap = as.integer(merging$min_overlap)
    )
  }
  if ("max_mismatch" %in% names(merging)) {
    checkmate::assert_number(merging$max_mismatch, lower = 0, finite = TRUE)
    if (do_dada2() && merging$max_mismatch < 1) {
      stop(
        "merging.max_mismatch must be an integer count when ",
        "denoising.method is 'dada2' (fractional mismatch rates are only ",
        "supported for vsearch/UNOISE merging).",
        call. = FALSE
      )
    }
    options(optimotu.pipeline.merge_max_mismatch = merging$max_mismatch)
  }
}

#' @rdname pipeline_options
#' @export
merge_min_overlap <- function() {
  stored <- getOption("optimotu.pipeline.merge_min_overlap", NULL)
  if (!is.null(stored)) {
    return(stored)
  }
  switch(
    denoising_method(),
    dada2 = 10L,
    unoise = 16L,
    10L
  )
}

#' @rdname pipeline_options
#' @export
merge_max_mismatch <- function() {
  stored <- getOption("optimotu.pipeline.merge_max_mismatch", NULL)
  if (!is.null(stored)) {
    return(stored)
  }
  switch(
    denoising_method(),
    dada2 = 1,
    unoise = 5,
    1
  )
}

#### filtering settings ####
#' @rdname parse_pipeline_options
#' @keywords internal
parse_filter_options <- function(pipeline_options) {
  checkmate::assert_list(pipeline_options$filtering, null.ok = TRUE)
  paired_keys <- paired_filter_option_names
  merged_keys <- merged_filter_option_names
  all_keys <- c(paired_keys, merged_keys)
  if (is.null(pipeline_options$filtering)) {
    message(
      "No 'filtering' options given in 'pipeline_options.yaml'\n",
      "Using defaults."
    )
    return(invisible())
  }
  filtering <- unnest_yaml_list(pipeline_options$filtering)
  checkmate::assert_names(names(filtering), subset.of = all_keys)
  present_paired <- intersect(names(filtering), paired_keys)
  present_merged <- intersect(names(filtering), merged_keys)
  if (do_unoise()) {
    if (length(present_paired) > 0L) {
      warning(
        "Ignoring paired-read filtering option(s) ",
        paste(present_paired, collapse = ", "),
        " because denoising method is 'unoise' (merged-read filtering is used).",
        call. = FALSE
      )
    }
    merged <- filtering[present_merged]
    if ("maxEE" %in% names(merged)) {
      checkmate::assert_number(
        merged$maxEE,
        lower = 0,
        finite = TRUE,
        null.ok = TRUE
      )
    }
    if ("maxEE_rate" %in% names(merged)) {
      checkmate::assert_number(
        merged$maxEE_rate,
        lower = 0,
        upper = 1,
        finite = TRUE,
        null.ok = TRUE
      )
    }
    if ("maxNs" %in% names(merged)) {
      checkmate::assert_count(merged$maxNs, null.ok = TRUE)
    }
    if ("maxLen" %in% names(merged)) {
      checkmate::assert_count(merged$maxLen, positive = TRUE, null.ok = TRUE)
    }
    if ("minLen" %in% names(merged)) {
      checkmate::assert_count(merged$minLen, null.ok = TRUE)
    }
    options(
      optimotu.pipeline.merged_filter_options = stats::update(
        merged_filter_options(),
        merged
      )
    )
  } else {
    if (length(present_merged) > 0L) {
      warning(
        "Ignoring merged-read filtering option(s) ",
        paste(present_merged, collapse = ", "),
        " because denoising method is 'dada2' (paired-read filtering is used).",
        call. = FALSE
      )
    }
    checkmate::assert_number(
      filtering$maxEE_R1,
      lower = 0,
      finite = TRUE,
      null.ok = TRUE
    )
    checkmate::assert_number(
      filtering$maxEE_R2,
      lower = 0,
      finite = TRUE,
      null.ok = TRUE
    )
    options(
      optimotu.pipeline.dada2_maxEE = stats::update(dada2_maxEE(), filtering)
    )
  }
}

#' @rdname pipeline_options
#' @export
dada2_maxEE <- function() {
  getOption("optimotu.pipeline.dada2_maxEE", dada2_filter_options(2, 2))
}

#### tag_jump settings ####

#' @rdname parse_pipeline_options
#' @keywords internal
parse_uncross_options <- function(pipeline_options) {
  tag_jump <- pipeline_options$tag_jump
  checkmate::assert(
    checkmate::check_null(tag_jump),
    checkmate::check_list(tag_jump),
    checkmate::check_false(tag_jump)
  )
  if (!is.null(tag_jump) && is.list(tag_jump)) {
    tag_jump <- unnest_yaml_list(tag_jump)
    checkmate::assert_names(names(tag_jump), subset.of = c("f", "p"))
    checkmate::assert_number(
      tag_jump$f,
      lower = 0,
      upper = 1,
      finite = TRUE,
      null.ok = TRUE
    )
    checkmate::assert_number(
      tag_jump$p,
      lower = 0,
      finite = TRUE,
      null.ok = TRUE
    )
    options(
      optimotu.pipeline.do_tag_jump = TRUE,
      optimotu.pipeline.tag_jump_f = tag_jump$f,
      optimotu.pipeline.tag_jump_p = tag_jump$p
    )
  }
}

#' @rdname pipeline_options
#' @export
do_tag_jump <- function() {
  getOption("optimotu.pipeline.do_tag_jump", FALSE)
}

#' @rdname pipeline_options
#' @export
tag_jump_f <- function() {
  getOption("optimotu.pipeline.tag_jump_f", 0.05)
}


#' @rdname pipeline_options
#' @export
tag_jump_p <- function() {
  getOption("optimotu.pipeline.tag_jump_p", 0.05)
}

#### LULU settings ####

#' @rdname parse_pipeline_options
#' @keywords internal
parse_lulu_options <- function(pipeline_options) {
  lulu_options <- pipeline_options$lulu
  if (is.null(lulu_options) || isFALSE(lulu_options)) {
    return()
  }
  checkmate::assert_list(lulu_options)
  lulu_options <- unnest_yaml_list(lulu_options)
  checkmate::assert_names(
    names(lulu_options),
    subset.of = c(
      "dist_type",
      "max_dist",
      "max_gap_length",
      "max_gap_total",
      "min_abundance_ratio",
      "min_cooccurrence_ratio",
      "use_mean_abundance_ratio",
      "dist_config"
    )
  )
  options(optimotu.pipeline.lulu_options = lulu_options)
  options(optimotu.pipeline.do_lulu = TRUE)
  ##### dist_type #####
  if ("dist_type" %in% names(lulu_options)) {
    checkmate::assert_string(lulu_options$dist_type)
    checkmate::assert_subset(
      lulu_options$dist_type,
      c("fraction", "score", "base")
    )
    options(optimotu.pipeline.lulu_dist_type = lulu_options$dist_type)
  }
  ##### max_dist #####
  if ("max_dist" %in% names(lulu_options)) {
    switch(
      lulu_options$dist_type,
      fraction = checkmate::assert_number(
        lulu_options$max_dist,
        lower = 0,
        upper = 1
      ),
      score = checkmate::assert_number(lulu_options$max_dist, lower = 0),
      base = checkmate::assert_integerish(lulu_options$max_dist, lower = 0)
    )
    options(optimotu.pipeline.lulu_max_dist = lulu_options$max_dist)
  }
  ##### max_gap_length #####
  if ("max_gap_length" %in% names(lulu_options)) {
    checkmate::assert_integerish(lulu_options$max_gap_length, lower = 0)
    options(optimotu.pipeline.lulu_max_gap_length = lulu_options$max_gap_length)
  }
  ##### max_gap_total #####
  if ("max_gap_total" %in% names(lulu_options)) {
    checkmate::assert(
      checkmate::check_integerish(lulu_options$max_gap_total, lower = 0),
      checkmate::check_number(lulu_options$max_gap_total, lower = 0, upper = 1)
    )
    options(optimotu.pipeline.lulu_max_gap_total = lulu_options$max_gap_total)
  }
  ##### min_abundance_ratio #####
  if ("min_abundance_ratio" %in% names(lulu_options)) {
    checkmate::assert_number(
      lulu_options$min_abundance_ratio,
      lower = 0,
      upper = 1
    )
    options(
      optimotu.pipeline.lulu_min_abundance_ratio = lulu_options$min_abundance_ratio
    )
  }
  ##### min_cooccurrence_ratio #####
  if ("min_cooccurrence_ratio" %in% names(lulu_options)) {
    checkmate::assert_number(
      lulu_options$min_cooccurrence_ratio,
      lower = 0,
      upper = 1
    )
    options(
      optimotu.pipeline.lulu_min_cooccurrence_ratio = lulu_options$min_cooccurrence_ratio
    )
  }
  ##### use_mean_abundance_ratio #####
  if ("use_mean_abundance_ratio" %in% names(lulu_options)) {
    checkmate::assert_flag(lulu_options$use_mean_abundance_ratio)
    options(
      optimotu.pipeline.lulu_use_mean_abundance_ratio = lulu_options$use_mean_abundance_ratio
    )
  }
  ##### dist_config #####
  lulu_dist_config <- resolve_section_dist_config(lulu_options$dist_config)
  options(optimotu.pipeline.lulu_dist_config = lulu_dist_config)
}

#' @rdname pipeline_options
#' @export
do_lulu <- function() {
  getOption("optimotu.pipeline.do_lulu", FALSE)
}

#' @rdname pipeline_options
#' @export
lulu_dist_type <- function() {
  getOption("optimotu.pipeline.lulu_dist_type", "fraction")
}

#' @rdname pipeline_options
#' @export
lulu_max_dist <- function() {
  getOption("optimotu.pipeline.lulu_max_dist", 0.1)
}

#' @rdname pipeline_options
#' @export
lulu_max_gap_length <- function() {
  getOption("optimotu.pipeline.lulu_max_gap_length", 1)
}

#' @rdname pipeline_options
#' @export
lulu_max_gap_total <- function() {
  getOption("optimotu.pipeline.lulu_max_gap_total", 3)
}


#' @rdname pipeline_options
#' @export
lulu_min_abundance_ratio <- function() {
  getOption("optimotu.pipeline.lulu_min_abundance_ratio", 1)
}

#' @rdname pipeline_options
#' @export
lulu_min_cooccurrence_ratio <- function() {
  getOption("optimotu.pipeline.lulu_min_cooccurrence_ratio", 1)
}

#' @rdname pipeline_options
#' @export
lulu_use_mean_abundance_ratio <- function() {
  getOption("optimotu.pipeline.lulu_use_mean_abundance_ratio", FALSE)
}


#' @rdname pipeline_options
#' @export
lulu_dist_config <- function() {
  getOption("optimotu.pipeline.lulu_dist_config", cluster_dist_config())
}


#### amplicon model settings ####

#' @rdname parse_pipeline_options
#' @export
#' @param pipeline_options (`list`) the pipeline options to parse
parse_amplicon_model_options <- function(pipeline_options) {
  amplicon_model_options <- pipeline_options$amplicon_model
  if (is.null(amplicon_model_options) || isFALSE(amplicon_model_options)) {
    return()
  }

  checkmate::assert_list(amplicon_model_options)
  amplicon_model_options <- unnest_yaml_list(amplicon_model_options)
  checkmate::assert_names(
    names(amplicon_model_options),
    must.include = "model_type"
  )
  ##### amplicon_model_type #####
  checkmate::assert_string(amplicon_model_options$model_type)
  checkmate::assert_subset(
    amplicon_model_options$model_type,
    c("CM", "HMM", "none")
  )
  options(
    optimotu.pipeline.amplicon_model_type = amplicon_model_options$model_type
  )

  if (!identical(amplicon_model_type, "none")) {
    # #### seed_aln ####
    # if ("seed_aln" %in% names(amplicon_model_options)) {
    #   seed_aln <<- amplicon_model_options$seed_aln
    #   checkmate::assert_file_exists(seed_aln, "r")
    #   do_generate_model <<- TRUE
    # }

    #### model_file ####
    checkmate::assert_names(
      names(amplicon_model_options),
      must.include = "model_file"
    )
    checkmate::assert_string(amplicon_model_options$model_file)
    # if (isFALSE(do_generate_model)) {
    checkmate::assert_file_exists(amplicon_model_options$model_file, "r")
    # } else {
    #   checkmate::assert_path_for_output(model_file, overwrite = TRUE)
    # }
    options(
      optimotu.pipeline.amplicon_model_file = amplicon_model_options$model_file
    )

    #### amplicon model filtering settings ####
    if (!is.null(amplicon_model_options$model_filter)) {
      parse_amplicon_model_filter_options(amplicon_model_options$model_filter)
    }

    #### amplicon alignment settings ###
    if (!is.null(amplicon_model_options$model_align)) {
      checkmate::assert_flag(amplicon_model_options$model_align)
      options(
        optimotu.pipeline.do_model_align = amplicon_model_options$model_align
      )
    }

    #### NuMt detection settings ####
    if ("numt_filter" %in% names(amplicon_model_options)) {
      checkmate::assert_logical(amplicon_model_options$numt_filter)
      if (isTRUE(amplicon_model_options$numt_filter)) {
        if (!identical(amplicon_model_type(), "HMM")) {
          stop("NuMt filter is only valid when HMM alignment is used")
        }
        options(optimotu.pipeline.do_numt_filter = TRUE)
      }
    }
  }
}

parse_amplicon_model_filter_options <- function(filter_options) {
  options(optimotu.pipeline.do_model_filter = TRUE)

  checkmate::assert_list(filter_options, min.len = 1)
  checkmate::assert_names(
    names(filter_options),
    subset.of = c("max_model_start", "min_model_end", "min_model_score")
  )
  if ("max_model_start" %in% names(filter_options)) {
    checkmate::assert_number(filter_options$max_model_start)
    options(optimotu.pipeline.max_model_start = filter_options$max_model_start)
  }

  if ("min_model_end" %in% names(filter_options)) {
    checkmate::assert_number(filter_options$min_model_end)
    options(optimotu.pipeline.min_model_end = filter_options$min_model_end)
  }

  if ("min_model_score" %in% names(filter_options)) {
    checkmate::assert_number(filter_options$min_model_score)
    options(optimotu.pipeline.min_model_score = filter_options$min_model_score)
  }
}

#' @rdname pipeline_options
#' @export
amplicon_model_type <- function() {
  getOption("optimotu.pipeline.amplicon_model_type", "none")
}

#' @rdname pipeline_options
#' @export
amplicon_model_file <- function() {
  getOption("optimotu.pipeline.amplicon_model_file")
}

#' @rdname pipeline_options
#' @export
max_model_start <- function() {
  getOption("optimotu.pipeline.max_model_start", Inf)
}

#' @rdname pipeline_options
#' @export
min_model_end <- function() {
  getOption("optimotu.pipeline.min_model_end", -Inf)
}

#' @rdname pipeline_options
#' @export
min_model_score <- function() {
  getOption("optimotu.pipeline.min_model_score", -Inf)
}

#' @rdname pipeline_options
#' @export
do_model_filter <- function() {
  getOption("optimotu.pipeline.do_model_filter", FALSE)
}

#' @rdname pipeline_options
#' @export
do_model_align <- function() {
  getOption("optimotu.pipeline.do_model_align", FALSE)
}

#' @rdname pipeline_options
#' @export
do_model_align_only <- function() {
  do_model_align() && !do_model_filter()
}

#' @rdname pipeline_options
#' @export
do_model_filter_only <- function() {
  do_model_filter() && !do_model_align()
}

#' @rdname pipeline_options
#' @export
do_model_both <- function() {
  do_model_filter() && do_model_align()
}

#' @rdname pipeline_options
#' @export
do_numt_filter <- function() {
  getOption("optimotu.pipeline.do_numt_filter", FALSE)
}

#### control sequence settings ####

#' @rdname parse_pipeline_options
#' @export
#' @param pipeline_options (`list`) the pipeline options to parse
parse_control_options <- function(pipeline_options) {
  if (!is.null(pipeline_options$control)) {
    checkmate::assert_list(pipeline_options$control)
    control_options <- unnest_yaml_list(pipeline_options$control)
    checkmate::assert_names(
      names(control_options),
      subset.of = c("spike", "positive")
    )

    if ("spike" %in% names(control_options)) {
      checkmate::assert(
        checkmate::check_file_exists(control_options$spike),
        checkmate::check_flag(control_options$spike, null.ok = TRUE)
      )
      if (isTRUE(control_options$spike)) {
        stop(
          "Option 'control':'spike' should be a file path, evaluate to ",
          "FALSE, or be left blank"
        )
      }
      if (is.character(control_options$spike)) {
        options(optimotu.pipeline.spike_file = control_options$spike)
      }
    }
    if ("positive" %in% names(control_options)) {
      checkmate::assert(
        checkmate::check_file_exists(control_options$positive),
        checkmate::check_flag(control_options$positive, null.ok = TRUE)
      )
      if (isTRUE(control_options$positive)) {
        stop(
          "Option 'control':'positive' should be a file path, evaluate ",
          "to FALSE, or be left blank"
        )
      }
      if (is.character(control_options$positive)) {
        options(optimotu.pipeline.pos_control_file = control_options$positive)
      }
    }
  }
}

#' @rdname pipeline_options
#' @export
spike_file <- function() {
  getOption("optimotu.pipeline.spike_file", NULL)
}

#' @rdname pipeline_options
#' @export
do_spike <- function() {
  is.character(spike_file())
}

#' @rdname pipeline_options
#' @export
pos_control_file <- function() {
  getOption("optimotu.pipeline.pos_control_file", NULL)
}

#' @rdname pipeline_options
#' @export
do_pos_control <- function() {
  is.character(pos_control_file())
}

#### taxonomic assignment settings ####

#' @rdname parse_pipeline_options
#' @export
parse_taxonomy_options <- function(pipeline_options) {
  taxonomy_options <- pipeline_options$taxonomy
  if (is.null(taxonomy_options) || isFALSE(taxonomy_options)) {
    return()
  }
  checkmate::assert_list(taxonomy_options)
  taxonomy_options <- unnest_yaml_list(taxonomy_options)
  classifier_names <- c("protax", "sintax", "bayesant", "epa")
  selected_classifier <- intersect(
    names(taxonomy_options),
    classifier_names
  )
  if (length(selected_classifier) > 1) {
    stop(
      "Only one of options 'taxonomy:protax', 'taxonomy:sintax',",
      " 'taxonomy:bayesant' and 'taxonomy:epa' may be given.",
      "(file: pipeline_options.yaml)"
    )
  }
  if (length(selected_classifier) == 0) {
    stop(
      "No classifier selected. Please select one of the following ",
      "classifiers:\n",
      "  - protax\n",
      "  - sintax\n",
      "  - bayesant\n",
      "  - epa\n",
      "(file: pipeline_options.yaml)"
    )
  }
  switch(
    selected_classifier,
    protax = parse_protax_options(taxonomy_options$protax),
    sintax = parse_sintax_options(taxonomy_options$sintax),
    bayesant = parse_bayesant_options(taxonomy_options$bayesant),
    epa = parse_epa_options(taxonomy_options$epa)
  )

  if ("ranks" %in% names(taxonomy_options)) {
    parse_taxonomy_ranks(taxonomy_options$ranks)
  }
}

##### taxonomic ranks #####

#' @rdname parse_pipeline_options
#' @export
#' @param rank_options (`list` or `character` vector) the taxonomic ranks to
#' use in the pipeline, in order from most inclusive (e.g., kingdom) to least
#' inclusive (e.g., species). Values may be either named or unnamed. When
#' named, the name is taken to be the rank, and the value is the "in-group"
#' taxon at that rank, i.e. the taxon for which results are desired. When
#' unnamed, the value is taken to be the rank. Example: `list(kingdom =
#' "Fungi", "phylum", "class", "order", "family", "genus", "species")`
#' @return `NULL`.  This function is called for its side effect, which is to
#' configure global options.
parse_taxonomy_ranks <- function(rank_options) {
  checkmate::assert(
    checkmate::check_list(
      rank_options,
      types = c("character", "list"),
      min.len = 1
    ),
    checkmate::check_character(
      rank_options,
      unique = TRUE,
      min.len = 1
    )
  )
  KNOWN_TAXA <- purrr::keep(
    rank_options,
    ~ dplyr::cumall(checkmate::test_list(.x))
  ) |>
    unlist()
  UNKNOWN_RANKS <- purrr::discard(
    rank_options,
    ~ dplyr::cumall(checkmate::test_list(.x))
  ) |>
    unlist()
  if (length(UNKNOWN_RANKS) == 0 || !is.null(names(UNKNOWN_RANKS))) {
    stop(
      "Option 'taxonomy':'ranks' should start from the most inclusive rank ",
      "(e.g. kingdom)\n",
      "  and continue to the least inclusive rank (e.g. species).  Optionally",
      " the first\n",
      "  rank(s) may be defined (e.g. '- kingdom: Fungi') but subsequent ranks",
      " must be undefined (e.g. '- class').\n",
      "  (file: pipeline_options.yaml)"
    )
  }
  if (length(KNOWN_TAXA) == 0) {
    KNOWN_TAXA <- c(rootrank = "root")
  } else {
    options(optimotu.pipeline.do_outgroup = TRUE)
  }
  set_known_ranks(names(KNOWN_TAXA))
  set_known_taxa(unname(KNOWN_TAXA))
  set_tax_ranks(c(known_ranks(), UNKNOWN_RANKS))
}

##### protax #####

#' @rdname parse_pipeline_options
#' @param protax_options (`list`) the protax options to parse
parse_protax_options <- function(protax_options) {
  options(optimotu.pipeline.do_protax = TRUE)
  checkmate::assert_list(protax_options)
  protax_options <- unnest_yaml_list(protax_options)
  if (length(protax_options) > 0L) {
    checkmate::assert_names(
      names(protax_options),
      subset.of = c("aligned", "location", "ranks", "added_reference")
    )
  }

  ##### protax version #####
  if ("aligned" %in% names(protax_options)) {
    checkmate::assert_flag(protax_options$aligned)
    if (protax_options$aligned && !do_model_align()) {
      stop(
        "Aligned Protax (taxonomy: protax: aligned: true) requires model",
        "alignment to be enabled (amplicon_model: model_align: true).\n",
        "(file: pipeline_options.yaml)"
      )
    }
    options(optimotu.pipeline.protax_aligned = protax_options$aligned)
  } else {
    message("Using unaligned protax by default.")
  }

  ##### protax location #####
  if ("location" %in% names(protax_options)) {
    checkmate::assert_directory_exists(protax_options$location)
    options("optimotu.pipeline.protax_location" = protax_options$location)
  } else {
    message("Using default protax directory: ", protax_location())
  }

  if ("ranks" %in% names(protax_options)) {
    parse_taxonomy_ranks(protax_options$ranks)
  }

  ##### added_reference #####
  if ("added_reference" %in% names(protax_options)) {
    apply_added_reference(protax_options$added_reference)
  }
}

#' @rdname pipeline_options
#' @export
do_protax <- function() {
  getOption("optimotu.pipeline.do_protax", FALSE)
}

#' @rdname pipeline_options
#' @export
protax_aligned <- function() {
  getOption("optimotu.pipeline.protax_aligned", FALSE)
}

#' @rdname pipeline_options
#' @export
protax_unaligned <- function() {
  do_protax() && !protax_aligned()
}

#' @rdname pipeline_options
#' @export
protax_location <- function() {
  getOption("optimotu.pipeline.protax_location", "protaxFungi")
}

##### sintax #####

#' @rdname parse_pipeline_options
#' @param sintax_options (`list`) the sintax options to parse
parse_sintax_options <- function(sintax_options) {
  checkmate::assert_list(sintax_options)
  sintax_options <- unnest_yaml_list(sintax_options)
  checkmate::assert_file_exists(sintax_options$reftax, "r")
  options(
    optimotu.pipeline.do_sintax = TRUE,
    optimotu.pipeline.sintax_ref = sintax_options$reftax
  )
}

#' @rdname pipeline_options
#' @export
do_sintax <- function() {
  getOption("optimotu.pipeline.do_sintax", FALSE)
}

#' @rdname pipeline_options
#' @export
sintax_ref <- function() {
  getOption("optimotu.pipeline.sintax_ref")
}

##### bayesant #####

#' @rdname parse_pipeline_options
#' @param bayesant_options (`list`) the bayesant options to parse
parse_bayesant_options <- function(bayesant_options) {
  checkmate::assert_list(bayesant_options)
  bayesant_options <- unnest_yaml_list(bayesant_options)
  checkmate::assert(
    checkmate::check_file_exists(bayesant_options$reftax, "r"),
    checkmate::check_file_exists(bayesant_options$model, "r")
  )
  checkmate::assert_flag(bayesant_options$aligned, null.ok = TRUE)
  options(
    optimotu.pipeline.do_bayesant = TRUE,
    optimotu.pipeline.bayesant_aligned = bayesant_options$aligned,
    optimotu.pipeline.bayesant_ref = bayesant_options$reftax,
    optimotu.pipeline.bayesant_model = bayesant_options$model
  )
}

#' @rdname pipeline_options
#' @export
do_bayesant <- function() {
  getOption("optimotu.pipeline.do_bayesant", FALSE)
}

#' @rdname pipeline_options
#' @export
bayesant_aligned <- function() {
  getOption("optimotu.pipeline.bayesant_aligned", FALSE)
}

#' @rdname pipeline_options
#' @export
bayesant_ref <- function() {
  getOption("optimotu.pipeline.bayesant_ref")
}

#' @rdname pipeline_options
#' @export
bayesant_model <- function() {
  getOption("optimotu.pipeline.bayesant_model")
}

##### epa-ng #####

#' @rdname parse_pipeline_options
#' @param epa_options (`list`) the epa options to parse
#' @export
parse_epa_options <- function(epa_options) {
  checkmate::assert_list(epa_options)
  epa_options <- unnest_yaml_list(epa_options)
  checkmate::assert_names(
    names(epa_options),
    must.include = c("reference", "taxonomy", "tree", "params", "outgroup")
  )
  checkmate::assert_file_exists(epa_options$reference, "r")
  checkmate::assert_file_exists(epa_options$taxonomy, "r")
  checkmate::assert_file_exists(epa_options$tree, "r")
  checkmate::assert_string(epa_options$params)
  checkmate::assert_character(epa_options$outgroup, min.chars = 1, min.len = 1)
  options(
    optimotu.pipeline.do_epa = TRUE,
    optimotu.pipeline.epa_ref = epa_options$reference,
    optimotu.pipeline.epa_taxonomy = epa_options$taxonomy,
    optimotu.pipeline.epa_tree = epa_options$tree,
    optimotu.pipeline.epa_params = epa_options$params,
    optimotu.pipeline.epa_outgroup = epa_options$outgroup
  )
}

#' @rdname pipeline_options
#' @export
do_epa <- function() {
  getOption("optimotu.pipeline.do_epa", FALSE)
}

#' @rdname pipeline_options
#' @export
epa_ref <- function() {
  getOption("optimotu.pipeline.epa_ref")
}

#' @rdname pipeline_options
#' @export
epa_taxonomy <- function() {
  getOption("optimotu.pipeline.epa_taxonomy")
}

#' @rdname pipeline_options
#' @export
epa_tree <- function() {
  getOption("optimotu.pipeline.epa_tree")
}

#' @rdname pipeline_options
#' @export
epa_params <- function() {
  getOption("optimotu.pipeline.epa_params")
}

#' @rdname pipeline_options
#' @export
epa_outgroup <- function() {
  getOption("optimotu.pipeline.epa_outgroup")
}

#### Outgroup reference settings ####

#' @rdname parse_pipeline_options
#' @export
parse_outgroup_options <- function(pipeline_options) {
  # if outgroup_reference is not set at all, then we defer to the setting
  # during parsing of taxonomic ranks
  if (is.null(pipeline_options$outgroup_reference)) {
    return()
  }
  # if outgroup_reference is explicitly FALSE, then we set the option.
  # this may lead to cases where there is a defined ingroup, but no outgroup
  # reference.
  if (isFALSE(pipeline_options$outgroup_reference)) {
    if (do_outgroup()) {
      warning(
        "Taxonomic rank definitions imply the possibility of outgroup ",
        "sequences but no outgroup reference file was provided.\n",
        "(file: pipeline_options.yaml)"
      )
    }
    options(
      optimotu.pipeline.do_outgroup = FALSE
    )
    return()
  }
  outgroup_options <- unnest_yaml_list(pipeline_options$outgroup_reference)
  checkmate::assert_names(
    names(outgroup_options),
    must.include = "sequences",
    subset.of = c("sequences", "taxonomy")
  )
  checkmate::assert_file_exists(outgroup_options$sequences, "r")
  options(
    optimotu.pipeline.outgroup_reference = outgroup_options$sequences,
    optimotu.pipeline.do_outgroup = TRUE
  )
  if (!is.null(outgroup_options$taxonomy)) {
    checkmate::assert_file_exists(outgroup_options$taxonomy, "r")
    options(
      optimotu.pipeline.outgroup_taxonomy = outgroup_options$taxonomy
    )
  }
}

#' @rdname pipeline_options
#' @export
do_outgroup <- function() {
  getOption("optimotu.pipeline.do_outgroup", FALSE)
}

#' @rdname pipeline_options
#' @export
outgroup_reference <- function() {
  getOption("optimotu.pipeline.outgroup_reference")
}

#' @rdname pipeline_options
#' @export
outgroup_taxonomy <- function() {
  getOption("optimotu.pipeline.outgroup_taxonomy")
}

#### clustering settings ####

#' Default clustering job-size thresholds for a distance method
#'
#' Hamming and USEARCH are treated as fast methods; WFA2, Edlib, and hybrid
#' alignment are ~100x slower per comparison, so both thresholds are 100x
#' smaller.
#'
#' @param method (`character` scalar) `dist_config` method name
#' @return named `list` with `min_parallel_ops` and `max_batch_ops`
#' @keywords internal
cluster_ops_defaults <- function(method) {
  checkmate::assert_string(method)
  if (method %in% c("wfa2", "edlib", "hybrid")) {
    list(min_parallel_ops = 1e4, max_batch_ops = 1e8)
  } else {
    list(min_parallel_ops = 1e6, max_batch_ops = 1e10)
  }
}

message_cluster_ops_defaults <- function(
  method,
  min_parallel_ops,
  max_batch_ops,
  min_is_default,
  max_is_default
) {
  if (!min_is_default && !max_is_default) {
    return(invisible())
  }
  fmt <- function(x) format(x, scientific = TRUE, digits = 1)
  lines <- paste0(
    "Using default clustering job sizing for dist_config method '",
    method,
    "':"
  )
  if (min_is_default) {
    lines <- c(
      lines,
      paste0(
        "  min_parallel_ops = ",
        fmt(min_parallel_ops),
        "  (parallel-efficiency cutoff)"
      )
    )
  }
  if (max_is_default) {
    lines <- c(
      lines,
      paste0(
        "  max_batch_ops = ",
        fmt(max_batch_ops),
        "  (target ~10-20 min per batch)"
      )
    )
  }
  lines <- c(
    lines,
    "Set clustering.min_parallel_ops and/or clustering.max_batch_ops",
    "in pipeline_options.yaml to override."
  )
  message(paste(lines, collapse = "\n"))
}

#' @rdname parse_pipeline_options
#' @keywords internal
parse_cluster_options <- function(pipeline_options) {
  checkmate::assert_list(pipeline_options$clustering, null.ok = TRUE)
  if (is.null(pipeline_options$clustering)) {
    checkmate::assert_string(
      pipeline_options$cluster_thresholds,
      null.ok = TRUE
    )
    if (!is.null(pipeline_options$cluster_thresholds)) {
      checkmate::assert_file_exists(
        pipeline_options$cluster_thresholds,
        access = "r"
      )
      options(
        optimotu.pipeline.clustering_thresholds = pipeline_options$cluster_thresholds
      )
    } else {
      message(
        "No clustering options given in 'pipeline_options.yaml'\n",
        "Using defaults."
      )
    }
  } else {
    clustering <- unnest_yaml_list(pipeline_options$clustering)
    checkmate::assert_names(
      names(clustering),
      subset.of = c(
        "thresholds",
        "measure",
        "dist_config",
        "force_denovo",
        "min_parallel_ops",
        "max_batch_ops"
      )
    )
    parse_cluster_thresholds(clustering$thresholds)
    checkmate::assert_string(clustering$measure, null.ok = TRUE)
    if (!is.null(clustering$measure)) {
      checkmate::assert_subset(
        clustering$measure,
        c("MCC", "RI", "ARI", "FMI", "MI", "AMI", "FM")
      )
    }
    dist_config <- resolve_section_dist_config(clustering$dist_config)
    if (
      !is.null(clustering$force_denovo) &&
        !isFALSE(clustering$force_denovo)
    ) {
      checkmate::assert_character(clustering$force_denovo, unique = TRUE)
      checkmate::assert_subset(clustering$force_denovo, unknown_ranks())
    } else {
      clustering$force_denovo <- character(0)
    }
    min_is_default <- is.null(clustering$min_parallel_ops)
    max_is_default <- is.null(clustering$max_batch_ops)
    method <- dist_config$method
    ops_defaults <- cluster_ops_defaults(method)
    min_parallel_ops <- clustering$min_parallel_ops
    if (min_is_default) {
      min_parallel_ops <- ops_defaults$min_parallel_ops
    }
    checkmate::assert_number(min_parallel_ops, lower = 1)
    max_batch_ops <- clustering$max_batch_ops
    if (max_is_default) {
      max_batch_ops <- ops_defaults$max_batch_ops
    }
    checkmate::assert_number(max_batch_ops, lower = 1)
    if (max_batch_ops < min_parallel_ops) {
      stop(
        "clustering max_batch_ops must be greater than or equal to ",
        "min_parallel_ops (file: pipeline_options.yaml)"
      )
    }
    message_cluster_ops_defaults(
      method = method,
      min_parallel_ops = min_parallel_ops,
      max_batch_ops = max_batch_ops,
      min_is_default = min_is_default,
      max_is_default = max_is_default
    )
    options(
      optimotu.pipeline.clustering_measure = clustering$measure,
      optimotu.pipeline.clustering_dist_config = dist_config,
      optimotu.pipeline.clustering_force_denovo = clustering$force_denovo,
      optimotu.pipeline.clustering_min_parallel_ops = min_parallel_ops,
      optimotu.pipeline.clustering_max_batch_ops = max_batch_ops
    )
  }
}

#' @rdname parse_pipeline_options
#' @keywords internal
parse_cluster_thresholds <- function(thresh_opts) {
  if (is.null(thresh_opts)) {
    message(
      "No threshold options given in 'pipeline_options.yaml'\n",
      "  Using default thresholds file:",
      cluster_thresholds(),
      "\n"
    )
    return()
  }
  if (identical(names(thresh_opts), "file", )) {
    thresh_opts <- thresh_opts$file
  }
  if (is.character(thresh_opts)) {
    checkmate::assert_file_exists(thresh_opts, "r")
    options(
      optimotu.pipeline.clustering_thresholds = thresh_opts
    )
  } else {
    checkmate::assert_list(thresh_opts)
    checkmate::assert_names(
      names(thresh_opts),
      subset.of = c(
        "file",
        "train_data",
        "min_conf",
        "dist_max",
        "dist_step",
        "min_taxa",
        "min_refseq"
      ),
      must.include = "train_data"
    )
    checkmate::assert_string(thresh_opts$train_data)
    checkmate::assert(
      checkmate::check_file_exists(thresh_opts$train_data, "r"),
      checkmate::check_choice(
        thresh_opts$train_data,
        c("self", "reference")
      )
    )
    checkmate::assert_number(
      thresh_opts$min_conf,
      lower = 0,
      upper = 1,
      null.ok = TRUE
    )
    checkmate::assert_number(
      thresh_opts$dist_max,
      lower = 0,
      upper = 1,
      null.ok = TRUE
    )
    checkmate::assert_number(
      thresh_opts$dist_step,
      lower = 0,
      upper = 1,
      null.ok = TRUE
    )
    if (
      !is.null(thresh_opts$dist_max) &&
        !is.null(thresh_opts$dist_step) &&
        thresh_opts$dist_step > thresh_opts$dist_max
    ) {
      stop(
        "Option 'dist_step' must be less than or equal to option",
        " 'dist_max'.\n",
        "(file: pipeline_options.yaml)"
      )
    }
    checkmate::assert_integerish(
      thresh_opts$min_taxa,
      lower = 2,
      null.ok = TRUE
    )
    checkmate::assert_integerish(
      thresh_opts$min_refseq,
      lower = 2,
      null.ok = TRUE
    )
    if (
      !is.null(thresh_opts$min_taxa) &&
        !is.null(thresh_opts$min_refseq) &&
        thresh_opts$min_taxa > thresh_opts$min_refseq
    ) {
      stop(
        "Option 'min_refseq' must be greater than or equal to option",
        " 'min_taxa' (and ideally should be >=2x 'min_taxa').\n",
        "(file: pipeline_options.yaml)"
      )
    }
    options(
      optimotu.pipeline.clustering_thresholds = thresh_opts$file,
      optimotu.pipeline.do_optimize_thresholds = TRUE,
      optimotu.pipeline.do_optimize_thresholds_self = thresh_opts$train_data ==
        "self",
      optimotu.pipeline.do_optimize_thresholds_reference = thresh_opts$train_data ==
        "reference",
      optimotu.pipeline.do_optimize_thresholds_file = !thresh_opts$train_data %in%
        c("self", "reference"),
      optimotu.pipeline.optimize_thresholds_file = if (
        thresh_opts$train_data %in% c("self", "reference")
      ) {
        NULL
      } else {
        thresh_opts$train_data
      },
      optimotu.pipeline.clustering_min_conf = thresh_opts$min_conf,
      optimotu.pipeline.clustering_dist_max = thresh_opts$dist_max,
      optimotu.pipeline.clustering_dist_step = thresh_opts$dist_step,
      optimotu.pipeline.clustering_min_taxa = thresh_opts$min_taxa,
      optimotu.pipeline.clustering_min_refseq = thresh_opts$min_refseq
    )
  }
}

#' @rdname parse_pipeline_options
#' @keywords internal
parse_top_dist_config <- function(pipeline_options) {
  dist_config <- pipeline_options$dist_config
  if (is.null(dist_config)) {
    return(invisible())
  }
  if (is.character(dist_config)) {
    dist_config <- list(method = dist_config)
  } else {
    checkmate::assert_list(dist_config)
    dist_config <- unnest_yaml_list(dist_config)
  }
  options(optimotu.pipeline.top_dist_config = dist_config)
}

#' @rdname pipeline_options
#' @export
top_dist_config <- function() {
  getOption("optimotu.pipeline.top_dist_config", NULL)
}

#' Merge top-level and section-level dist_config lists
#'
#' Section keys always win. Top-level keys are inherited only when the section
#' omits `method` or specifies the same `method`, so parameters for one
#' distance method never leak into another.
#'
#' @param section_dist (`list`, `character`, or `NULL`) section-level config
#' @param top_dist (`list` or `NULL`) top-level config
#' @return (`list`) merged raw dist_config
#' @noRd
merge_dist_config_lists <- function(
  section_dist,
  top_dist = top_dist_config()
) {
  if (is.character(section_dist)) {
    section_dist <- list(method = section_dist)
  } else if (!is.null(section_dist)) {
    checkmate::assert_list(section_dist)
    section_dist <- unnest_yaml_list(section_dist)
  }
  if (is.null(top_dist) || length(top_dist) == 0L) {
    if (is.null(section_dist) || length(section_dist) == 0L) {
      return(list(method = "usearch"))
    }
    if (is.null(section_dist$method)) {
      section_dist$method <- "usearch"
    }
    return(section_dist)
  }
  if (is.null(section_dist) || length(section_dist) == 0L) {
    if (is.null(top_dist$method)) {
      top_dist$method <- "usearch"
    }
    return(top_dist)
  }
  section_method <- section_dist$method
  top_method <- top_dist$method
  if (
    is.null(section_method) ||
      (!is.null(top_method) && identical(section_method, top_method))
  ) {
    merged <- top_dist
    for (nm in names(section_dist)) {
      merged[[nm]] <- section_dist[[nm]]
    }
    if (is.null(merged$method)) {
      merged$method <- "usearch"
    }
    return(merged)
  }
  # Different methods: do not inherit top-level parameters
  if (is.null(section_dist$method)) {
    section_dist$method <- "usearch"
  }
  section_dist
}

#' Resolve a section dist_config against the top-level default
#' @param section_dist (`list`, `character`, or `NULL`)
#' @return result of [parse_dist_config()]
#' @noRd
resolve_section_dist_config <- function(section_dist) {
  parse_dist_config(merge_dist_config_lists(section_dist))
}

#' @rdname parse_pipeline_options
#' @keywords internal
parse_dist_config <- function(dist_config) {
  if (is.null(dist_config) || length(dist_config) == 0L) {
    dist_config <- list(method = "usearch")
  } else if (is.character(dist_config)) {
    dist_config <- list(method = dist_config)
  } else {
    checkmate::assert_list(dist_config)
    dist_config <- unnest_yaml_list(dist_config)
    if (is.null(dist_config$method)) {
      dist_config$method <- "usearch"
    }
  }
  dist_config <- c(list(quote(optimotu::dist_config)), dist_config)
  if (dist_config$method == "usearch" && !"usearch" %in% names(dist_config)) {
    dist_config$usearch <- find_usearch()
  }
  if (dist_config$method == "hamming" && !do_model_align()) {
    stop(
      "Hamming distance requires model alignment to be enabled ",
      "(amplicon_model: model_align: true).\n",
      "(file: pipeline_options.yaml)"
    )
  }
  eval(as.call(dist_config))
}

#' @rdname pipeline_options
#' @export
cluster_thresholds <- function() {
  getOption(
    "optimotu.pipeline.clustering_thresholds",
    "metadata/GSSP_thresholds.tsv"
  )
}

#' @rdname pipeline_options
#' @export
do_optimize_thresholds <- function() {
  getOption("optimotu.pipeline.do_optimize_thresholds", FALSE)
}

#' @rdname pipeline_options
#' @export
do_optimize_thresholds_self <- function() {
  getOption("optimotu.pipeline.do_optimize_thresholds_self", FALSE)
}

#' @rdname pipeline_options
#' @export
do_optimize_thresholds_reference <- function() {
  getOption("optimotu.pipeline.do_optimize_thresholds_reference", FALSE)
}

#' @rdname pipeline_options
#' @export
do_optimize_thresholds_file <- function() {
  getOption("optimotu.pipeline.do_optimize_thresholds_file", FALSE)
}

#' @rdname pipeline_options
#' @export
optimize_thresholds_file <- function() {
  getOption("optimotu.pipeline.optimize_thresholds_file", NULL)
}

#' @rdname pipeline_options
#' @export
cluster_measure <- function() {
  getOption("optimotu.pipeline.clustering_measure", NULL)
}

#' @rdname pipeline_options
#' @export
cluster_dist_config <- function() {
  getOption(
    "optimotu.pipeline.clustering_dist_config",
    optimotu::dist_usearch()
  )
}

#' @rdname pipeline_options
#' @export
cluster_force_denovo <- function() {
  getOption("optimotu.pipeline.clustering_force_denovo", character(0))
}

#' @rdname pipeline_options
#' @export
cluster_min_parallel_ops <- function() {
  getOption(
    "optimotu.pipeline.clustering_min_parallel_ops",
    cluster_ops_defaults(cluster_dist_config()$method)$min_parallel_ops
  )
}

#' @rdname pipeline_options
#' @export
cluster_max_batch_ops <- function() {
  getOption(
    "optimotu.pipeline.clustering_max_batch_ops",
    cluster_ops_defaults(cluster_dist_config()$method)$max_batch_ops
  )
}

#' @rdname pipeline_options
#' @export
cluster_min_conf <- function() {
  getOption("optimotu.pipeline.clustering_min_conf", 0.5)
}

#' @rdname pipeline_options
#' @export
cluster_dist_max <- function() {
  getOption("optimotu.pipeline.clustering_dist_max", 0.4)
}

#' @rdname pipeline_options
#' @export
cluster_dist_step <- function() {
  getOption("optimotu.pipeline.clustering_dist_step", 0.001)
}

#' @rdname pipeline_options
#' @export
cluster_min_taxa <- function() {
  getOption("optimotu.pipeline.clustering_min_taxa", 5)
}

#' @rdname pipeline_options
#' @export
cluster_min_refseq <- function() {
  getOption("optimotu.pipeline.clustering_min_refseq", 2 * cluster_min_taxa())
}

#### guilds settings ####
# TODO: this needs additional customization options

#' @rdname parse_pipeline_options
#' @keywords internal
parse_guilds_options <- function(pipeline_options) {
  if (!is.null(pipeline_options$guilds)) {
    checkmate::assert_flag(pipeline_options$guilds)
    options(
      optimotu.pipeline.do_guilds = pipeline_options$guilds
    )
  }
}

#' @rdname pipeline_options
#' @export
do_guilds <- function() {
  getOption("optimotu.pipeline.do_guilds", FALSE)
}

#### output / OTU table settings ####

#' Normalize and validate output format names from pipeline options
#' @param formats (`character` or `list`) format names from YAML
#' @return normalized lowercase `character` vector
#' @keywords internal
normalize_output_formats <- function(formats) {
  if (is.null(formats)) {
    return(c("rds", "tsv"))
  }
  if (is.list(formats) && !is.character(formats)) {
    if (all(vapply(formats, is.atomic, logical(1)))) {
      formats <- unlist(formats, use.names = FALSE)
    } else {
      formats <- unlist(unnest_yaml_list(formats), use.names = FALSE)
    }
  }
  checkmate::assert_character(formats, min.len = 1)
  formats <- unique(tolower(as.character(formats)))
  formats[formats == "qd"] <- "qdata"
  allowed <- c(
    "rds",
    "tsv",
    "csv",
    "xlsx",
    "fst",
    "feather",
    "parquet",
    "qs2",
    "qdata",
    "rdata",
    "rda"
  )
  bad <- setdiff(formats, allowed)
  if (length(bad) > 0L) {
    stop(
      "Unknown output format(s) in 'pipeline_options.yaml': ",
      paste(bad, collapse = ", "),
      ".\nAllowed formats: ",
      paste(allowed, collapse = ", "),
      call. = FALSE
    )
  }
  formats
}

#' Parse output options from pipeline_options.yaml
#' @param pipeline_options (`list`) parsed YAML options
#' @keywords internal
parse_output_options <- function(pipeline_options) {
  output_opts <- pipeline_options$output
  if (!is.null(output_opts) && is.list(output_opts)) {
    output_opts <- unnest_yaml_list(output_opts)
  }

  wide_table <- pipeline_options$wide_table
  dense_table <- pipeline_options$dense_table
  formats <- NULL

  if (!is.null(output_opts)) {
    if (!is.null(output_opts$wide_table)) {
      wide_table <- output_opts$wide_table
    }
    if (!is.null(output_opts$dense_table)) {
      dense_table <- output_opts$dense_table
    }
    if (!is.null(output_opts$formats)) {
      formats <- output_opts$formats
    }
  }

  if (!is.null(wide_table) && !is.null(dense_table)) {
    stop(
      "Only one of 'wide_table' and 'dense_table' may be set in ",
      "'pipeline_options.yaml'.\n",
      "Please remove one of them.",
      call. = FALSE
    )
  }

  checkmate::assert_flag(wide_table, null.ok = TRUE)
  checkmate::assert_flag(dense_table, null.ok = TRUE)

  wide_on <- isTRUE(wide_table) || isTRUE(dense_table)
  if (wide_on) {
    options(optimotu.pipeline.wide_table = TRUE)
  }

  options(
    optimotu.pipeline.output_formats = normalize_output_formats(formats)
  )
}

#' @rdname parse_pipeline_options
#' @export
parse_otu_table_options <- function(pipeline_options) {
  parse_output_options(pipeline_options)
}

#' @rdname pipeline_options
#' @export
do_wide_otu_table <- function() {
  getOption("optimotu.pipeline.wide_table", FALSE)
}

#' Output format names requested in pipeline_options.yaml
#'
#' Includes per-file formats and `rdata` when a bundled `.RData` file is
#' requested.
#'
#' @return `character` vector of format names (lowercase). Defaults to
#'   `c("rds", "tsv")` when not configured.
#' @export
output_formats <- function() {
  getOption("optimotu.pipeline.output_formats", c("rds", "tsv"))
}

#' Per-file tabular output formats (excludes `rdata`)
#'
#' @return `character` vector of format names for [write_tabular_outputs()]
#' @export
output_table_formats <- function() {
  formats <- output_formats()
  formats[!formats %in% c("rdata", "rda")]
}

#' Whether to write a bundled RData file of tabular outputs
#'
#' @return `logical` scalar
#' @export
do_output_rdata <- function() {
  any(output_formats() %in% c("rdata", "rda"))
}

#### rarefaction settings ####
#' @rdname parse_pipeline_options
#' @export
parse_rarefy_options <- function(pipeline_options) {
  rarefy_options <- pipeline_options$rarefy
  if (!is.null(rarefy_options)) {
    rarefy_options <- unnest_yaml_list(rarefy_options)
    checkmate::assert_names(
      names(rarefy_options),
      subset.of = c("numerator", "denominator", "number")
    )
    if ("number" %in% names(rarefy_options)) {
      if (any(c("numerator", "denominator") %in% names(rarefy_options))) {
        stop(
          "Option 'rarefy:number' cannot be used in conjunction with ",
          "option 'rarefy:numerator' or 'rarefy:denominator'.\n",
          "(file: pipeline_options.yaml)"
        )
      }
      checkmate::assert_integerish(rarefy_options$number, lower = 1)
      options(
        optimotu.pipeline.rarefy_number = rarefy_options$number
      )
    } else if ("numerator" %in% names(rarefy_options)) {
      if (!"denominator" %in% names(rarefy_options)) {
        stop(
          "Option 'rarefy:numerator' requires option 'rarefy:denominator'.\n",
          "(file: pipeline_options.yaml)"
        )
      }
      checkmate::assert_integerish(rarefy_options$numerator, lower = 1)
      checkmate::assert_integerish(rarefy_options$denominator, lower = 1)
      options(
        optimotu.pipeline.rarefy_numerator = rarefy_options$numerator,
        optimotu.pipeline.rarefy_denominator = rarefy_options$denominator
      )
    } else if ("denominator" %in% names(rarefy_options)) {
      stop(
        "Option 'rarefy:denominator' requires option 'rarefy:numerator'.\n",
        "(file: pipeline_options.yaml)"
      )
    } else {
      warning(
        "Empty rarefy options given in 'pipeline_options.yaml'\n",
        "No rarefaction will be performed."
      )
    }
  }
}

#' @rdname pipeline_options
#' @export
rarefy_number <- function() {
  getOption("optimotu.pipeline.rarefy_number", NULL)
}

#' @rdname pipeline_options
#' @export
rarefy_numerator <- function() {
  getOption("optimotu.pipeline.rarefy_numerator", NULL)
}

#' @rdname pipeline_options
#' @export
rarefy_denominator <- function() {
  getOption("optimotu.pipeline.rarefy_denominator", NULL)
}

#' @rdname pipeline_options
#' @export
do_rarefy <- function() {
  !is.null(rarefy_number()) || !is.null(rarefy_numerator())
}

#' Mapping variables for rarefaction
#' @param dots (`logical` scalar) whether to prefix the variable names with a
#' dot; if `TRUE`, the variable names will be prefixed with a dot, e.g.
#' `.numerator` instead of `numerator`
#' @return a `data.frame` giving the variables to map over for rarefaction;
#' when rarefaction is to be performed (as determined by `do_rarefy()`) then
#' these will always include `.rarefy_text` (`character)`), and will also
#' include either `.numerator` (`integer`) and `.denominator` (`integer`) or
#' `.number` (`integer`).  Alternatively, if `do_rarefy()` returns `FALSE`, an
#' empty data.frame.
#' @export
rarefy_meta <- function(dots = TRUE) {
  # avoid R CMD check NOTE about global variables due to NSE
  numerator <- denominator <- number <- NULL
  if (do_rarefy()) {
    (if (is.null(rarefy_number())) {
      tibble::tibble(
        numerator = rarefy_numerator(),
        denominator = rarefy_denominator(),
        rarefy_text = sprintf("%d_per_%d", numerator, denominator)
      )
    } else {
      tibble::tibble(
        number = rarefy_number(),
        rarefy_text = sprintf("%d_reads", number)
      )
    }) |>
      dplyr::bind_rows(
        tibble::tibble(rarefy_text = "full")
      ) |>
      dplyr::rename_with(\(x) if (isTRUE(dots)) paste0(".", x) else x)
  } else {
    tibble::tibble()
  }
}

#### main options function ####
#' @rdname parse_pipeline_options
#' @export
parse_pipeline_options <- function() {
  if (file.exists("pipeline_options.yaml")) {
    pipeline_options <- yaml::read_yaml("pipeline_options.yaml")
  } else {
    warning(
      "Options file 'pipeline_options.yaml' is missing!\n",
      "Using defaults for all parameters."
    )
    pipeline_options <- list()
  }

  parse_project_name(pipeline_options)
  parse_file_extension(pipeline_options)
  parse_orient(pipeline_options)
  parse_duplicate_policy(pipeline_options)
  parse_custom_sample_table(pipeline_options)
  parse_supplemental_asv_options(pipeline_options)
  parse_executable_options(pipeline_options)
  parse_top_dist_config(pipeline_options)
  parse_parallel_options(pipeline_options)
  parse_forward_primer(pipeline_options)
  parse_reverse_primer(pipeline_options)
  parse_trim_options(pipeline_options)
  parse_denoising_options(pipeline_options)
  parse_merge_options(pipeline_options)
  parse_filter_options(pipeline_options)
  parse_uncross_options(pipeline_options)
  parse_amplicon_model_options(pipeline_options)
  parse_lulu_options(pipeline_options)
  parse_output_options(pipeline_options)
  parse_control_options(pipeline_options)
  if (!is.null(pipeline_options$protax)) {
    parse_protax_options(pipeline_options$protax)
  }
  parse_taxonomy_options(pipeline_options)
  parse_added_reference(pipeline_options)
  parse_outgroup_options(pipeline_options)
  parse_cluster_options(pipeline_options)
  parse_guilds_options(pipeline_options)
  parse_rarefy_options(pipeline_options)
  options(optimotu.pipeline.did_options = TRUE)
}

#' @rdname pipeline_options
#' @export
did_pipeline_options <- function() {
  getOption("optimotu.pipeline.did_options", FALSE)
}
