#' Flatten a nested list of length-1 lists into a single list
#'
#' This function is useful for unnesting lists originally imported from YAML files
#' formatted like:
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
    all(vapply(x, length, 1L) == 1)
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
  if (!("project_name" %in% names(pipeline_options))
      || length(pipeline_options$project_name) == 0) {
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
    stop("Project name should consist of alphanumeric characters, '_', and",
         " '-'. (file:pipeline_options.yaml)")
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
    options(optimotu.pipeline.read_file_extension = pipeline_options$file_extension)
  }
}

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
    message("No read orientation specified in 'pipeline_options.yaml'.\n",
            "Using default: 'forward'")
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

#### added_reference ####
#' @rdname parse_pipeline_options
#' @keywords internal
parse_added_reference <- function(pipeline_options) {
  if (!is.null(pipeline_options$added_reference)) {
    checkmate::assert_list(pipeline_options$added_reference)
    added_reference <-
      unnest_yaml_list(pipeline_options$added_reference)
    checkmate::assert_string(added_reference$fasta, null.ok = TRUE)
    checkmate::assert_string(added_reference$table, null.ok = TRUE)

    if (xor(is.null(added_reference$fasta),
            is.null(added_reference$table))) {
      stop(
        "If one of 'added_reference_fasta' and 'added_reference_table' is given ",
        "in 'pipeline_options.yaml', then both must be given."
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

#### parallelism ####
#' @rdname parse_pipeline_options
#' @keywords internal
parse_parallel_options <- function(pipeline_options) {
  checkmate::assert_count(pipeline_options$local_threads, positive = TRUE, null.ok = TRUE)
  if (!is.null(pipeline_options$local_threads)) {
    options(optimotu_num_threads = pipeline_options$local_threads)
  }

  checkmate::assert_count(pipeline_options$max_batchsize, na.ok = TRUE, null.ok = TRUE)
  max_batchsize <- NULL
  if (checkmate::test_count(pipeline_options$max_batchsize, positive = TRUE))
    max_batchsize <- pipeline_options$max_batchsize

  workers_per_seqrun <- 1L
  checkmate::assert_count(pipeline_options$workers_per_seqrun, positive = TRUE, null.ok = TRUE)
  checkmate::assert_count(pipeline_options$jobs_per_seqrun, positive = TRUE, null.ok = TRUE)
  if (!is.null(pipeline_options$workers_per_seqrun)) {
    if (!is.null(pipeline_options$jobs_per_seqrun)) {
      warning(
        "both 'workers_per_seqrun' and 'jobs_per_seqrun' (deprecated) were given",
        " in 'pipeline_options.yaml'. Using 'workers_per_seqrun' (=",
        pipeline_options$workers_per_seqrun, ")")
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
  checkmate::assert_int(pipeline_options$min_workers, lower = 1L, null.ok = TRUE)
  if (!is.null(pipeline_options$min_workers))
    min_workers <- pipeline_options$min_workers

  max_workers <- Inf
  checkmate::assert_int(pipeline_options$max_workers, lower = min_workers, null.ok = TRUE)
  if (!is.null(pipeline_options$max_workers))
    max_workers <- pipeline_options$max_workers

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
    message("Forward primer string missing (file: pipeline_options.yaml)\n",
            "Using default: GCATCGATGAAGAACGCAGC")
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
    message("Reverse primer string missing (file: pipeline_options.yaml)\n",
            "Using default: TCCTCCGCTTATTGATATGC")
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
    message("No 'trimming' options given in 'pipeline_options.yaml'\n",
            "Using defaults.")
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

#### filtering settings ####
#' @rdname parse_pipeline_options
#' @keywords internal
parse_filter_options <- function(pipeline_options) {
  checkmate::assert_list(pipeline_options$filtering, null.ok = TRUE)
  if (is.null(pipeline_options$filtering)) {
    message("No 'filtering' options given in 'pipeline_options.yaml'\n",
            "Using defaults.")
  } else {
    filtering <- unnest_yaml_list(pipeline_options$filtering)
    checkmate::assert_names(
      names(filtering),
      subset.of = c("maxEE_R1", "maxEE_R2")
    )
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
      optimotu.pipeline.dada2_maxEE = update(dada2_maxEE(), filtering)
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
  if (is.null(lulu_options) || isFALSE(lulu_options)) return()
  checkmate::assert_list(lulu_options)
  lulu_options <- unnest_yaml_list(lulu_options)
  checkmate::assert_names(
    names(lulu_options),
    subset.of = c("dist_type", "max_dist", "max_gap_length", "max_gap_total",
                  "min_abundance_ratio", "min_cooccurrence_ratio",
                  "use_mean_abundance_ratio", "dist_config")
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
      fraction = checkmate::assert_number(lulu_options$max_dist,
                                          lower = 0, upper = 1),
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
      optimotu.pipeline.lulu_min_abundance_ratio =
        lulu_options$min_abundance_ratio
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
      optimotu.pipeline.lulu_min_cooccurrence_ratio =
        lulu_options$min_cooccurrence_ratio
    )
  }
  ##### use_mean_abundance_ratio #####
  if ("use_mean_abundance_ratio" %in% names(lulu_options)) {
    checkmate::assert_flag(lulu_options$use_mean_abundance_ratio)
    options(
      optimotu.pipeline.lulu_use_mean_abundance_ratio =
        lulu_options$use_mean_abundance_ratio
    )
  }
  ##### dist_config #####
  lulu_dist_config <- parse_dist_config(lulu_options$dist_config)
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
  getOption("optimotu.pipeline.lulu_dist_config", dist_config())
}


#### amplicon model settings ####

#' @rdname parse_pipeline_options
#' @export
#' @param amplicon_model_options (`list`) the amplicon model options to parse
parse_amplicon_model_options <- function(pipeline_options) {
  amplicon_model_options <- pipeline_options$amplicon_model
  if (is.null(amplicon_model_options) || isFALSE(amplicon_model_options)) return()

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
  options(optimotu.pipeline.amplicon_model_type = amplicon_model_options$model_type)

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
    options(optimotu.pipeline.amplicon_model_file = amplicon_model_options$model_file)

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
        if (!identical(amplicon_model_type(), "HMM"))
          stop("NuMt filter is only valid when HMM alignment is used")
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
#' @title Functions to access pipeline-wide options
#'
#' @description When used in a targets plan, these should always be
#' pre-evaluated with `!!` or `!!!` to ensure proper dependency tracking.
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
        stop("Option 'control':'spike' should be a file path, evaluate to ",
             "FALSE, or be left blank")
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
        stop("Option 'control':'positive' should be a file path, evaluate ",
             "to FALSE, or be left blank")
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
    stop("Only one of options 'taxonomy:protax', 'taxonomy:sintax',",
         " 'taxonomy:bayesant' and 'taxonomy:epa' may be given.",
         "(File: pipeline_options.yaml)")
  }
  if (length(selected_classifier) == 0) {
    stop("No classifier selected. Please select one of the following classifiers:\n",
         "  - protax\n",
         "  - sintax\n",
         "  - bayesant\n",
         "  - epa\n",
         "(File: pipeline_options.yaml)")
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
#' @param rank_options (`list` or `character` vector) the taxonomic ranks to use in the
#' pipeline, in order from most inclusive (e.g., kingdom) to least inclusive
#' (e.g., species). Values may be either named or unnamed. When named, the name is
#' taken to be the rank, and the value is the "in-group" taxon at that rank, i.e.
#' the taxon for which results are desired. When unnamed, the value is taken to
#' be the rank. Example: `list(kingdom = "Fungi", "phylum", "class", "order", "family", "genus", "species")`
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
      "Option 'taxonomy':'ranks' should start from the most inclusive rank (e.g. kingdom)\n",
      "  and continue to the least inclusive rank (e.g. species).  Optionally the first\n",
      "  rank(s) may be defined (e.g. '- kingdom: Fungi') but subsequent ranks must be \n",
      "  undefined (e.g. '- class')."
    )
  }
  if (length(KNOWN_TAXA) == 0) {
    KNOWN_TAXA = c(rootrank = "root")
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

  ##### protax version #####
  if ("aligned" %in% names(protax_options)) {
    checkmate::assert_flag(protax_options$aligned)
    if (protax_options$aligned && !do_model_align()) {
      stop("Aligned Protax (taxonomy: protax: aligned: true) requires model",
           "alignment to be enabled (amplicon_model: model_align: true).\n",
           "(file: pipeline_options.yaml)")
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
  checkmate::assert_file_exists(bayesant_options$reftax, "r")
  options(
    optimotu.pipeline.do_bayesant = TRUE,
    optimotu.pipeline.bayesant_ref = bayesant_options$reftax
  )
}

#' @rdname pipeline_options
#' @export
do_bayesant <- function() {
  getOption("optimotu.pipeline.do_bayesant", FALSE)
}

#' @rdname pipeline_options
#' @export
bayesant_ref <- function() {
  getOption("optimotu.pipeline.bayesant_ref")
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
  if (is.null(pipeline_options$outgroup_reference)) return()
  # if outgroup_reference is explicitly FALSE, then we set the option.
  # this may lead to cases where there is a defined ingroup, but no outgroup
  # reference.
  if (isFALSE(pipeline_options$outgroup_reference)) {
    if (do_outgroup()) {
      warning(
        "Taxonomic rank definitions imply the possibility of outgroup sequences",
        "but no outgroup reference file was provided.\n",
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
      message("No clustering options given in 'pipeline_options.yaml'\n",
              "Using defaults.")
    }
  } else {
    clustering <- unnest_yaml_list(pipeline_options$clustering)
    checkmate::assert_names(
      names(clustering),
      subset.of = c("thresholds", "measure", "dist_config")
    )
    checkmate::assert_file_exists(clustering$thresholds, "r")
    checkmate::assert_string(clustering$measure, null.ok = TRUE)
    if (!is.null(clustering$measure)) {
      checkmate::assert_subset(
        clustering$measure,
        c("MCC", "RI", "ARI", "FMI", "MI", "AMI", "FM")
      )
    }
    dist_config <- parse_dist_config(clustering$dist_config)
    options(
      optimotu.pipeline.clustering_thresholds = clustering$thresholds,
      optimotu.pipeline.clustering_measure = clustering$measure,
      optimotu.pipeline.clustering_dist_config = dist_config
    )
  }
}

parse_dist_config <- function(dist_config) {
    if (is.character(dist_config)) {
      dist_config <- list(method = dist_config)
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
    eval(dist_config)
}

#' @rdname pipeline_options
#' @export
cluster_thresholds <- function() {
  getOption("optimotu.pipeline.clustering_thresholds",
            "metadata/GSSP_thresholds.tsv")
}

#' @rdname pipeline_options
#' @export
cluster_measure <- function() {
  getOption("optimotu.pipeline.clustering_measure", NULL)
}

#' @rdname pipeline_options
#' @export
cluster_dist_config <- function() {
  getOption("optimotu.pipeline.clustering_dist_config", optimotu::dist_usearch())
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

#### OTU table settings ####
#' @rdname parse_pipeline_options
#' @keywords internal
parse_otu_table_options <- function(pipeline_options) {
  if (!is.null(pipeline_options$wide_table) && !is.null(pipeline_options$dense_table)) {
    stop(
      "Only one of 'wide_table' and 'dense_table' may be set in 'pipeline_options.yaml'.\n",
      "Please remove one of them."
    )
  }
  checkmate::assert_flag(pipeline_options$wide_table, null.ok = TRUE)
  checkmate::assert_flag(pipeline_options$dense_table, null.ok = TRUE)
  if (all(pipeline_options$wide_table, pipeline_options$dense_table)) {
    options(
      optimotu.pipeline.wide_table = TRUE
    )
  }
}

#' @rdname pipeline_options
#' @export
do_wide_otu_table <- function() {
  getOption("optimotu.pipeline.wide_table", FALSE)
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
      warning("Empty rarefy options given in 'pipeline_options.yaml'\n",
              "No rarefaction will be performed.")
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
#' @return a `data.frame` giving the variables to map over for rarefaction;
#' when rarefaction is to be performed (as determined by `do_rarefy()`) then
#' these will always include `.rarefy_text` (`character)`), and will also
#' include either `.numerator` (`integer`) and `.denominator` (`integer`) or
#' `.number` (`integer`).  Alternatively, if `do_rarefy()` returns `FALSE`, an
#' empty data.frame.
#' @export
rarefy_meta <- function(dots = TRUE) {
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
  parse_added_reference(pipeline_options) #TODO: this should go to taxonomy?
  parse_parallel_options(pipeline_options)
  parse_forward_primer(pipeline_options)
  parse_reverse_primer(pipeline_options)
  parse_trim_options(pipeline_options)
  parse_filter_options(pipeline_options)
  parse_uncross_options(pipeline_options)
  parse_amplicon_model_options(pipeline_options)
  parse_lulu_options(pipeline_options)
  parse_otu_table_options(pipeline_options)
  parse_control_options(pipeline_options)
  if (!is.null(pipeline_options$protax)) {
    parse_protax_options(pipeline_options$protax)
  }
  parse_taxonomy_options(pipeline_options)
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
