#' Try to find the cutadapt executable
#' @return a `character` string giving the path of the cutadapt executable
#' @export
find_cutadapt <- function() {
  find_executable("cutadapt")
}

#' Options for cutadapt filtering and trimming
#' @rdname cutadapt_paired_options
#' @export
cutadapt_paired_option_names <- c(
  "max_err",
  "min_overlap",
  "action",
  "discard_untrimmed",
  "max_n",
  "max_ee",
  "min_length",
  "max_length",
  "truncQ_R1",
  "truncQ_R2",
  "cut_R1",
  "cut_R2"
)

#' Options for cutadapt paired-end filtering and trimming
#' @param max_err (`numeric`) maximum error rate for adapter matching
#' @param min_overlap (`integer`) minimum overlap between read and adapter
#' @param action (`character`) action to take on the read; options are "trim",
#' "retain", "mask", "lowercase", "none"
#' @param discard_untrimmed (`logical`) discard reads that do not match the
#' adapter
#' @param max_n (`integer`) maximum number of Ns allowed in a read
#' @param max_ee (`numeric`) maximum expected errors in a read after trimming
#' @param min_length (`integer`) minimum length of a read after trimming
#' @param max_length (`integer`) maximum length of a read after trimming
#' @param truncQ_R1 (`integer`) quality score threshold for truncating the 3'
#' end of read 1; or a vector of two integers for the 5' and 3' ends
#' @param truncQ_R2 (`integer`) quality score threshold for truncating the 3'
#' end of read 2; or a vector of two integers for the 5' and 3' ends
#' @param cut_R1 (`integer`) number of bases to unconditionally cut from the 5'
#' (if positive) or 3' (if negative) end of read 1 prior to adapter trimming;
#' or a vector of two integers with different signs for both ends
#' @param cut_R2 (`integer`) number of bases to unconditinoally cut from the 5'
#' (if positive) or 3' (if negative) end of read 2 prior to adapter trimming;
#' or a vector of two integers with different signs for both ends
#' @return an object of class `cutadapt_paired_options`
#' @export
cutadapt_paired_options <- function(
    max_err = 0.2,
    min_overlap = 10L,
    action = "retain",
    discard_untrimmed = TRUE,
    max_n = 0L,
    max_ee = NULL,
    min_length = NULL, max_length = NULL,
    truncQ_R1 = 2, truncQ_R2 = 2,
    cut_R1 = NULL, cut_R2 = NULL
) {
  checkmate::assert_number(max_err, lower = 0, null.ok = TRUE)
  checkmate::assert_count(min_overlap, positive = TRUE, null.ok = TRUE)
  checkmate::assert_choice(
    action,
    c("trim", "retain", "mask", "lowercase", "none"),
    null.ok = TRUE
  )
  checkmate::assert_flag(discard_untrimmed)
  checkmate::assert_count(max_n, null.ok = TRUE)
  checkmate::assert_number(max_ee, lower = 0, finite = TRUE, null.ok = TRUE)
  checkmate::assert_count(min_length, positive = TRUE, null.ok = TRUE)
  checkmate::assert_integerish(truncQ_R1, min.len = 1, max.len = 2, null.ok = TRUE)
  checkmate::assert_integerish(truncQ_R2, min.len = 1, max.len = 2, null.ok = TRUE)
  checkmate::assert_integerish(cut_R1, min.len = 1, max.len = 2, null.ok = TRUE)
  checkmate::assert_integerish(cut_R2, min.len = 1, max.len = 2, null.ok = TRUE)
  structure(
    list(
      max_err = max_err,
      min_overlap = min_overlap,
      action = action,
      discard_untrimmed = discard_untrimmed,
      max_n = max_n,
      max_ee = max_ee,
      min_length = min_length,
      truncQ_R1 = truncQ_R1,
      truncQ_R2 = truncQ_R2,
      cut_R1 = cut_R1,
      cut_R2 = cut_R2
    ),
    class = "cutadapt_paired_options"
  )
}

#' Run cutadapt on paired-end reads
#' @param file_R1 (`character`) path to the R1 input file
#' @param file_R2 (`character`) path to the R2 input file
#' @param primer_R1 (`character`) primer sequence for R1
#' @param primer_R2 (`character`) primer sequence for R2
#' @param trim_R1 (`character`) path to the trimmed R1 output file
#' @param trim_R2 (`character`) path to the trimmed R2 output file
#' @param options ([`cutadapt_paired_options`](`cutadapt_paired_options()`))
#' options for cutadapt
#' @param ncpu (`integer`) number of CPUs to use
#' @param cutadapt (`character`) path to the cutadapt executable
#' @param logfile (`character`) path to a log file which will receive the
#' cutadapt stdout.
#' @param ... additional arguments; currently ignored
#' @return a character vector of the trimmed output files
#' @export
cutadapt_paired_filter_trim <- function(
  file_R1, file_R2,
  primer_R1, primer_R2,
  trim_R1, trim_R2,
  options = cutadapt_paired_options(),
  ncpu = local_cpus(),
  cutadapt = find_cutadapt(),
  logfile = NULL,
  ...
) {
  checkmate::assert_class(options, "cutadapt_paired_options")
  args <- c(
    "-g", primer_R1,
    "-G", primer_R2,
    "-o", trim_R1,
    "-p", trim_R2
  )
  if (!is.null(options$max_err)) {
    args <- c(args, "-e", options$max_err)
  }
  if (!is.null(options$min_overlap)) {
    args <- c(args, "-O", options$min_overlap)
  }
  if (!is.null(options$action)) {
    args <- c(args, paste0("--action=", options$action))
  }
  if (isTRUE(options$discard_untrimmed)) {
    args <- c(args, "--discard-untrimmed")
  }
  if (!is.null(options$max_n)) {
    args <- c(args, "--max-n", options$max_n)
  }
  if (!is.null(options$max_ee)) {
    args <- c(args, "--max-ee", options$max_ee)
  }
  if (!is.null(options$min_length)) {
    args <- c(args, "-m", options$min_length)
  }
  if (!is.null(options$max_length)) {
    args <- c(args, "-M", options$max_length)
  }
  if (!is.null(options$truncQ_R1)) {
    args <- c(args, "-q", paste(round(options$truncQ_R1), collapse = ","))
  }
  if (!is.null(options$truncQ_R2)) {
    args <- c(args, "-Q", paste(round(options$truncQ_R2), collapse = ","))
  }
  if (!is.null(options$cut_R1)) {
    args <- c(args, t(data.frame("-u", round(options$cut_R1))))
  }
  if (!is.null(options$cut_R2)) {
    args <- c(args, t(data.frame("-U", round(options$cut_R2))))
  }
  if (!is.null(ncpu)) {
    checkmate::assert_count(ncpu, positive = TRUE)
    args <- c(args, "-j", ncpu)
  }

  checkmate::assert_string(file_R1)
  checkmate::assert_string(file_R2)

  file_R1 <- strsplit(file_R1, ",")[[1]]
  file_R2 <- strsplit(file_R2, ",")[[1]]
  if (length(file_R1) > 1 || length(file_R2) > 1) {
    seqs_R1 <- Biostrings::readQualityScaledDNAStringSet(file_R1)
    file_R1 <- withr::local_tempfile(fileext = ".fastq.gz")
    Biostrings::writeQualityScaledXStringSet(seqs_R1, file_R1, compress = TRUE)
    seqs_R2 <- Biostrings::readQualityScaledDNAStringSet(file_R2)
    file_R2 <- withr::local_tempfile(fileext = ".fastq.gz")
    Biostrings::writeQualityScaledXStringSet(seqs_R2, file_R2, compress = TRUE)
  }

  args <- c(args, file_R1, file_R2)
  out <- processx::run(
    cutadapt,
    args = args,
    error_on_status = TRUE
  )
  if (!is.null(logfile)) writeLines(out$stdout, logfile)
  c(trim_R1, trim_R2)
}

#' @rdname cutadapt_options
#' @export
cutadapt_option_names <- c(
  "max_err",
  "min_overlap",
  "action",
  "discard_untrimmed",
  "max_n",
  "max_ee",
  "min_length",
  "max_length",
  "truncQ",
  "cut"
)

#' Options for cutadapt single-read filtering and trimming
#' @param max_err (`numeric`) maximum error rate for adapter matching
#' @param min_overlap (`integer`) minimum overlap between read and adapter
#' @param action (`character`) action to take on the read; options are "trim",
#' "retain", "mask", "lowercase", "none"
#' @param discard_untrimmed (`logical`) discard reads that do not match the
#' adapter
#' @param max_n (`integer`) maximum number of Ns allowed in a read after trimming
#' @param max_ee (`numeric`) maximum expected errors in a read after trimming
#' @param min_length (`integer`) minimum length of a read after trimming
#' @param max_length (`integer`) maximum length of a read after trimming
#' @param truncQ (`integer`) quality score threshold for truncating the 3' end
#' of the read; or a vector of two integers for the 5' and 3' ends
#' @param cut (`integer`) number of bases to unconditionally cut from the 5' (if
#' positive) or 3' (if negative) end of the read prior to adapter trimming; or
#' a vector of two integers with different signs for both ends
#' @return an object of class `cutadapt_options`
#' @export
cutadapt_options <- function(
    max_err = 0.2,
    min_overlap = 10L,
    action = "retain",
    discard_untrimmed = TRUE,
    max_n = 0L,
    max_ee = NULL,
    min_length = NULL, max_length = NULL,
    truncQ = NULL,
    cut = NULL
) {
  checkmate::assert_number(max_err, lower = 0, null.ok = TRUE)
  checkmate::assert_count(min_overlap, positive = TRUE, null.ok = TRUE)
  checkmate::assert_choice(
    action,
    c("trim", "retain", "mask", "lowercase", "none"),
    null.ok = TRUE
  )
  checkmate::assert_flag(discard_untrimmed)
  checkmate::assert_count(max_n, null.ok = TRUE)
  checkmate::assert_number(max_ee, lower = 0, finite = TRUE, null.ok = TRUE)
  checkmate::assert_count(min_length, positive = TRUE, null.ok = TRUE)
  checkmate::assert_count(max_length, positive = TRUE, null.ok = TRUE)
  checkmate::assert_integerish(truncQ, min.len = 1, max.len = 2, null.ok = TRUE)
  checkmate::assert_integerish(cut, min.len = 1, max.len = 2, null.ok = TRUE)
  structure(
    list(
      max_err = max_err,
      min_overlap = min_overlap,
      action = action,
      discard_untrimmed = discard_untrimmed,
      max_n = max_n,
      max_ee = max_ee,
      min_length = min_length,
      max_length = max_length,
      truncQ = truncQ,
      cut = cut
    ),
    class = "cutadapt_options"
  )
}

#' Replace existing cutadapt options with new values
#'
#' @param object ([`cutadapt_paired_options`](`cutadapt_paired_options()`) or
#' [`cutadapt_options`](`cutadapt_options()`)) existing cutadapt options to
#' update
#' @param new_options (`data.frame`, named `list`, or named `character` vector)
#' new values for the options. If a `data.frame`, then the values in the
#' `data.frame` should be all the same.
#' @param ... additional arguments; currently ignored
#' @return updated options
#' @exportS3Method stats::update
#' @rdname update_cutadapt_options
update.cutadapt_paired_options <- function(object, new_options, ...) {
  checkmate::assert(
    checkmate::check_list(new_options, null.ok = TRUE),
    checkmate::check_data_frame(new_options, null.ok = TRUE),
    checkmate::check_character(new_options, null.ok = TRUE)
  )
  new_options <-
    new_options[intersect(names(new_options), cutadapt_paired_option_names)]
  if (is.data.frame(new_options) && ncol(new_options) > 0) {
    new_options <- unique(new_options)
    if (nrow(new_options) > 1L)
      stop(
        "'new_options' must be the same for all samples in a batch. \n",
        "Current batch has ", nrow(new_options),
        "unique combinations of options."
      )
  }
  new_options <- lapply(unclass(new_options), unlist)
  if (length(new_options) > 0) {
    object[names(new_options)] <- new_options
    do.call(cutadapt_paired_options, object)
  } else {
    object
  }
}

#' @rdname update_cutadapt_options
#' @exportS3Method stats::update
update.cutadapt_options <- function(object, new_options, ...) {
  checkmate::assert(
    checkmate::check_list(new_options, null.ok = TRUE),
    checkmate::check_data_frame(new_options, null.ok = TRUE),
    checkmate::check_character(new_options, null.ok = TRUE)
  )
  new_options <-
    new_options[intersect(names(new_options), cutadapt_option_names)]
  if (is.data.frame(new_options) && ncol(new_options) > 0) {
    new_options <- unique(new_options)
    if (nrow(new_options) > 1L)
      stop(
        "'new_options' must be the same for all samples in a batch. \n",
        "Current batch has ", nrow(new_options),
        "unique combinations of options."
      )
  }
  new_options <- lapply(unclass(new_options), unlist)
  if (length(new_options) > 0) {
    object[names(new_options)] <- new_options
    do.call(cutadapt_options, object)
  } else {
    object
  }
}

#' Run cutadapt on single-end reads
#' @param file (`character`) path to the input file
#' @param primer (`character`) primer sequence
#' @param trim (`character`) path to the trimmed output file
#' @param options ([`cutadapt_options`](`cutadapt_options()`)) options for
#' cutadapt
#' @param ncpu (`integer`) number of CPUs to use
#' @param cutadapt (`character`) path to the cutadapt executable
#' @param logfile (`character`) path to a log file which will receive the
#' cutadapt stdout.
#' @param ... additional arguments; currently ignored
#' @return a `character` string giving the trimmed output file
#' @export
cutadapt_filter_trim <- function(
  file,
  primer,
  trim,
  options = cutadapt_options(),
  ncpu = local_cpus(),
  cutadapt = find_cutadapt(),
  logfile = NULL,
  ...
) {
  args <- c(
    "-g", primer,
    "-o", trim
  )
  if (!is.null(options$max_err)) {
    args <- c(args, "-e", options$max_err)
  }
  if (!is.null(options$min_overlap)) {
    args <- c(args, "-O", round(options$min_overlap))
  }
  if (!is.null(options$action)) {
    args <- c(args, paste0("--action=", options$action))
  }
  if (isTRUE(options$discard_untrimmed)) {
    args <- c(args, "--discard-untrimmed")
  }
  if (!is.null(options$max_n)) {
    args <- c(args, "--max-n", options$max_n)
  }
  if (!is.null(options$max_ee)) {
    args <- c(args, "--max-ee", options$max_ee)
  }
  if (!is.null(options$min_length)) {
    args <- c(args, "-m", options$min_length)
  }
  if (!is.null(options$max_length)) {
    args <- c(args, "-M", options$max_length)
  }
  if (!is.null(options$truncQ)) {
    args <- c(args, "-q", paste(round(options$truncQ), collapse = ","))
  }
  if (!is.null(options$cut)) {
    args <- c(args, t(data.frame("-u", round(options$cut))))
  }
  if (!is.null(ncpu)) {
    checkmate::assert_count(ncpu, positive = TRUE)
    args <- c(args, "-j", ncpu)
  }
  args <- c(args, file)
  out <- processx::run(
    cutadapt,
    args = args,
    error_on_status = TRUE
  )
  if (!is.null(logfile)) writeLines(out$stdout, logfile)
  trim
}

#' Trim primers from a set of sequences
#' @param seqs ([`DNAStringSet`][Biostrings::XStringSet-class], `character`
#' vector, or `data.frame`) sequences to trim
#' @param primer (`character`) primer sequence
#' @param ... additional arguments to `cutadapt_filter_trim()`
#' @return a `tibble` with columns `seq_id` and `seq`
#' @export
trim_primer <- function(seqs, primer, ...) {
  tempseqs <- withr::local_tempfile(fileext = ".fasta")
  write_sequence(seqs, tempseqs)
  temptrimmed <- withr::local_tempfile(fileext = ".fasta")
  cutadapt_filter_trim(
    file = tempseqs,
    primer = primer,
    trim = temptrimmed,
    ...
  )
  Biostrings::readDNAStringSet(temptrimmed) |>
    as.character() |>
    tibble::enframe(name = "seq_id", value = "seq")
}


#' Top-level function to trim a batch of raw paired-end reads
#' @param pairs_meta (`data.frame`) metadata for the paired-end reads
#' @param seqrun (`character`) seqrun identifier for the batch
#' @param orient (`character`) orientation of the reads for the batch; either
#' "fwd" or "rev"
#' @param trim_options ([`cutadapt_paired_options`](`cutadapt_paired_options()`)
#' options for cutadapt
#' @param primer_R1 (`character`) primer specification for R1
#' @param primer_R2 (`character`) primer specification for R2
#' @param ncpu (`integer`) number of CPU threads to use
#' @param ... additional arguments for dependency tracking
#' @return a `character` vector with the path to the trimmed output files
#' @export
trim_raw_pairs <- function(
    pairs_meta,
    seqrun,
    orient,
    trim_options,
    primer_R1,
    primer_R2,
    ncpu = local_cpus(),
    ...
) {
  logfile_name <- sprintf("logs/trim_%s_%s.log", seqrun, orient)
  logfile <- withr::local_connection(file(logfile_name, "w"))

  dplyr::group_by(
    pairs_meta,
    dplyr::pick(any_of(c("seqrun", cutadapt_paired_option_names)))
  ) |>
    dplyr::group_map(
      ~ purrr::pmap(
        dplyr::transmute(
          .x,
          file_R1 = fastq_R1,
          file_R2 = fastq_R2,
          trim_R1 = trim_R1,
          trim_R2 = trim_R2
        ),
        cutadapt_paired_filter_trim,
        primer_R1 = ifelse(orient == "fwd", primer_R1, primer_R2),
        primer_R2 = ifelse(orient == "fwd", primer_R2, primer_R1),
        options = update(trim_options, .y),
        ncpu = ncpu,
        logfile = logfile
      ) |>
        unlist(),
      .keep = TRUE
    ) |>
    unlist()
}
