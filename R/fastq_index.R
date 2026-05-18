#' @title Indexing and random access extraction of gzipped FASTA and FASTQ files
#'
#' @param file (`character` filename) file to create an index for
#'
#' These functions are wrappers around the
#' [`FastqIndEx`](https://dkfz-odcf.github.io/FastqIndEx/) utility.
#'
#'
#' @return file name of the created index
#' @describeIn fastx_gz Generate a gzip index file
#' @export
fastx_gz_index <- function(file) {
  fastqindex <- find_executable("fastqindex")
  checkmate::assert_file_exists(fastqindex, access = "x")
  index <- sprintf("%s.fqi", file)
  args <- c(
    "index",
    sprintf("-f=%s", file),
    sprintf("-i=%s", index),
    "-w"
  )
  out <- system2(fastqindex, args)
  stopifnot(out == 0)
  checkmate::assert_file_exists(index, "r")
  index
}

#' Resolve `seq_idx` for indexed extraction (`NULL` means all records, in order)
#'
#' @param index_obj A [`fastqindexr_index`][fastqindexr::create_index()] object.
#' @param seq_idx (`integer` or `NULL`) subset to extract; `integer()` requests
#'   an empty extraction.
#'
#' @return Integer vector of 1-based indices into the logically concatenated
#'   records.
#' @noRd
resolve_indexed_seq_idx <- function(index_obj, seq_idx) {
  checkmate::assert_class(index_obj, "fastqindexr_index")
  if (is.null(seq_idx)) {
    n <- as.integer(index_obj$n_records)
    if (n < 1L) {
      return(integer())
    }
    return(seq_len(n))
  }
  checkmate::assert_integerish(seq_idx, lower = 1L, any.missing = FALSE)
  as.integer(seq_idx)
}

#' Resolve `seq_idx` against a fixed record count (in-memory or flat file read)
#'
#' @param n_records Total number of sequences in order (`integer` scalar).
#' @param seq_idx (`NULL` or `integer`) same semantics as
#'   resolve_indexed_seq_idx semantics but bounds-checked against `n_records`.
#'
#' @return Integer vector of 1-based indices (possibly empty).
#' @noRd
resolve_linear_seq_idx <- function(n_records, seq_idx) {
  checkmate::assert_count(n_records)
  if (is.null(seq_idx)) {
    if (n_records < 1L) {
      return(integer())
    }
    return(seq_len(n_records))
  }
  checkmate::assert_integerish(
    seq_idx,
    lower = 1L,
    upper = n_records,
    any.missing = FALSE
  )
  as.integer(seq_idx)
}

#' Split a vector into up to `ncpu` contiguous blocks of nearly equal size
#'
#' Used to parallelize tools over sequence batches. Empty `vec` yields a
#' single empty block.
#'
#' @param vec Atomic vector (often integer indices).
#' @param ncpu Maximum number of blocks (`integer` scalar, at least 1).
#'
#' @return `list` of pieces of `vec`, length `min(ncpu, length(vec))` when
#'   `length(vec) > 0`, else a one-element list holding an empty integer vector.
#' @noRd
partition_vector_equal_ncpu <- function(vec, ncpu) {
  checkmate::assert_count(ncpu)
  m <- length(vec)
  if (m == 0L) {
    return(list(integer()))
  }
  k <- min(as.integer(ncpu), m)
  if (k < 2L) {
    return(list(vec))
  }
  base <- m %/% k
  rem <- m %% k
  group_sizes <- rep.int(base, k)
  if (rem > 0L) {
    group_sizes[seq_len(rem)] <- group_sizes[seq_len(rem)] + 1L
  }
  ends <- cumsum(group_sizes)
  starts <- c(1L, ends[-k] + 1L)
  lapply(seq_len(k), function(i) vec[starts[i]:ends[i]])
}

normalize_fastx_extract_inputs <- function(infile, index) {
  checkmate::assert_character(infile, min.len = 1L, any.missing = FALSE)
  checkmate::assert_file_exists(infile, "r")
  if (inherits(index, "fastqindexr_index")) {
    if (length(infile) != length(index$files)) {
      stop(
        "`infile` must have the same length as the files in `index`.",
        call. = FALSE
      )
    }
    return(
      list(
        index = index,
        file = infile
      )
    )
  }
  checkmate::assert_character(index, min.len = 1L, any.missing = FALSE)
  checkmate::assert_file_exists(index, "r")
  if (length(index) != length(infile)) {
    stop("`index` and `infile` must have the same length.", call. = FALSE)
  }
  list(
    index = fastqindexr::read_fqi_index(
      fqi_path = index,
      files = infile,
      type = "auto"
    ),
    file = NULL
  )
}

#' @param infile (`character` filename) gzipped fasta or fastq file
#' @param index (`character` filename(s) or `fastqindexr_index` object) index
#'   for `infile`
#' @param i (`integer` vector) indices to extract
#' @param outfile (`character` filename) file to write the extracted sequences
#'   to. If it ends in ".gz", the output will be gzipped.
#' @param renumber (`logical` flag) if `TRUE`, replace the sequence names with
#'   integers, starting at 0.
#' @param append (`logical` flag) if `TRUE`, append to `outfile` if it already
#'   exists, rather than overwriting.
#' @param hash (`character` scalar) md5 hash of the infile; ignored but included
#'   as a parameter for dependency tracking
#'
#' @return filename of the output file
#' @describeIn fastx_gz Extract sequences from a gzipped FASTA or FASTQ file
#' @export
fastx_gz_extract <- function(
  infile,
  index,
  i,
  outfile,
  renumber = FALSE,
  append = FALSE,
  hash = NULL
) {
  checkmate::assert_integerish(i, lower = 1)
  checkmate::assert_string(outfile)
  checkmate::assert_flag(renumber)
  checkmate::assert_flag(append)
  input <- normalize_fastx_extract_inputs(infile = infile, index = index)
  if (file.exists(outfile) && !append) {
    unlink(outfile)
  }
  ensure_directory(outfile)
  if (length(i) < 1L && !file.exists(outfile)) {
    file.create(outfile)
    return(outfile)
  }
  fastqindexr::extract_sequences_to_file(
    index = input$index,
    seq_idx = i,
    file = input$file,
    outfile = outfile,
    type = "auto",
    append = append,
    compress = endsWith(outfile, ".gz"),
    renumber = if (isTRUE(renumber)) "zero_based" else "none"
  )
  outfile
}

#' @param infile (`character` filename) gzipped fasta or fastq file
#' @param index (`character` filename(s) or `fastqindexr_index` object) index
#'   for `infile`
#' @param i (`integer` vector) indices to extract
#' @param outfile (`character` filename) file to write the extracted sequences
#'   to. If it ends in ".gz", the output will be gzipped.
#' @param max_gap (`integer` scalar) maximum number of consecutive missing
#'   sequences to allow in an extraction batch. This parameter is for
#'   performance tuning and does not affect the results.
#' @param ncpu (`integer` scalar) number of threads to use.  Ignored in current
#'   implementation, but included for backwards compatibility.
#' @describeIn fastx_gz Extract sequences from a gzipped FASTA or FASTQ file
#' @export
fastx_gz_random_access_extract <- function(
  infile,
  index,
  i,
  outfile = NULL,
  renumber = FALSE,
  append = FALSE,
  hash = NULL,
  max_gap = 100L,
  ncpu = local_cpus()
) {
  checkmate::assert_integerish(i, lower = 1)
  checkmate::assert_string(outfile, null.ok = TRUE)
  checkmate::assert_flag(renumber)
  checkmate::assert_flag(append)
  checkmate::assert_integerish(max_gap, lower = 1)
  input <- normalize_fastx_extract_inputs(infile = infile, index = index)

  if (is.null(outfile)) {
    fastqindexr::extract_sequences_dnastringset(
      index = input$index,
      seq_idx = i,
      file = input$file,
      renumber = if (isTRUE(renumber)) "zero_based" else "none"
    )
  } else {
    if (!append && file.exists(outfile)) {
      unlink(outfile)
    }
    fastqindexr::extract_sequences_to_file(
      index = input$index,
      seq_idx = i,
      file = input$file,
      outfile = outfile,
      type = "fasta",
      append = append,
      compress = endsWith(outfile, ".gz"),
      renumber = if (isTRUE(renumber)) "zero_based" else "none"
    )
  }
}

#' Generate MD5 hash of a subset of sequences in a gzipped FASTA or FASTQ file
#' @inheritParams fastx_gz_extract
#' @param start (`integer` scalar) one-based index to start hashing
#' @param n (`integer` scalar) number of sequences to hash
#'
#' @return (`character`) md5 hash
#' @export
fastx_gz_hash <- function(infile, index, start, n) {
  checkmate::assert_file_exists(infile, "r")
  checkmate::assert(
    checkmate::check_file_exists(index, "r"),
    checkmate::check_class(index, "fastqindexr_index"),
    combine = "or"
  )
  checkmate::assert_integerish(start, lower = 1)
  checkmate::assert_integerish(n, lower = 1)

  tmp_fifo <- withr::local_tempfile(fileext = ".fifo")
  system2("mkfifo", tmp_fifo)
  md5_run <- processx::process$new(
    command = "md5sum",
    args = tmp_fifo,
    stdout = "|"
  )
  fastqindexr::extract_sequences_to_file(
    index = index,
    seq_idx = seq(start, start + n - 1L),
    file = infile,
    outfile = tmp_fifo,
    type = "auto",
    append = FALSE,
    compress = FALSE,
    collapse_sequence_lines = FALSE,
    renumber = "none"
  )
  md5_run$wait()
  stopifnot(md5_run$get_exit_status() == 0)
  c(strtrim(md5_run$read_all_output_lines(), 32))
}
