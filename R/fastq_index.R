#' @title Indexing and random access extraction of gzipped FASTA and FASTQ files
#'
#' @param file (`character` filename) file to create an index for
#' @param fastqindex (`character`) path to the `fastqindex` executable
#'
#' These functions are wrappers around the
#' [`FastqIndEx`](https://dkfz-odcf.github.io/FastqIndEx/) utility.
#'
#'
#' @return file name of the created index
#' @describeIn fastx_gz Generate a gzip index file
#' @export
fastx_gz_index <- function(file, fastqindex = find_executable("fastqindex")) {
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

#' Whether character paths are existing `.fqi` or `.qs2` index files
#'
#' @param x Candidate index path(s).
#' @return `TRUE` if `x` is a non-missing character vector of existing paths
#'   that all end in `.fqi` or all end in `.qs2` (case-insensitive).
#' @noRd
is_fastqindexr_index_path_set <- function(x) {
  if (!is.character(x) || anyNA(x) || length(x) < 1L) {
    return(FALSE)
  }
  if (!all(file.exists(x))) {
    return(FALSE)
  }
  is_fqi <- grepl("\\.fqi$", x, ignore.case = TRUE)
  is_qs2 <- grepl("\\.qs2$", x, ignore.case = TRUE)
  (all(is_fqi) && !any(is_qs2)) || (all(is_qs2) && !any(is_fqi))
}

#' Load a `fastqindexr_index` from `.fqi` and/or `.qs2` path(s)
#'
#' A single `.qs2` file holds one index object (possibly covering multiple
#' FASTA/FASTQ files). Multiple `.fqi` paths are passed to
#' [fastqindexr::read_fqi_index()]. Mixing `.fqi` and `.qs2`, or multiple
#' `.qs2` paths, is not supported.
#'
#' @param index Character path(s) to index file(s).
#' @param files Optional FASTA/FASTQ path override(s) for the index.
#' @return A `fastqindexr_index` object.
#' @noRd
load_fastqindexr_index <- function(index, files = NULL) {
  checkmate::assert_character(index, min.len = 1L, any.missing = FALSE)
  checkmate::assert_file_exists(index, "r")
  checkmate::assert_character(files, null.ok = TRUE)
  if (all(grepl("\\.qs2$", index, ignore.case = TRUE))) {
    if (length(index) != 1L) {
      stop(
        "Multiple .qs2 index paths are not supported; use one .qs2 file ",
        "or .fqi path(s).",
        call. = FALSE
      )
    }
    if (!requireNamespace("qs2", quietly = TRUE)) {
      stop(
        "Package 'qs2' is required to read .qs2 index files.",
        call. = FALSE
      )
    }
    idx <- qs2::qs_read(index)
    if (!inherits(idx, "fastqindexr_index")) {
      stop(
        "Expected a fastqindexr_index object in ",
        index,
        ".",
        call. = FALSE
      )
    }
    return(idx)
  }
  if (!all(grepl("\\.fqi$", index, ignore.case = TRUE))) {
    stop(
      "`index` paths must all end in .fqi or a single path ending in .qs2.",
      call. = FALSE
    )
  }
  fastqindexr::read_fqi_index(
    fqi_path = index,
    files = files,
    type = "auto"
  )
}

#' Write a `fastqindexr_index` to a qs2 file
#'
#' Drops the live `._cache` environment (native `index_ptr`) before
#' serializing. Reloaded indexes rebuild cache on first use.
#'
#' @param index A [`fastqindexr_index`][fastqindexr::create_index()] object.
#' @param file (`character` scalar) output path; typically ends in `.qs2`.
#'
#' @return `file`, invisibly via [write_and_return_file()].
#' @export
write_fastqindexr_index <- function(index, file) {
  checkmate::assert_class(index, "fastqindexr_index")
  checkmate::assert_string(file)
  index$`._cache` <- NULL
  write_and_return_file(index, file, type = "qs2")
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
  if (all(grepl("\\.qs2$", index, ignore.case = TRUE))) {
    idx <- load_fastqindexr_index(index)
    if (length(infile) != length(idx$files)) {
      stop(
        "`infile` must have the same length as the files in `index`.",
        call. = FALSE
      )
    }
    return(
      list(
        index = idx,
        file = infile
      )
    )
  }
  if (length(index) != length(infile)) {
    stop("`index` and `infile` must have the same length.", call. = FALSE)
  }
  list(
    index = load_fastqindexr_index(index, files = infile),
    file = NULL
  )
}

#' @param infile (`character` filename) gzipped fasta or fastq file
#' @param index (`character` filename(s) or `fastqindexr_index` object) index
#'   for `infile`. Character paths may be `.fqi` file(s) or a single `.qs2`
#'   file from [write_fastqindexr_index()].
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

  input <- normalize_fastx_extract_inputs(infile = infile, index = index)
  tmp_fifo <- withr::local_tempfile(fileext = ".fifo")
  system2("mkfifo", tmp_fifo)
  md5_run <- processx::process$new(
    command = "md5sum",
    args = tmp_fifo,
    stdout = "|"
  )
  fastqindexr::extract_sequences_to_file(
    index = input$index,
    seq_idx = seq(start, start + n - 1L),
    file = input$file,
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
