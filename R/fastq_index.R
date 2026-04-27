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

write_fastx_dnastringset <- function(seqs, outfile, append = FALSE) {
  ensure_directory(outfile)
  Biostrings::writeXStringSet(
    x = seqs,
    filepath = outfile,
    append = append,
    compress = endsWith(outfile, ".gz")
  )
  outfile
}

extract_fastx_as_dnastringset <- function(index, file, i, renumber = FALSE) {
  seqs <- fastqindexr::extract_sequences(
    index = index,
    seq_idx = i,
    file = file,
    return = "seq"
  )
  out <- Biostrings::DNAStringSet(unname(seqs))
  names(out) <- names(seqs)
  if (renumber) {
    names(out) <- as.character(seq_along(out))
  }
  out
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
  if (isTRUE(renumber)) {
    seqs <- extract_fastx_as_dnastringset(
      index = input$index,
      file = input$file,
      i = i,
      renumber = FALSE
    )
    names(seqs) <- as.character(seq_along(seqs) - 1L)
    return(write_fastx_dnastringset(seqs, outfile = outfile, append = append))
  }
  fastqindexr::extract_sequences_to_file(
    index = input$index,
    seq_idx = i,
    file = input$file,
    outfile = outfile,
    type = "auto",
    append = append,
    compress = endsWith(outfile, ".gz")
  )
  outfile
}

#' @param infile (`character` filename) gzipped fasta or fastq file
#' @param index (`character` filename(s) or `fastqindexr_index` object) index
#'   for `infile`
#' @param i (`integer` vector) indices to extract
#' @param outfile (`character` filename) file to write the extracted sequences
#'   to. If it ends in ".gz", the output will be gzipped.
#' @inheritParams fastx_gz_extract
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
  seqs <- extract_fastx_as_dnastringset(
    index = input$index,
    file = input$file,
    i = i,
    renumber = renumber
  )
  if (is.null(outfile)) {
    seqs
  } else {
    if (!append && file.exists(outfile)) {
      unlink(outfile)
    }
    if (isTRUE(renumber)) {
      return(write_fastx_dnastringset(
        seqs = seqs,
        outfile = outfile,
        append = append
      ))
    }
    fastqindexr::extract_sequences_to_file(
      index = input$index,
      seq_idx = i,
      file = input$file,
      outfile = outfile,
      type = "fasta",
      append = append,
      compress = endsWith(outfile, ".gz")
    )
    outfile
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
  checkmate::assert_file_exists(index, "r")
  checkmate::assert_integerish(start, lower = 1)
  checkmate::assert_integerish(n, lower = 1)
  fastqindex <- find_executable("fastqindex")
  checkmate::assert_file_exists(fastqindex, access = "x")
  is_fastq <- endsWith(infile, "fastq.gz") || endsWith(infile, "fq.gz")
  command <- sprintf(
    "%s extract -s=%i -n=%i -e=%i -f=%s -i=%s | tail -n+8 | md5sum",
    fastqindex,
    start - 1L,
    n,
    if (is_fastq) 4 else 2,
    infile,
    index
  )
  result <- system(command, intern = TRUE)
  stopifnot(attr(result, "status") == 0)
  c(strtrim(result, 32))
}
