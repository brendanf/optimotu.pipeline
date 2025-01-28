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
fastx_gz_index <- function(file) {
  index <- sprintf("%s.fqi", file)
  args <- c(
    "index",
    sprintf("-f=%s", file),
    sprintf("-i=%s", index),
    "-w"
  )
  out = system2("bin/fastqindex_0.9.0b", args)
  stopifnot(out == 0)
  checkmate::assert_file_exists(index, "r")
  index
}

#' @param infile (`character` filename) gzipped fasta or fastq file
#' @param index (`character` filename) index file for `infile`
#' @param i (`integer` vector) indices to extract
#' @param outfile (`character` filename) file to write the extracted sequences to
#' @param renumber (`logical` flag) if `TRUE`, replace the sequence names with
#'   integers, starting at 0.
#' @param append (`logical` flag) if `TRUE`, append to `outfile` if it already
#'   exists, rather than overwriting.
#'
#' @return filename of the output file
#' @describeIn fastx_gz Extract sequences from a gzipped FASTA or FASTQ file
fastx_gz_extract <- function(infile, index, i, outfile, renumber = FALSE, append = FALSE, hash = NULL) {
  checkmate::assert_file_exists(infile, "r")
  checkmate::assert_file_exists(index, "r")
  checkmate::assert_integerish(i, lower = 1)
  checkmate::assert_string(outfile)
  checkmate::assert_flag(renumber)
  checkmate::assert_flag(append)
  if (file.exists(outfile) && !append) unlink(outfile)
  ensure_directory(outfile)
  if (!file.exists(outfile)) file.create(outfile)
  start <- which(i != dplyr::lag(i, 1, -1) + 1L)
  end <- c(start[-1] - 1L, length(i))
  is_fastq <- endsWith(infile, "fastq.gz") || endsWith(infile, "fq.gz")
  command <- sprintf(
    "bin/fastqindex_0.9.0b extract -s=%i -n=%i -e=%i -f=%s -i=%s | tail -n+8",
    i[start] - 1L,
    end - start + 1L,
    if (is_fastq) 4 else 2,
    infile,
    index
  )
  if (renumber) {
    command = paste(
      command,
      sprintf(
        "| awk -v n=%i 'NR%%%i==1{print \"%c\" n; n++; next}; {print}'",
        cumsum(dplyr::lag(end - start + 1L, 1L, 0L)),
        if (is_fastq) 4 else 2,
        if (is_fastq) "@" else ">"
      )
    )
  }
  if (endsWith(outfile, ".gz")) {
    command = paste(command, "| gzip -c -")
  }
  command = paste(command, ">>", outfile)
  result <- vapply(command, system, 0L)
  stopifnot(all(result == 0))
  outfile
}

#' @param infile (`character` filename) gzipped fasta or fastq file
#' @param index (`character` filename) index file for `infile`
#' @param i (`integer` vector) indices to extract
#' @param outfile (`character` filename) file to write the extracted sequences to
#' @param renumber (`logical` flag) if `TRUE`, replace the sequence names with
#' integers, starting at 0.
#' @param append (`logical` flag) if `TRUE`, append to `outfile` if it already
#' exists, rather than overwriting.
#' @param hash (`character` scalar) md5 hash of the infile; ignored but included
#' as a parameter for dependency tracking
#' @param max_gap (`integer` scalar) maximum gap between consecutive indices for
#' a single extraction process
#' @param ncpu (`integer` scalar) number of processes to use
#'
#' @return filename of the output file
#' @describeIn fastx_gz Extract sequences from a gzipped FASTA or FASTQ file
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
  checkmate::assert_file_exists(infile, "r")
  checkmate::assert_file_exists(index, "r")
  checkmate::assert_integerish(i, lower = 1)
  checkmate::assert_string(outfile, null.ok = TRUE)
  checkmate::assert_flag(renumber)
  checkmate::assert_flag(append)
  checkmate::assert_integerish(max_gap, lower = 1)
  isort <- sort(unique(i))
  start <- which(isort > dplyr::lag(isort, 1, -max_gap) + max_gap)
  end <- c(start[-1] - 1L, length(i))
  is_fastq <- endsWith(infile, "fastq.gz") || endsWith(infile, "fq.gz")
  tmpfile <- replicate(length(start), withr::local_tempfile(fileext = ".fasta"))
  processes <- vector("list", length(start))
  for (j in seq_len(min(ncpu, length(start)))) {
    processes[[j]] <- processx::process$new(
      command = "bin/fastqindex_0.9.0b",
      arg = c(
        "extract",
        sprintf("-s=%d", isort[start[j]] - 1L),
        sprintf("-n=%d", isort[end[j]] - isort[start[j]] + 1L),
        sprintf("-e=%d", if (is_fastq) 4 else 2),
        sprintf("-f=%s", infile),
        sprintf("-i=%s", index),
        sprintf("-o=%s", tmpfile[j])
      ),
      stderr = ""
    )
  }
  j <- 1
  fastqindex_return <- 0
  while (j <= length(start) && !is.null(processes[[j]])) {
    fastqindex_return <- fastqindex_return + processes[[j]]$wait()$get_exit_status()
    if (fastqindex_return == 0 && (k <- j + ncpu) <= length(start)) {
      processes[[k]] <- processx::process$new(
        command = "bin/fastqindex_0.9.0b",
        arg = c(
          "extract",
          sprintf("-s=%d", isort[start[k]] - 1L),
          sprintf("-n=%d", isort[end[k]] - isort[start[k]] + 1L),
          sprintf("-e=%d", if (is_fastq) 4 else 2),
          sprintf("-f=%s", infile),
          sprintf("-i=%s", index),
          sprintf("-o=%s", tmpfile[k])
        ),
        stderr = ""
      )
    }
    j <- j + 1
  }
  stopifnot(fastqindex_return == 0)

  included <- unlist(mapply("seq", isort[start], isort[end]))
  selected <- match(i, included)
  seqs <- Biostrings::readBStringSet(tmpfile, seek.first.rec = TRUE)[selected]
  if (renumber) {
    names(seqs) <- as.character(seq_along(seqs))
  }
  if (is.null(outfile)) {
    seqs
  } else {
    write_sequence(seqs, outfile, compress = endsWith(outfile, ".gz"))
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
  is_fastq <- endsWith(infile, "fastq.gz") || endsWith(infile, "fq.gz")
  command <- sprintf(
    "bin/fastqindex_0.9.0b extract -s=%i -n=%i -e=%i -f=%s -i=%s | tail -n+8 | md5sum",
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
