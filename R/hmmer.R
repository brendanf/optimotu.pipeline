#' Try to find the hmmalign executable
#' @return a `character` string giving the path to the hmmalign executable
#' @export
find_hmmalign <- function() {
  find_executable("hmmalign")
}

#' Try to find the hmmsearch executable
#' @return a `character` string giving the path to the hmmsearch executable
#' @export
find_hmmsearch <- function() {
  find_executable("hmmsearch")
}

#' Try to find the nhmmer executable
#' @return a `character` string giving the path to the nhmmer executable
#' @export
find_nhmmer <- function() {
  find_executable("nhmmer")
}

#' Align query sequences to an HMM
#' @param seqs ([`XStringSet`][Biostrings::XStringSet-class], `character`
#'   vector, `data.frame`, or
#'   [`fastqindexr_index`][fastqindexr::create_index()] object) query sequences.
#'   A `character` vector may contain literal sequences (in which case it should
#'   be named) or file paths to one or more FASTA files (possibly gzipped) or
#'   `.fqi` files for [fastqindexr::read_fqi_index()].
#' @param hmm (`character` string giving a file name) HMM for alignment
#' @param outfile (`character` string or vector) output path(s). Normally a
#'   single file; if one path is given, chunk outputs are combined in order.
#'   If more than one path is given, its length must equal `ncpu` (one slot per
#'   potential worker); if fewer chunks are needed (typically because of very
#'   small or empty inputs), only the necessary number of paths are written
#'   When there are no input sequences, an empty alignment is written to the
#'   **first** path in `outfile` (additional paths are ignored).
#' @param outformat (`character` string) output format for alignment; options
#' are `"A2M"` and `"AFA"`, and lower-case versions of these.
#' @param compress (`logical`) if TRUE, compress the output file with gzip
#' @param files (`NULL` or `character`) when `seqs` is a `fastqindexr_index` or
#'   `.fqi` paths, optional per-file paths overriding those stored in the index
#'   (same length as the index file list). Useful after moving inputs or for
#'   `targets` dependency tracking.
#' @param seq_idx (`NULL` or `integer`) optional 1-based indices into the
#'   logical sequence stream (`NULL` means all sequences in order). Applies
#'   after concatenating multiple FASTA inputs, and supports duplicates and
#'   request order. Indexed inputs use
#'   `fastqindexr::extract_sequences_to_file()`
#'   with a live index; plain FASTA paths (no index object) use the same
#'   function with `mode = "sequential"` so records are streamed without
#'   building an index.
#' @param ncpu (`integer`) maximum number of parallel `hmmalign` processes.
#'   The selected sequences are split into up to `ncpu` contiguous chunks of
#'   nearly equal size (fewer when there are fewer sequences than `ncpu`).
#' @param hmmalign (`character`) path to the `hmmalign` executable
#' @param ... ignored; reserved for `targets` dependency tracking (e.g. file
#'   hashes) without changing behavior.
#' @return `character` vector of output file path(s) actually written: length
#'   1 when `outfile` has length 1, otherwise length equal to the number of
#'   materialized input chunks (at least 1, at most `ncpu`).
#' @export
hmmalign <- function(
  seqs,
  hmm,
  outfile,
  outformat = "A2M",
  compress = unique(endsWith(outfile, ".gz")),
  files = NULL,
  seq_idx = NULL,
  ncpu = local_cpus(),
  hmmalign = find_hmmalign(),
  ...
) {
  checkmate::assert_string(hmm)
  checkmate::assert_file_exists(hmm, access = "r")
  checkmate::assert_character(
    outfile,
    min.len = 1L,
    unique = TRUE,
    any.missing = FALSE
  )
  checkmate::assert_count(ncpu, positive = TRUE)
  checkmate::assert(
    checkmate::test_character(outfile, len = 1L),
    checkmate::test_character(outfile, len = ncpu)
  )
  checkmate::assert_choice(outformat, c("A2M", "a2m", "afa", "AFA"))
  checkmate::assert_flag(compress)
  checkmate::assert_count(ncpu)
  checkmate::assert_file_exists(hmmalign, access = "x")
  indexed_like <- inherits(seqs, "fastqindexr_index") ||
    seq_batch_is_fqi_path_set(seqs)
  if (!is.null(files) && !indexed_like) {
    stop(
      "`files` is only valid when `seqs` is a fastqindexr_index or .fqi/.qs2 paths.",
      call. = FALSE
    )
  }
  tmp_parent <- environment()
  tseqs <- seq_batch_make_chunk_files(
    seqs = seqs,
    files = files,
    seq_idx = seq_idx,
    ncpu = ncpu,
    local_envir = tmp_parent
  )
  n <- length(tseqs)
  if (n == 0L) {
    outfile <- outfile[1L]
    ensure_directory(outfile)
    checkmate::assert_path_for_output(outfile, overwrite = TRUE)
    write_sequence(Biostrings::DNAStringSet(), outfile, compress = compress)
    return(outfile)
  }
  if (length(outfile) > n) {
    outfile <- outfile[seq_len(n)]
  }
  ensure_directory(outfile)
  checkmate::assert_path_for_output(outfile, overwrite = TRUE)
  tout <- if (length(outfile) == n && isFALSE(compress)) {
    outfile
  } else {
    replicate(
      n,
      withr::local_tempfile(fileext = ".fasta", .local_envir = tmp_parent)
    )
  }
  mout <- replicate(
    n,
    withr::local_tempfile(fileext = ".fasta", .local_envir = tmp_parent)
  )
  for (i in seq_len(n)) {
    processx::run("mkfifo", mout[i])
  }

  # fmt: skip
  args <- data.frame(
      "--outformat", outformat,
      "--trim",
      "-o", mout,
      hmm,
      tseqs
    )
  args <- as.matrix(args)
  hmmer <- vector("list", n)
  deline <- vector("list", n)
  for (i in seq_len(n)) {
    hmmer[[i]] <- processx::process$new(
      command = hmmalign,
      args = args[i, ],
      supervise = TRUE
    )
    deline[[i]] <- processx::process$new(
      command = "awk",
      args = 'BEGIN{ORS=""};NR>1&&/^>/{print "\\n"};{print};/^>/{print "\\n"};END{print "\\n"}',
      stdin = mout[i],
      stdout = tout[i],
      supervise = TRUE
    )
  }
  hmmer_return <- integer()
  for (i in seq_len(n)) {
    hmmer[[i]]$wait()
    hmmer_return <- union(hmmer[[i]]$get_exit_status(), hmmer_return)
  }
  stopifnot(identical(hmmer_return, 0L))
  for (i in seq_len(n)) {
    deline[[i]]$wait()
  }

  if (isTRUE(compress) && length(outfile) == n) {
    gzip <- vector("list", n)
    for (i in seq_len(n)) {
      gzip[[i]] <- processx::process$new(
        command = "gzip",
        args = c("-c", tout[i]),
        stdout = outfile[i]
      )
    }
    gzip_return <- integer()
    for (i in seq_len(n)) {
      gzip_return <- union(gzip[[i]]$wait()$get_exit_status(), gzip_return)
    }
    stopifnot(identical(gzip_return, 0L))
  } else if (length(outfile) == 1L && n > 1L) {
    fastx_combine(tout, outfile)
  }
  outfile
}

#' Open a HMMER fixed-width "tblout" file
#'
#' This function determines the column widths using the header lines of the
#' @param file (`character` string giving file name, or a connection) file to
#' read
#' @param col_names (`character` vector) names to apply to columns; passed to
#' `readr::read_fwf()`
#' @param col_types (`character` string or object returned by `readr::cols()`)
#' column typed; passed to `readr::read_fwf()`
#' @return a [`tibble`][tibble::tibble()] giving the contents of the file
#' @export
read_hmmer_tblout <- function(file, col_names, col_types) {
  # avoid R CMD check NOTE for undeclared global variables due to NSE
  text <- is_widths <- part <- NULL

  tibble::tibble(
    text = readLines(file),
    is_widths = grepl("^#[- ]+$", text),
    part = cumsum(is_widths)
  ) |>
    dplyr::filter(is_widths | !startsWith(text, "#")) |>
    dplyr::group_split(part, .keep = FALSE) |>
    purrr::discard(\(x) nrow(x) == 1) |>
    purrr::map_dfr(
      \(x) {
        paste(x$text, collapse = "\n") |>
          readr::read_fwf(
            col_positions = stringr::str_locate_all(x$text[1], "#?-+")[[1]] |>
              tibble::as_tibble() |>
              tibble::add_column(col_names = col_names) |>
              do.call(readr::fwf_positions, args = _),
            skip = 1,
            col_types = col_types
          )
      }
    )
}

#' @describeIn read_hmmer_tblout Read a HMMER domain hits file
#' @export
read_domtblout <- function(file) {
  read_hmmer_tblout(
    file,
    col_names = c(
      "seq_name",
      "seq_accno",
      "seq_length",
      "hmm_name",
      "hmm_accno",
      "hmm_length",
      "Evalue",
      "full_score",
      "full_bias",
      "hit_num",
      "total_hits",
      "c_Evalue",
      "i_Evalue",
      "hit_score",
      "hit_bias",
      "hmm_from",
      "hmm_to",
      "seq_from",
      "seq_to",
      "env_from",
      "env_to",
      "acc",
      "description"
    ),
    col_types = "cciccidddiiddddiiiiiidc"
  )
}

#' @describeIn read_hmmer_tblout Read a NHMMER hits file
#' @export
read_dna_tblout <- function(file) {
  read_hmmer_tblout(
    file,
    col_names = c(
      "seq_name",
      "seq_accno",
      "hmm_name",
      "hmm_accno",
      "hmm_from",
      "hmm_to",
      "seq_from",
      "seq_to",
      "env_from",
      "env_to",
      "seq_len",
      "strand",
      "Evalue",
      "bit_score",
      "bias",
      "description"
    ),
    col_types = "cccciiiiiiicdddc"
  )
}

empty_domtblout <- function() {
  tibble::tibble(
    seq_name = character(),
    seq_accno = character(),
    seq_length = integer(),
    hmm_name = character(),
    hmm_accno = character(),
    hmm_length = integer(),
    Evalue = numeric(),
    full_score = numeric(),
    full_bias = numeric(),
    hit_num = integer(),
    total_hits = integer(),
    c_Evalue = numeric(),
    i_Evalue = numeric(),
    hit_score = numeric(),
    hit_bias = numeric(),
    hmm_from = integer(),
    hmm_to = integer(),
    seq_from = integer(),
    seq_to = integer(),
    env_from = integer(),
    env_to = integer(),
    acc = numeric(),
    description = character()
  )
}

empty_dna_tblout <- function() {
  tibble::tibble(
    seq_name = character(),
    seq_accno = character(),
    hmm_name = character(),
    hmm_accno = character(),
    hmm_from = integer(),
    hmm_to = integer(),
    seq_from = integer(),
    seq_to = integer(),
    env_from = integer(),
    env_to = integer(),
    seq_len = integer(),
    strand = character(),
    Evalue = numeric(),
    bit_score = numeric(),
    bias = numeric(),
    description = character()
  )
}

#' Search for subsequences matching one or more HMMs in a set of sequences
#'
#' @param seqs ([`XStringSet`][Biostrings::XStringSet-class], `character`
#'   vector, `data.frame`, or
#'   [`fastqindexr_index`][fastqindexr::create_index()] object) query sequences.
#'   A `character` vector may contain literal sequences (in which case it should
#'   be named) or file paths to one or more FASTA files (possibly gzipped) or
#'   `.fqi` files for [fastqindexr::read_fqi_index()].
#' @param hmm (`character` file name) path to HMM(s) to search for
#' @param files (`NULL` or `character`) when `seqs` is a `fastqindexr_index` or
#'   `.fqi` paths, optional per-file paths overriding those stored in the index
#'   (same length as the index file list). Useful after moving inputs or for
#'   `targets` dependency tracking.
#' @param seq_idx (`NULL` or `integer`) optional 1-based indices into the
#'   logical sequence stream (`NULL` means all sequences in order). Applies
#'   after concatenating multiple FASTA inputs, and supports duplicates and
#'   request order.
#' @param ncpu (`integer`) maximum number of parallel `hmmsearch` processes.
#'   The selected sequences are split into up to `ncpu` contiguous chunks of
#'   nearly equal size (fewer when there are fewer sequences than `ncpu`).
#' @param hmmsearch (`character`) path to the `hmmsearch` executable
#' @param ... ignored; reserved for `targets` dependency tracking (e.g. file
#'   hashes) without changing behavior.
#' @return a [`tibble`][tibble::tibble()] listing the HMM hits. Empty sequence
#' inputs return an empty tibble with the standard HMMER columns.
#' @export
hmmsearch <- function(
  seqs,
  hmm,
  files = NULL,
  seq_idx = NULL,
  ncpu = local_cpus(),
  hmmsearch = find_hmmsearch(),
  ...
) {
  checkmate::assert_string(hmm)
  checkmate::assert_file_exists(hmm, access = "r")
  checkmate::assert_count(ncpu)
  checkmate::assert_file_exists(hmmsearch, access = "x")
  indexed_like <- inherits(seqs, "fastqindexr_index") ||
    seq_batch_is_fqi_path_set(seqs)
  if (!is.null(files) && !indexed_like) {
    stop(
      "`files` is only valid when `seqs` is a fastqindexr_index or .fqi/.qs2 paths.",
      call. = FALSE
    )
  }
  tmp_parent <- environment()
  tseqs <- seq_batch_make_chunk_files(
    seqs = seqs,
    files = files,
    seq_idx = seq_idx,
    ncpu = ncpu,
    local_envir = tmp_parent
  )
  n <- length(tseqs)
  if (n == 0L) {
    return(empty_domtblout())
  }
  outfile <- replicate(
    n,
    withr::local_tempfile(fileext = ".hmmout", .local_envir = tmp_parent)
  )
  # fmt: skip
  args <- data.frame(
    "--noali",
    "--notextw",
    "--domtblout", outfile,
    hmm,
    tseqs
  )
  args <- as.matrix(args)

  hmmer <- vector("list", n)
  for (i in seq_len(n)) {
    hmmer[[i]] <- processx::process$new(
      command = hmmsearch,
      args = args[i, ],
      supervise = TRUE
    )
  }
  hmmer_return <- integer()
  for (i in seq_len(n)) {
    hmmer[[i]]$wait()
    hmmer_return <- union(hmmer[[i]]$get_exit_status(), hmmer_return)
    stopifnot(identical(hmmer_return, 0L))
  }
  purrr::map_dfr(outfile, read_domtblout)
}

#' Search for subsequences matching one or more nucleotide HMMs in a set of
#' sequences
#'
#' @param seqs ([`XStringSet`][Biostrings::XStringSet-class], `character`
#'   vector, `data.frame`, or
#'   [`fastqindexr_index`][fastqindexr::create_index()] object) query sequences.
#'   A `character` vector may contain literal sequences (in which case it should
#'   be named) or file paths to one or more FASTA files (possibly gzipped) or
#'   `.fqi` files for [fastqindexr::read_fqi_index()].
#' @param hmm (`character` file name) path to HMM(s) to search for
#' @param files (`NULL` or `character`) when `seqs` is a `fastqindexr_index` or
#'   `.fqi` paths, optional per-file paths overriding those stored in the index
#'   (same length as the index file list). Useful after moving inputs or for
#'   `targets` dependency tracking.
#' @param seq_idx (`NULL` or `integer`) optional 1-based indices into the
#'   logical sequence stream (`NULL` means all sequences in order). Applies
#'   after concatenating multiple FASTA inputs, and supports duplicates and
#'   request order.
#' @param ncpu (`integer`) number of threads passed to nhmmer (`--cpu`).
#'   Parallelism is only inside nhmmer; the query is not split across
#'   processes.
#' @param nhmmer (`character`) path to the `nhmmer` executable
#' @param ... ignored; reserved for `targets` dependency tracking (e.g. file
#'   hashes) without changing behavior.
#' @return a [`tibble`][tibble::tibble()] like [hmmsearch()] tblout parsing.
#'   Empty sequence inputs return an empty tibble with the DNA tblout columns.
#' @export
nhmmer <- function(
  seqs,
  hmm,
  files = NULL,
  seq_idx = NULL,
  ncpu = local_cpus(),
  nhmmer = find_nhmmer(),
  ...
) {
  checkmate::assert_string(hmm)
  checkmate::assert_file_exists(hmm, access = "r")
  checkmate::assert_count(ncpu)
  checkmate::assert_file_exists(nhmmer, access = "x")
  if (is.list(seq_idx) && length(seq_idx) > 1L) {
    stop(
      "`seq_idx` must not be a list with more than one partition.",
      call. = FALSE
    )
  }
  indexed_like <- inherits(seqs, "fastqindexr_index") ||
    seq_batch_is_fqi_path_set(seqs)
  if (!is.null(files) && !indexed_like) {
    stop(
      "`files` is only valid when `seqs` is a fastqindexr_index or .fqi/.qs2 paths.",
      call. = FALSE
    )
  }
  tmp_parent <- environment()
  qfiles <- seq_batch_make_chunk_files(
    seqs = seqs,
    files = files,
    seq_idx = seq_idx,
    ncpu = 1L,
    local_envir = tmp_parent
  )
  if (length(qfiles) == 0L) {
    return(empty_dna_tblout())
  }
  tseqs <- qfiles[[1L]]
  if (length(Biostrings::fasta.seqlengths(tseqs)) == 0L) {
    return(empty_dna_tblout())
  }
  outfile <- withr::local_tempfile(
    fileext = ".hmmout",
    .local_envir = tmp_parent
  )
  # fmt: skip
  args <- c(
    "--noali",
    "--notextw",
    "--tblout", outfile,
    "--watson",
    "--cpu", ncpu,
    hmm,
    tseqs
  )
  processx::run(
    command = nhmmer,
    args = args,
    error_on_status = TRUE
  )
  read_dna_tblout(outfile)
}
