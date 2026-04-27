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
#' string giving a FASTA file name (possibly gzipped), `character` vector, or
#' `data.frame`) sequences to align
#' @param hmm (`character` string giving a file name) HMM for alignment
#' @param outfile (`character` string) file name for output
#' @param outformat (`character` string) output format for alignment; options
#' are `"A2M"` and `"AFA"`, and lower-case versions of these.
#' @param compress (`logical`) if TRUE, compress the output file with gzip
#' @return a `character` string giving the output file. Empty sequence inputs
#' produce empty alignment output files with matching format/compression.
#' @export
hmmalign <- function(
  seqs,
  hmm,
  outfile,
  outformat = "A2M",
  compress = endsWith(outfile, ".gz")
) {
  checkmate::assert_string(hmm)
  checkmate::assert_file_exists(hmm, access = "r")
  ensure_directory(outfile)
  checkmate::assert_path_for_output(outfile, overwrite = TRUE)
  checkmate::assert_choice(outformat, c("A2M", "a2m", "afa", "AFA"))
  checkmate::assert_flag(compress)
  exec <- find_hmmalign()
  checkmate::assert_file_exists(exec, access = "x")
  if (checkmate::test_file_exists(seqs, "r")) {
    tseqs <- seqs
    n <- length(seqs)
    is_gz <- endsWith(tseqs, ".gz")
    if (any(is_gz)) {
      tseqs_gz <- tseqs[is_gz]
      tseqs_ungz <- replicate(
        length(tseqs_gz),
        withr::local_tempfile(
          fileext = ".fasta"
        )
      )
      for (i in seq_along(tseqs_gz)) {
        write_sequence(
          Biostrings::readDNAStringSet(tseqs_gz[[i]]),
          tseqs_ungz[[i]]
        )
      }
      tseqs[is_gz] <- tseqs_ungz
    }
  } else if (
    checkmate::test_list(
      seqs,
      types = c("character", "XStringSet", "data.frame")
    )
  ) {
    n <- length(seqs)
    tseqs <- replicate(n, withr::local_tempfile(fileext = ".fasta"))
    purrr::pwalk(list(seq = seqs, fname = tseqs), write_sequence)
  } else {
    checkmate::assert_multi_class(
      seqs,
      c("data.frame", "character", "XStringSet")
    )
    n <- 1
    tseqs <- withr::local_tempfile(fileext = ".fasta")
    write_sequence(seqs, tseqs)
  }
  checkmate::assert(
    length(outfile) == 1,
    length(outfile) == n
  )
  if (length(outfile) == 1 && n > 1) {
    tout <- replicate(n, withr::local_tempfile(fileext = ".fasta"))
  } else {
    tout <- outfile
  }
  if (compress) {
    if (identical(tout, outfile)) {
      tout <- replicate(n, withr::local_tempfile(fileext = ".fasta"))
    }
  }
  is_empty <- vapply(
    tseqs,
    function(x) length(Biostrings::fasta.seqlengths(x)) == 0L,
    logical(1)
  )
  if (any(is_empty)) {
    purrr::walk(
      tout[is_empty],
      function(x) write_sequence(Biostrings::DNAStringSet(), x)
    )
  }
  run_idx <- which(!is_empty)
  if (length(run_idx) > 0L) {
    mout <- replicate(
      length(run_idx),
      withr::local_tempfile(fileext = ".fasta")
    )
    for (i in seq_along(run_idx)) {
      processx::run("mkfifo", mout[i])
    }

    # fmt: skip
    args <- data.frame(
      "--outformat", outformat,
      "--trim",
      "-o", mout,
      hmm,
      tseqs[run_idx]
    )
    args <- as.matrix(args)
    hmmer <- vector("list", length(run_idx))
    deline <- vector("list", length(run_idx))
    for (i in seq_along(run_idx)) {
      hmmer[[i]] <- processx::process$new(
        command = exec,
        args = args[i, ],
        supervise = TRUE
      )
      deline[[i]] <- processx::process$new(
        command = "awk",
        args = 'BEGIN{ORS=""};NR>1&&/^>/{print "\\n"};{print};/^>/{print "\\n"};END{print "\\n"}',
        stdin = mout[i],
        stdout = tout[run_idx[i]],
        supervise = TRUE
      )
    }
    hmmer_return <- integer()
    for (i in seq_along(run_idx)) {
      hmmer[[i]]$wait()
      hmmer_return <- union(hmmer[[i]]$get_exit_status(), hmmer_return)
    }
    stopifnot(identical(hmmer_return, 0L))
    for (i in seq_along(run_idx)) {
      deline[[i]]$wait()
    }
  }

  if (compress && length(outfile) == n) {
    gzip <- vector("list", n)
    for (i in seq_len(n)) {
      gzip[[i]] <- processx::process$new(
        command = "gzip",
        error_on_status = TRUE,
        args = c("-c", tout[i]),
        stdout = outfile[i]
      )
    }
    gzip_return <- integer()
    for (i in seq_len(n)) {
      gzip_return <- union(gzip[[i]]$wait()$status, gzip_return)
    }
    stopifnot(identical(gzip_return, 0L))
  } else if (length(outfile) < n) {
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
#' file name of a FASTA file (possibly gzipped), `character` vector, or
#' `data.frame`) sequences to search
#' @param hmm (`character` file name) path to HMM(s) to search for
#' @return a [`tibble`][tibble::tibble()] listing the HMM hits. Empty sequence
#' inputs return an empty tibble with the standard HMMER columns.
#' @export
hmmsearch <- function(seqs, hmm) {
  checkmate::assert_string(hmm)
  checkmate::assert_file_exists(hmm, access = "r")
  exec <- find_hmmsearch()
  checkmate::assert_file_exists(exec, access = "x")
  if (checkmate::test_file_exists(seqs, "r")) {
    tseqs <- seqs
    n <- length(seqs)
  } else if (
    checkmate::test_list(
      seqs,
      types = c("character", "XStringSet", "data.frame")
    )
  ) {
    n <- length(seqs)
    tseqs <- replicate(n, withr::local_tempfile(fileext = ".fasta"))
    purrr::pwalk(list(seq = seqs, fname = tseqs), write_sequence)
  } else {
    checkmate::assert_multi_class(
      seqs,
      c("data.frame", "character", "XStringSet")
    )
    n <- 1
    tseqs <- withr::local_tempfile(fileext = ".fasta")
    write_sequence(seqs, tseqs)
  }
  is_empty <- vapply(
    tseqs,
    function(x) length(Biostrings::fasta.seqlengths(x)) == 0L,
    logical(1)
  )
  run_idx <- which(!is_empty)
  if (length(run_idx) == 0L) {
    return(empty_domtblout())
  }
  outfile <- replicate(n, withr::local_tempfile(fileext = ".hmmout"))
  # fmt: skip
  args <- data.frame(
    "--noali",
    "--notextw",
      "--domtblout", outfile[run_idx],
    hmm,
      tseqs[run_idx]
  )
  args <- as.matrix(args)

  hmmer <- vector("list", length(run_idx))
  for (i in seq_along(run_idx)) {
    hmmer[[i]] <- processx::process$new(
      command = exec,
      args = args[i, ],
      supervise = TRUE
    )
  }
  hmmer_return <- integer()
  for (i in seq_along(run_idx)) {
    hmmer[[i]]$wait()
    hmmer_return <- union(hmmer[[i]]$get_exit_status(), hmmer_return)
    stopifnot(identical(hmmer_return, 0L))
  }
  purrr::map_dfr(
    outfile[run_idx],
    read_domtblout
  )
}

#' Search for subsequences matching one or more nucleotide HMMs in a set of
#' sequences
#'
#' @inherit hmmsearch params return
#' @param ncpu (`integer`) number of threads to use for searching
#' @export
nhmmer <- function(seqs, hmm, ncpu = local_cpus()) {
  checkmate::assert_string(hmm)
  checkmate::assert_file_exists(hmm, access = "r")
  checkmate::assert_count(ncpu)
  exec <- find_nhmmer()
  checkmate::assert_file_exists(exec, access = "x")
  if (length(seqs) == 1 && checkmate::test_file_exists(seqs, "r")) {
    tseqs <- seqs
  } else {
    checkmate::assert_multi_class(
      seqs,
      c("data.frame", "character", "XStringSet")
    )
    tseqs <- withr::local_tempfile(fileext = ".fasta")
    write_sequence(seqs, tseqs)
  }
  if (length(Biostrings::fasta.seqlengths(tseqs)) == 0L) {
    return(empty_dna_tblout())
  }
  outfile <- withr::local_tempfile(fileext = ".hmmout")
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
    command = exec,
    args = args,
    error_on_status = TRUE
  )
  read_dna_tblout(outfile)
}
