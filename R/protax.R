#' Taxonomically identify sequences using the PERL implementation of ProtaxFungi
#'
#' @param seqs ([`DNAStringSet`][Biostrings::XStringSet-class], `character`
#' path to a FASTA file, `character` vector, or `data.frame`) sequences to
#' taxonomically identify
#' @param outdir (`character` string) directory to write output to. If the'
#' directory exists, all existing files in it will be deleted.
#' @param modeldir (`character` string) directory containing trained model files
#' for ProtaxFungi
#' @param ncpu (`integer`) number of threads to use for calls to USEARCH
#' @param script (`character` string) path to the ProtaxFungi script
run_protax <- function(seqs, outdir, modeldir, ncpu = local_cpus(),
                       script = "scripts/runprotax") {
  checkmate::assert_directory_exists(modeldir)
  checkmate::assert_integer(ncpu, lower = 1)
  checkmate::assert_file(script, access = "x")
  checkmate::assert_character(outdir)
  if (dir.exists(outdir)) unlink(outdir, recursive = TRUE)
  dir.create(outdir)
  if (length(seqs) == 1 && file.exists(seqs)) {
    if (seqs != file.path(outdir, "all.fa"))
      file.copy(seqs, file.path(outdir, "all.fa"))
  } else {
    write_sequence(seqs, file.path(outdir, "all.fa"))
  }
  status <- system2(
    script,
    c(outdir, modeldir, ncpu)
  )
  stopifnot(status == 0)
  list.files(outdir, full.names = TRUE)
}

#' Parse the output of ProtaxFungi
#' @param nameprob (`character` string) path to the "nameprob*" files output by
#' ProtaxFungi
#' @param id_is_int (`logical`) if `TRUE`, sequence labels should be interpreted
#' as integers
#' @return a `data.frame` with columns `seq_id` (or `seq_idx` if `id_is_int` is
#' `TRUE`), `rank`, `taxonomy`, `prob`
#' @export
parse_protax_nameprob <- function(nameprob, id_is_int = FALSE) {
  # avoid R CMD check NOTE about global variables due to NSE
  name <- value <- taxon <- prob <- NULL

  checkmate::assert_flag(id_is_int)
  id_col <- if (isTRUE(id_is_int)) "seq_idx" else "seq_id"
  id_col_name <- as.symbol(id_col)
  `names<-`(nameprob, basename(nameprob)) |>
    lapply(readLines) |>
    tibble::enframe() |>
    tidyr::extract(
      name,
      into = "rank",
      regex = "query(\\d+)\\.nameprob",
      convert = TRUE
    ) |>
    tidyr::unchop(value) |>
    dplyr::mutate(
      rank = rank2factor(tax_ranks()[rank]),
      value = gsub("([^\t]+)\t([0-9.]+)", "\\1:\\2", value) |>
        gsub("(:[0-9.]+)\t", "\\1;", x = _)
    ) |>
    tidyr::separate(value, into = c(id_col, "nameprob"), sep = "\t", fill = "right") |>
    tidyr::separate_rows(nameprob, sep = ";") |>
    tidyr::separate(nameprob, into = c("name", "prob"), sep = ":", convert = TRUE) |>
    tidyr::extract(name, into = c("parent_taxonomy", "taxon"), regex = "(.+),([^,]+)$") |>
    dplyr::mutate(
      !!id_col_name := if (id_is_int) as.integer(!!id_col_name) else !!id_col_name,
      prob = ifelse(is.na(taxon), 0, prob)
    ) |>
    dplyr::arrange(!!id_col_name, dplyr::desc(rank), dplyr::desc(prob))
}

#' Taxonomically identify sequences using ProtaxAnimal
#'
#' vectorized on aln_seqs
#' attempts to run them ALL in parallel, be careful!
#'
#' @param aln_seqs (`character` vector) file name(s) of FASTA files containing
#' sequences to taxonomically identify. Sequences should all be aligned to the
#' reference alignment used to train Protax, for instance using `hmmalign()`.
#' If they end in `".gz"` then they are assumed to be compressed with gzip.
#' @param modeldir (`character` string) directory containing trained ProtaxA
#' model files.  These should include `"taxonomy.priors"`, `"refs.aln"`,
#' `"model.rseqs.numeric"`, `"model.pars"`, `"model.scs"`.
#' @param min_p (`numeric`) minimum probability to recurse to the next rank
#' @param rep_p (`numeric`) minimum probability to include an identification in
#' the output
#' @param strip_inserts (`logical`) if `TRUE`, `aln_seqs` may include insertions
#' in lower-case, so strip these.
#' @param id_is_int (`logical`) if `TRUE`, sequence labels should be interpreted
#' as integers
#' @param info (`logical`) if `TRUE`, return information about which reference
#' sequences were the top two matches responsible for each identification
#' @param options (`character` vector) additional command line options to
#' Protax
#' @return a `data.frame` with columns `seq_id`, `rank`, `taxonomy`, `prob`,
#' and if `info` is `TRUE`, `best_id`, `best_dist`, `second_id`, `second_dist`
#' @export
run_protax_animal <- function(aln_seqs, modeldir, min_p = 0.1, rep_p = 0.01,
                              strip_inserts = FALSE, id_is_int = FALSE,
                              info = FALSE, options = character()) {
  checkmate::assert_file_exists(aln_seqs, access = "r")
  checkmate::assert_directory_exists(modeldir)
  priors <- file.path(modeldir, "taxonomy.priors")
  checkmate::assert_file_exists(priors, access = "r")
  refs <- file.path(modeldir, "refs.aln")
  checkmate::assert_file_exists(refs, access = "r")
  rseqs <- file.path(modeldir, "model.rseqs.numeric")
  checkmate::assert_file_exists(rseqs, access = "r")
  pars <- file.path(modeldir, "model.pars")
  checkmate::assert_file_exists(pars, access = "r")
  scs <- file.path(modeldir, "model.scs")
  checkmate::assert_file_exists(scs, access = "r")
  checkmate::assert_flag(info)
  executable <- find_executable(if (info) "classify_info" else "classify_v2")
  checkmate::check_number(min_p, lower = 0, upper = 1, finite = TRUE)
  checkmate::check_number(rep_p, lower = 0, upper = min_p, finite = TRUE)
  checkmate::assert_flag(strip_inserts)
  checkmate::assert_flag(id_is_int)
  checkmate::assert_character(options)
  args <- c("-t", rep_p, options, priors, refs, rseqs, pars, scs, as.character(min_p))

  is_gz <-endsWith(aln_seqs, ".gz")
  stopifnot(all(is_gz) | all(!is_gz))
  is_gz <- all(is_gz)

  n <- length(aln_seqs)
  if (strip_inserts | is_gz) {
    prepipe <- vector("list", n)
    protax_in <- replicate(n, withr::local_tempfile(fileext = ".fasta"))
    pipecommand <- if (is_gz) "zcat" else "cat"
    pipecommand <- paste(pipecommand, aln_seqs)
    if (strip_inserts) pipecommand <- paste(pipecommand, "| tr -d 'acgt'")
    pipecommand <- paste(pipecommand, ">", protax_in)
  } else {
    protax_in <- aln_seqs
  }

  protax <- vector("list", n)
  outfiles <- replicate(n, withr::local_tempfile())
  for (i in seq_len(n)) {
    if (strip_inserts | is_gz) {
      system(pipecommand[i])
    }
    protax[[i]] <- processx::process$new(
      command = executable,
      args = c(args, protax_in[i]),
      stdout = outfiles[i]
    )
  }
  protax_exit_status = 0L
  output <- vector("list", n)
  for (i in seq_len(n)) {
    protax[[i]]$wait()
    protax_exit_status <- max(protax_exit_status, protax[[i]]$get_exit_status())
    stopifnot(protax_exit_status == 0L)
    output[[i]] <- readr::read_delim(
      outfiles[i],
      col_names = c(
        if (id_is_int) "seq_idx" else "seq_id",
        "rank",
        "taxonomy",
        "prob",
        if (info) c("best_id", "best_dist", "second_id", "second_dist") else NULL
      ),
      col_types = paste0(
        if (id_is_int) "i" else "c",
        "icd",
        if (info) "-id-id-" else ""
      )
    )
  }
  dplyr::bind_rows(output)
}

#' Find the best hit for each query sequence in a reference alignment
#'
#' @param aln_query (`character` vector) file name(s) of FASTA files containing
#' query sequences. If they end in `".gz"` then they are assumed to be
#' compressed with gzip.
#' @param aln_ref (`character` vector) file name(s) of FASTA files containing
#' reference sequences. If they end in `".gz"` then they are assumed to be
#' compressed with gzip.
#' @param options (`character` vector) additional command line options to
#' Protax
#' @param query_id_is_int (`logical`) if `TRUE`, query sequence labels should be
#' interpreted as integers
#' @param ref_id_is_int (`logical`) if `TRUE`, reference sequence labels should
#' be interpreted as integers
#' @param command (`character`) the Protax command to run, either `"dist_best"`
#' or `"dist_bipart"`
#' @return a `data.frame` with columns `seq_id` (or `seq_idx` if
#' `query_id_is_int` is `TRUE`), `ref_id` (or `ref_idx` if `ref_id_is_int` is
#' `TRUE`), and `dist`.  When command in `"dist_best"`, each value in
#' `seq_id`/`seq_idx` occurs only once, while for `"dist_bipart"` each
#' `seq_id`/`seq_idx` may occur multiple times.
#' @export
run_protax_besthit <- function(aln_query, aln_ref, options = character(),
                               query_id_is_int = TRUE, ref_id_is_int = TRUE,
                               command = "dist_best") {
  checkmate::assert_file_exists(aln_query, access = "r")
  n_query <- length(aln_query)
  checkmate::assert_file_exists(aln_ref, access = "r")
  n_ref <- length(aln_ref)
  executable <- find_executable(command)
  checkmate::assert_character(options)
  checkmate::assert_flag(query_id_is_int)
  checkmate::assert_flag(ref_id_is_int)

  is_gz_ref <- endsWith(aln_ref, ".gz")
  stopifnot(all(is_gz_ref) | all(!is_gz_ref))
  is_gz_ref <- all(is_gz_ref)

  if (n_ref > 1) {
    ref <- replicate(n_query, withr::local_tempfile(fileext = ".fasta"))
    refcommand <- if (is_gz_ref) "zcat" else "cat"
    for (i in seq_len(n_query)) {
      system2("mkfifo", args = ref[i])
      system2(
        command = refcommand,
        args = aln_ref,
        stdout = ref[i],
        wait = FALSE
      )
    }
  } else {
    ref <- rep_len(aln_ref, n_query)
  }

  besthit <- vector("list", n_query)
  outfiles <- replicate(n_query, withr::local_tempfile())
  for (i in seq_len(n_query)) {
    besthit[[i]] <- processx::process$new(
      command = executable,
      args = c(options, ref[i], aln_query[i]),
      stdout = outfiles[i]
    )
  }
  besthit_exit_status = 0L
  output <- vector("list", n_query)
  for (i in seq_len(n_query)) {
    besthit[[i]]$wait()
    besthit_exit_status <- max(besthit_exit_status, besthit[[i]]$get_exit_status())
    stopifnot(besthit_exit_status == 0L)
    besthit[[i]] <- readr::read_delim(
      outfiles[i],
      col_names = c(
        if (query_id_is_int) "seq_idx" else "seq_id",
        if (ref_id_is_int) "ref_idx" else "ref_id",
        "dist"
      ),
      col_types = paste0(
        if (query_id_is_int) "i" else "c",
        if (ref_id_is_int) "i" else "c",
        "d"
      ),
      delim = " "
    )
  }
  dplyr::bind_rows(besthit)
}

#' @describeIn run_protax_besthit Generate a sparse bipartite distance matrix
#' @param max_d (`numeric`) maximum distance to include in the distance matrix
run_protax_bipart <- function(aln_query, aln_ref, max_d = 0.2,
                              query_id_is_int = TRUE, ref_id_is_int = TRUE) {
  checkmate::assert_number(max_d, lower = 0, upper = 1)
  run_protax_besthit(
    aln_query = aln_ref,
    aln_ref = aln_query,
    options = as.character(max_d),
    query_id_is_int = query_id_is_int,
    ref_id_is_int = ref_id_is_int,
    command = "dist_bipart"
  )
}

#' Closed-reference single-linkage clustering of aligned sequences
#' @param infile (`character`) file name of gzipped FASTA file containing
#' sequences to cluster
#' @param index (`character`) file name of index for `infile`, as created by
#' `fastx_gz_index()`
#' @param i (`integer`) integer indices of sequences to cluster
#' @param unknowns (`logical`) logical vector indicating which sequences in `i`
#' are considered unknowns (the rest are references)
#' @param thresh (`numeric`) clustering similarity threshold as a percentage,
#' i.e., theoretically 0.0 to 100.0, but only values between 50.0 and 100.0 are
#' accepted, to avoid accidentally supplying distances (likely to be less than
#' 50) or fractions (always between 0.0 and 1.0).
#' @param seq_width (`integer`) width of alignment
#' @param ncpu (`integer`) number of threads to use
#' @param max_gap (`integer`) passed to `fastx_gz_random_access_extract()`
#' @return a `data.frame` with columns `seq_id`, `cluster`, and `dist`
protax_besthit_closedref <- function(infile, index, i, unknowns, thresh,
                                     seq_width, ncpu = local_cpus(),
                                     max_gap = 100L) {
  # avoid R CMD check NOTE about global variables due to NSE
  dist <- ref_id <- NULL

  # check arguments
  checkmate::assert_file(infile, "r")
  checkmate::assert_file(index, "r")
  checkmate::assert_integerish(i)
  checkmate::assert_logical(unknowns, len = length(i))
  checkmate::assert_number(thresh, lower = 50, upper = 100)
  checkmate::assert_integer(seq_width, lower = 1)
  checkmate::assert_integer(ncpu, lower = 1)
  checkmate::assert_integer(max_gap, lower = 1)

  seqs <- fastx_gz_random_access_extract(
    infile = infile,
    index = index,
    i = i,
    ncpu = ncpu,
    max_gap = max_gap
  )

  queries <- names(seqs)[unknowns]
  refs <- names(seqs)[!unknowns]

  out <- NULL
  while (length(queries) > 0 && length(refs) > 0) {
    query_file <- write_sequence(
      seqs[queries],
      withr::local_tempfile(fileext=".fasta")
    ) |>
      fastx_split(
        # only parallelize if it is likely to be worth it.
        n = if (length(queries) > 1000) local_cpus() else 1L,
        outroot = tempfile(tmpdir = withr::local_tempdir()),
        compress = FALSE
      )

    ref_file <- write_sequence(
      seqs[refs],
      withr::local_tempfile(fileext=".fasta")
    )

    newout <-
      run_protax_besthit(
        aln_query = query_file,
        aln_ref = ref_file,
        options = c(
          "-m", "100",
          "-l", seq_width,
          "-r", as.character(length(refs))
        ),
        query_id_is_int = FALSE,
        ref_id_is_int = FALSE
      ) |>
      dplyr::filter(dist <= 1 - thresh/100)

    refs <- unique(newout$seq_id)
    queries <- setdiff(queries, refs)
    if (is.null(out)) {
      out <- dplyr::rename(newout, cluster = ref_id)
    } else {
      newout <- dplyr::left_join(newout, out, by = c("ref_id" = "seq_id")) |>
        dplyr::select(-ref_id)
      out <- dplyr::bind_rows(out, newout)
    }
  }
  out
}

#' Single-linkage clustering of aligned sequences
#'
#' @param aln_seq (`character`) file name of gzipped FASTA file containing
#' sequences to cluster
#' @param aln_index (`character`) file name of index for `aln_seq`, as created
#' by `fastx_gz_index()`
#' @param which (`integer`) integer indices of sequences to cluster
#' @param thresh (`numeric`) clustering similarity threshold as a percentage,
#' i.e., 0.0 to 100.0
#' @param aln_len (`integer`) width of alignment
#' @return a `data.frame` with columns `seq_id`, `cluster`, and `dist`
#' @export
seq_cluster_protax <- function(aln_seq, aln_index, which, thresh, aln_len) {
  nslice <- floor(sqrt(local_cpus()-1))
  allseq <- fastx_gz_random_access_extract(
    infile = aln_seq,
    index = aln_index,
    i = which
  )
  names(allseq) <- as.character(seq_along(allseq) - 1L)
  if (length(which) > 1000) {
    seq <- character(nslice)
    spl <- sort(rep_len(seq_len(nslice), length(which)))
    i <- split(seq_len(length(which)), spl)
    for (j in seq_len(nslice)) {
      seq[j] <- write_sequence(
        allseq[i[[j]]],
        fname = withr::local_tempfile(fileext = ".fasta")
      )
    }
    i_half <- lapply(i, \(x) split(x, rep(c(1, 2), length.out = length(x))))
    i_half <- do.call(c, args = i_half)
    seq_half <- character(nslice*2)
    for (j in seq(3, nslice*2)) {
      seq_half[j] <- write_sequence(
        allseq[i_half[[j]]],
        fname = tempfile(tmpdir = withr::local_tempdir())
      )
    }
  } else {
    nslice <- 1
    i <- list(seq_along(allseq))
    seq <- write_sequence(
      allseq,
      fname = withr::local_tempfile(fileext = ".fasta")
    )
  }
  distmx <- withr::local_tempfile()
  system2("mkfifo", distmx)
  for (j in seq_len(nslice)) {
    system2(
      find_executable("dist_matrix"),
      args = c(
        "-l", aln_len,
        "-i", length(i[[j]]),
        "-m", 100,
        max(thresh),
        seq[j]
      ),
      stdout = distmx,
      wait = FALSE
    )
    if (j < nslice) {
      for (k in (2*j+1):(2*nslice)) {
        system2(
          find_executable("dist_bipart"),
          args = c(
            "-l", aln_len,
            "-i", length(i[[j]]),
            "-m", 100,
            max(thresh),
            seq_half[k],
            seq[j]
          ),
          stdout = distmx,
          wait = FALSE
        )
      }
    }
  }
  optimotu::distmx_cluster(
    distmx = distmx,
    names = as.character(which),
    threshold_config = optimotu::threshold_set(thresh),
    clust_config = optimotu::clust_tree(),
    parallel_config = optimotu::parallel_concurrent(max(1, local_cpus() - nslice))
  )
}

#' Parse the output of ProtaxAnimal
#' @param x (`character`) raw output of ProtaxAnimal, one line per vector element
#' @return a `data.frame` with columns `seq_id`, `rank`, `parent_taxonomy`,
#' `taxon`, and `prob`
#' @export
parse_protaxAnimal_output <- function(x) {
  # avoid R CMD check NOTE about global variables due to NSE
  output <- assignment <- taxonomy <- parent_taxonomy <- prob <- taxon <-
    seq_id <- NULL
  out <-
    tibble::tibble(output = x) |>
    tidyr::separate(
      output,
      into = c("seq_id", "assignment"),
      extra = "merge",
      fill = "right",
      convert = TRUE
    ) |>
    dplyr::mutate(
      assignment = gsub("([^ ]+) ([0-9.]+)", "\\1\x1f\\2", assignment)
    ) |>
    tidyr::separate_longer_delim(assignment, " ") |>
    tidyr::separate(
      assignment,
      into = c("taxonomy", "prob"),
      sep = "\x1f",
      convert = TRUE
    ) |>
    dplyr::mutate(
      rank = int2rankfactor(
        stringr::str_count(taxonomy, ",") + 1L + length(known_ranks())
      )
    ) |>
    tidyr::extract(taxonomy, into = c("parent_taxonomy", "taxon"), regex = "(?:(.+),)?([^,]+)$") |>
    dplyr::mutate(
      parent_taxonomy = paste(
        paste(known_taxa(), collapse = ","),
        parent_taxonomy,
        sep = ","
      ) |>
        trimws(whitespace = ","),
      taxon = dplyr::na_if(taxon, "unk")
    ) |>
    dplyr::select(seq_id, rank, parent_taxonomy, taxon, prob) |>
    dplyr::arrange(seq_id, rank)

  if (is.integer(out$seq_id)) out <- dplyr::rename(out, seq_idx = seq_id)
  out
}
