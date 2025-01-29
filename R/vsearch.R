#' Try to find the vsearch executable
#' @return (`character` string) the full path to the vsearch executable.
#' @export
find_vsearch <- function() {
  find_executable("vsearch")
}

#' "usearch_global" function of vsearch
#'
#' @param query (`data.frame`, [`DNAStringSet`][Biostrings::XStringSet-class], `character` vector, or
#' file name) query sequences
#' @param ref (`data.frame`, [`DNAStringSet`][Biostrings::XStringSet-class], `character` vector, or
#' file name) reference sequences
#' @param threshold (`numeric` scalar) identity threshold, in range 0.0-1.0
#' @param global (`logical` flag) if `TRUE`, end gaps and internal gaps are
#' penalized equally.  Otherwise end gaps are not penalized.
#' @param ncpu (`integer` count) number of threads to use
#' @param id_is_int (`logical` flag) if `TRUE`, return the sequence IDs as
#' integers
#'
#' @return `tibble::tibble` with columns `seq_id`, `clust`, and `dist`, where `seq_id` is
#' the name of a sequence from `query`, `clust` is the closest match to that
#' sequence in `ref`, and `dist` is the distance between them
#' @export
vsearch_usearch_global <- function(query, ref, threshold, global = TRUE,
                                   ncpu = local_cpus(), id_is_int = FALSE) {
  checkmate::check_flag(id_is_int)
  if (is.character(query) && length(query) == 1 && file.exists(query)) {
    tquery <- query
  } else {
    tquery <- withr::local_tempfile(pattern = "query", fileext = ".fasta")
    write_sequence(query, tquery)
  }
  if (is.character(ref) && length(ref) == 1 && file.exists(ref)) {
    tref <- ref
  } else {
    tref <- withr::local_tempfile(pattern = "ref", fileext = ".fasta")
    write_sequence(ref, tref)
  }
  checkmate::assert_flag(global)
  gap <- if (global) "1" else "1I/0E"
  uc = system(
    paste(
      find_vsearch(),
      "--usearch_global", tquery,
      "--db", tref,
      "--id", threshold,
      "--uc", "-",
      "--maxaccepts", "100",
      "--top_hits_only",
      "--threads", ncpu,
      "--gapopen", gap,
      "--gapext", gap,
      "--match", "1",
      "--mismatch", "-1",
      "| awk '$1==\"H\" {print $9,$10,$4}'"
    ),
    intern = TRUE
  )
  stopifnot(attr(uc, "status") == 0)
  if (length(uc) > 0) {
    readr::read_delim(
      I(uc),
      col_names = c(if (id_is_int) "seq_idx" else "seq_id", "cluster", "dist"),
      delim = " ",
      col_types = if (id_is_int) "icd" else"ccd"
    )
  } else if (id_is_int) {
    tibble::tibble(
      seq_idx = integer(),
      cluster = character(),
      dist = numeric()
    )
  } else {
    tibble::tibble(
      seq_id = character(),
      cluster = character(),
      dist = numeric()
    )
  }
}

#' "uchime_ref" function of vsearch
#' @param query (`data.frame`, [`DNAStringSet`][Biostrings::XStringSet-class], `character` vector,
#' or file name) query sequences
#' @param ref (`data.frame`, [`DNAStringSet`][Biostrings::XStringSet-class], `character` vector, or
#' file name) reference sequences
#' @param ncpu (`integer` count) number of threads to use
#' @param id_only (`logical` flag) if `TRUE`, return only the sequence IDs
#' @param id_is_int (`logical` flag) if `TRUE`, return the sequence IDs as
#' integers
#' @return if `id_only` is FALSE, a `tibble::tibble` with columns `seq_id` (or
#' `seq_idx` if `id_is_int` is TRUE) and `seq`, where `seq_id` (`seq_idx`) is
#' the name of a sequence from `query` that is a chimera, and `seq` is the
#' sequence itself.  Alternatively, if `id_only` is TRUE, a `character` vector
#' (or `integer` vector if `id_is_int` is TRUE) of the sequence IDs which are
#' chimeras.
vsearch_uchime_ref <- function(query, ref, ncpu = local_cpus(), id_only = FALSE,
                               id_is_int = FALSE) {
  # avoid R CMD check NOTE for undeclared global variables due to NSE
  seq_id <- NULL

  # check arguments
  checkmate::assert_integerish(ncpu, lower = 1, min.len = 1, max.len = 1)
  checkmate::assert_flag(id_only)
  checkmate::assert_flag(id_is_int)

  if (checkmate::test_file_exists(query, "r")) {
    tquery <- query
  } else {
    tquery <- withr::local_tempfile(pattern = "query", fileext = ".fasta")
    write_sequence(query, tquery)
  }
  if (checkmate::test_file_exists(ref, "r")) {
    tref <- ref
  } else {
    tref <- withr::local_tempfile(pattern = "ref", fileext = ".fasta")
    write_sequence(ref, tref)
  }
  tchimeras <- withr::local_tempfile(pattern = "chimeras", fileext = ".fasta")
  vs <- system2(
    find_vsearch(),
    args = c(
      "--uchime_ref", tquery,
      "--db", tref,
      "--chimeras", tchimeras,
      "--threads", ncpu
    )
  )
  stopifnot(vs == 0L)
  if (id_only) {
    out <- names(Biostrings::fasta.seqlengths(tchimeras))
    if (id_is_int) {
      as.integer(out)
    } else {
      out
    }
  } else {
    out <- Biostrings::readDNAStringSet(tchimeras) |>
      as.character() |>
      tibble::enframe(name = "seq_id", value = "seq")
    if (id_is_int) {
      dplyr::transmute(out, seq_idx = as.integer(seq_id), seq)
    } else {
      out
    }
  }
}

#' Perform closed-reference clustering using vsearch
#' @param query (`data.frame`, [`Biostrings::DNAStringSet`][Biostrings::XStringSet-class], `character` vector,
#' or file name) query sequences
#' @param ref (`data.frame`, `Biostrings::DNAStringSet`[`Biostrings::DNAStringSet`][Biostrings::XStringSet-class], `character` vector, or
#' file name) reference sequences
#' @param threshold (`numeric` scalar) identity threshold, in range 0.0-1.0
#' @param ... additional arguments to pass to `vsearch_usearch_global()`
#' @return `tibble::tibble` with columns `seq_id` and `cluster`, where `seq_id`
#' is the name of a sequence from `query`, and `cluster` is the closest match
#' to that sequence in `ref`
#' @export
vsearch_usearch_global_closed_ref <- function(query, ref, threshold, ...) {
  # avoid R CMD check NOTE for undeclared global variables due to NSE
  seq_id <- cluster <- NULL
  out <- tibble::tibble(seq_id = character(0), cluster = character(0))
  while(sequence_size(query) > 0 && sequence_size(ref) > 0) {
    result <- vsearch_usearch_global(query, ref, threshold, ...)
    if (nrow(out) > 0) {
      result <- dplyr::left_join(
        result,
        out,
        by = c("cluster" = "seq_id"),
        suffix = c(".orig", "")
      ) |>
        dplyr::select(seq_id, cluster)
    }
    out <- dplyr::bind_rows(out, result)
    ref <- select_sequence(query, result$seq_id)
    query <- select_sequence(query, result$seq_id, negate = TRUE)
  }
  out
}

#' "cluster_smallmem" function of vsearch
#' @param seq (`data.frame`, [`DNAStringSet`][Biostrings::XStringSet-class], `character` vector, or
#' file name) sequences to cluster
#' @param threshold (`numeric` scalar) identity threshold, in range 0.0-1.0
#' @param ncpu (`integer` count) number of threads to use
#' @return `tibble::tibble` with columns `query` and `hit`, where `query` is
#' the name of a sequence from `seq`, and `hit` is the name of the sequence
#' which is the centroid of the cluster containing `query`
#' @export
vsearch_cluster_smallmem <- function(seq, threshold = 1, ncpu = local_cpus()) {
  if (is.character(seq) && length(seq) == 1 && file.exists(seq)) {
    tout <- seq
  } else {
    tout <- withr::local_tempfile(pattern = "data", fileext = ".fasta")
    write_sequence(seq, tout)
  }
  uc = system(
    paste(
      find_vsearch(),
      "--cluster_smallmem", tout,
      "--usersort",
      "--id", threshold,
      "--uc -",
      "--threads", ncpu,
      "| awk '$1==\"H\" {print $9,$10}'"
    ),
    intern = TRUE
  )
  stopifnot(attr(uc, "status") == 0)
  if (length(uc) > 0) {
    readr::read_delim(
      I(uc),
      col_names = c("query", "hit"),
      delim = " ",
      col_types = "cc"
    )
  } else {
    tibble::tibble(query = character(), hit = character())
  }
}

#' Collapse ASVs which are identical in their overlaps using VSEARCH
#' @param seqtab (integer `matrix` ) DADA2-style sequence table as returned by
#' `dada2::makeSequenceTable()`
#' @param ... additional arguments to pass to `vsearch_cluster_smallmem()`
#' @param ncpu (`integer` count) number of threads to use
#' @return `integer` matrix with identical sequences collapsed
collapseNoMismatch_vsearch <- function(seqtab, ..., ncpu = local_cpus()) {
  seqs <- colnames(seqtab)
  names(seqs) <- seq_along(seqs)
  matches <- vsearch_cluster_smallmem(seqs, ncpu = ncpu)
  map <- tibble::tibble(
    seq_idx_in = seq_len(ncol(seqtab)),
    seq_idx_out = seq_len(ncol(seqtab))
  )
  if (nrow(matches) > 0) {
    matches$query <- as.integer(matches$query)
    matches$hit <- as.integer(matches$hit)
    matches <- matches[order(matches$query),]
    for (i in unique(matches$hit)) {
      seqtab[,i] <- seqtab[,i] +
        as.integer(rowSums(seqtab[,matches$query[matches$hit == i], drop = FALSE]))
    }
    seqtab <- seqtab[,-matches$query]
    map$seq_idx_out[matches$query] <- matches$hit
    map$seq_idx_out = map$seq_idx_out - findInterval(map$seq_idx_out, matches$query)
  }
  attr(seqtab, "map") <- map
  return(seqtab)
}

#' Collapse ASVs which are identical in their overlaps using VSEARCH
#' @param seqtab (`data.frame`) long sequence table as returned by
#' `make_long_sequence_table()` or `make_mapped_sequence_table()`
#' @param seqs (`NULL`, `data.frame`,
#' [`DNAStringSet`][Biostrings::XStringSet-class], `character` vector, or
#' file name) if seqtab is a mapped sequence table, then the master sequence
#' list, otherwise should be `NULL`.
#' @param abund_col (`character` string) name of the column in `seqtab` which
#' contains the abundance of each sequence
#' @param fastx_index (`NULL` of file name) if `seqs` is a file name pointing to
#' a gzipped FASTA file, then the (optional) index file for that file
#' @param ... additional arguments to pass to `vsearch_cluster_smallmem()`
#' @param ncpu (`integer` count) number of threads to use
#' @return `data.frame` with columns `seq_idx_in` and `seq_idx_out`, where
#' `seq_idx_in` is the index of a sequence in `seqtab`, and `seq_idx_out` is
#' the index of the sequence which is the centroid of the cluster containing
#' `seq_idx_in`
#' @importFrom dplyr everything
nomismatch_hits_vsearch <- function(seqtab, seqs = NULL,
                                    abund_col = "nread",
                                    fastx_index = NULL,
                                    ...,
                                    ncpu = local_cpus()) {
  if (is.null(seqs)) {
    checkmate::assert_names(
      names(seqtab),
      must.include = c("seq", abund_col),
      disjunct.from = "seq_idx"
    )
    seqs <- sort_seq_table(seqtab, abund_col = abund_col, ...)
  } else {
    checkmate::assert(
      checkmate::check_names(
        names(seqtab),
        must.include = c("seq_idx", abund_col),
        disjunct.from = "seq"
      ),
      checkmate::check_names(
        names(seqtab),
        must.include = c("seq_id", abund_col),
        disjunct.from = "seq"
      )
    )
    o <- sort_seq_table(seqtab, seqs = seqs, abund_col = abund_col, ...)
    if (checkmate::check_file_exists(seqs)) {
      # no easy way to re-order without reading it all into memory
      seqs <- Biostrings::readDNAStringSet(seqs)
    }
    if (is.character(o)) o <- match(o, names(seqs))
    seqs <- seqs[o]
    names(seqs) <- as.character(o)
  }
  nseqs <- sequence_size(seqs)
  vsearch_cluster_smallmem(seqs, ncpu = ncpu) |>
    dplyr::mutate(dplyr::across(everything(), as.integer))
}
