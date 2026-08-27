#' Try to find the vsearch executable
#' @return (`character` string) the full path to the vsearch executable.
#' @export
find_vsearch <- function() {
  find_executable("vsearch")
}

#' "usearch_global" function of vsearch
#'
#' @param query (`data.frame`, [`DNAStringSet`][Biostrings::XStringSet-class],
#' `character` vector, or file name) query sequences
#' @param ref (`data.frame`, [`DNAStringSet`][Biostrings::XStringSet-class],
#' `character` vector, or file name) reference sequences
#' @param threshold (`numeric` scalar) identity threshold, in range 0.0-1.0
#' @param global (`logical` flag) if `TRUE`, end gaps and internal gaps are
#' penalized equally.  Otherwise end gaps are not penalized.
#' @param ncpu (`integer` count) number of threads to use
#' @param id_is_int (`logical` flag) if `TRUE`, return the sequence IDs as
#' integers
#' @param vsearch (`character` string) path to the vsearch executable
#'
#' @return `tibble::tibble` with columns `seq_id`, `clust`, and `dist`, where
#' `seq_id` is the name of a sequence from `query`, `clust` is the closest match
#' to that sequence in `ref`, and `dist` is the distance between them
#' @export
vsearch_usearch_global <- function(
  query,
  ref,
  threshold,
  global = TRUE,
  ncpu = local_cpus(),
  id_is_int = FALSE,
  vsearch = find_vsearch()
) {
  checkmate::check_flag(id_is_int)
  checkmate::assert_file_exists(vsearch, "x")
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
  # fmt: skip
  uc <- system(
    paste(
      vsearch,
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
      col_types = if (id_is_int) "icd" else "ccd"
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
#' @param query (`data.frame`, [`DNAStringSet`][Biostrings::XStringSet-class],
#' `character` vector, or file name) query sequences
#' @param ref (`data.frame`, [`DNAStringSet`][Biostrings::XStringSet-class],
#' `character` vector, or file name) reference sequences
#' @param ncpu (`integer` count) number of threads passed to vsearch
#'   (`--threads`).
#' @param id_only (`logical` flag) if `TRUE`, return only the sequence IDs
#' @param id_is_int (`logical` flag) if `TRUE`, return the sequence IDs as
#' integers
#' @param files (`character` vector) optional per-file paths overriding those
#' stored in the index, if `query` is a
#' [`fastqindexr_index`][fastqindexr::create_index()] object or `.fqi`
#' path(s). Useful after moving inputs or for `targets` dependency tracking.
#' @param seq_idx (`integer` vector) optional 1-based indices into the logical
#' sequence stream (`NULL` means all sequences in order). Applies after
#' concatenating multiple FASTA inputs, and supports duplicates and
#' reordering.
#' @param vsearch (`character`) path to the `vsearch` executable
#' @param ... currently unused; reserved for future extensions.
#' @return if `id_only` is FALSE, a `tibble::tibble` with columns `seq_id` (or
#' `seq_idx` if `id_is_int` is TRUE) and `seq`, where `seq_id` (`seq_idx`) is
#' the name of a sequence from `query` that is a chimera, and `seq` is the
#' sequence itself.  Alternatively, if `id_only` is TRUE, a `character` vector
#' (or `integer` vector if `id_is_int` is TRUE) of the sequence IDs which are
#' chimeras.
#' @export
vsearch_uchime_ref <- function(
  query,
  ref,
  ncpu = local_cpus(),
  id_only = FALSE,
  id_is_int = FALSE,
  files = NULL,
  seq_idx = NULL,
  vsearch = find_vsearch(),
  ...
) {
  # avoid R CMD check NOTE for undeclared global variables due to NSE
  seq_id <- NULL

  # check arguments
  checkmate::assert_integerish(ncpu, lower = 1, min.len = 1, max.len = 1)
  checkmate::assert_flag(id_only)
  checkmate::assert_flag(id_is_int)
  checkmate::assert_file_exists(vsearch, "x")

  if (is.list(seq_idx) && length(seq_idx) > 1L) {
    stop(
      "`seq_idx` must not be a list with more than one partition.",
      call. = FALSE
    )
  }

  indexed_like <- inherits(query, "fastqindexr_index") ||
    seq_batch_is_fqi_path_set(query)
  if (!is.null(files) && !indexed_like) {
    stop(
      "`files` is only valid when `query` is a fastqindexr_index or .fqi/.qs2 paths.",
      call. = FALSE
    )
  }
  tmp_parent <- environment()
  qfiles <- seq_batch_make_chunk_files(
    seqs = query,
    files = files,
    seq_idx = seq_idx,
    ncpu = 1L,
    local_envir = tmp_parent
  )
  if (length(qfiles) == 0L) {
    if (id_only) {
      if (id_is_int) {
        return(integer())
      }
      return(character())
    }
    if (id_is_int) {
      return(tibble::tibble(seq_idx = integer(), seq = character()))
    }
    return(tibble::tibble(seq_id = character(), seq = character()))
  }
  tquery <- qfiles[[1L]]
  if (checkmate::test_file_exists(ref, "r")) {
    tref <- ref
  } else {
    tref <- withr::local_tempfile(pattern = "ref", fileext = ".fasta")
    write_sequence(ref, tref)
  }
  tchimeras <- withr::local_tempfile(pattern = "chimeras", fileext = ".fasta")
  vs <- system2(
    vsearch,
    # fmt: skip
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
#' @param query (`data.frame`,
#' [`Biostrings::DNAStringSet`][Biostrings::XStringSet-class], `character`
#' vector, or file name) query sequences
#' @param ref (`data.frame`,
#' [`Biostrings::DNAStringSet`][Biostrings::XStringSet-class], `character`
#' vector, or file name) reference sequences
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
  while (sequence_size(query) > 0 && sequence_size(ref) > 0) {
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
#' @param seq (`data.frame`, [`DNAStringSet`][Biostrings::XStringSet-class],
#' `character` vector, or file name) sequences to cluster
#' @param threshold (`numeric` scalar) identity threshold, in range 0.0-1.0
#' @param ncpu (`integer` count) number of threads to use
#' @param vsearch (`character` string) path to the vsearch executable
#' @return `tibble::tibble` with columns `query` and `hit`, where `query` is
#' the name of a sequence from `seq`, and `hit` is the name of the sequence
#' which is the centroid of the cluster containing `query`
#' @export
vsearch_cluster_smallmem <- function(
  seq,
  threshold = 1,
  ncpu = local_cpus(),
  vsearch = find_vsearch()
) {
  checkmate::assert_file_exists(vsearch, "x")
  if (is.character(seq) && length(seq) == 1 && file.exists(seq)) {
    tout <- seq
  } else {
    tout <- withr::local_tempfile(pattern = "data", fileext = ".fasta")
    write_sequence(seq, tout)
  }
  # fmt: skip
  uc <- system(
    paste(
      vsearch,
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
#' @export
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
    matches <- matches[order(matches$query), ]
    for (i in unique(matches$hit)) {
      seqtab[, i] <- seqtab[, i] +
        as.integer(
          rowSums(
            seqtab[, matches$query[matches$hit == i], drop = FALSE]
          )
        )
    }
    seqtab <- seqtab[, -matches$query]
    map$seq_idx_out[matches$query] <- matches$hit
    map$seq_idx_out <- map$seq_idx_out -
      findInterval(map$seq_idx_out, matches$query)
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
#' list (plain or gzipped FASTA), otherwise should be `NULL`.
#' @param abund_col (`character` string) name of the column in `seqtab` which
#' contains the abundance of each sequence
#' @param fastx_index (`NULL` of file name) if `seqs` is a file name pointing to
#' a gzipped FASTA file, then the (optional) index file for that file
#' @param ... additional arguments to pass to `vsearch_cluster_smallmem()`
#' @param ncpu (`integer` count) number of threads to use
#' @return `data.frame` with columns `seq_idx_in` and `seq_idx_out`, where
#' `seq_idx_in` is the index of a sequence in `seqtab`, and `seq_idx_out` is
#' the index of the sequence which is the centroid of the cluster containing
#' `seq_idx_in`. Empty sequence inputs return an empty two-column mapping table.
#' @importFrom dplyr everything
nomismatch_hits_vsearch <- function(
  seqtab,
  seqs = NULL,
  abund_col = "nread",
  fastx_index = NULL,
  ...,
  ncpu = local_cpus()
) {
  empty_hits <- tibble::tibble(
    query = integer(),
    hit = integer()
  )
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
      if (length(Biostrings::fasta.seqlengths(seqs)) == 0L) {
        return(empty_hits)
      }
      # no easy way to re-order without reading it all into memory
      seqs <- Biostrings::readDNAStringSet(seqs)
    }
    if (is.character(o)) {
      o <- match(o, names(seqs))
    }
    seqs <- seqs[o]
    names(seqs) <- as.character(o)
  }
  if (length(seqs) == 0L) {
    return(empty_hits)
  }
  vsearch_cluster_smallmem(seqs, ncpu = ncpu) |>
    dplyr::mutate(dplyr::across(everything(), as.integer))
}

#' Perform taxonomic classification using SINTAX in VSEARCH
#' @param query (`data.frame`, [`DNAStringSet`][Biostrings::XStringSet-class],
#' `character` vector, or file name) query sequences. File inputs may be plain
#' FASTA or gzipped FASTA.
#' @param ref (`data.frame`, [`DNAStringSet`][Biostrings::XStringSet-class],
#' `character` vector, or file name) reference sequences
#' @param ncpu (`integer` count) number of threads passed to vsearch
#'   (`--threads`), or `NULL` to omit the flag.
#' @param id_is_int (`logical` flag) if `TRUE`, return the sequence IDs as
#' integers
#' @param hash (`character` string) hash value for the queries; ignored (but
#' used by `targets` for dependency tracking)
#' @param files (`character` vector) optional per-file paths overriding those
#' stored in the index, if `query` is a
#' [`fastqindexr_index`][fastqindexr::create_index()] object or `.fqi`
#' path(s). Useful after moving inputs or for `targets` dependency tracking.
#' @param seq_idx (`integer` vector) optional 1-based indices into the logical
#' sequence stream (`NULL` means all sequences in order). Applies after
#' concatenating multiple FASTA inputs, and supports duplicates and
#' reordering.
#' @param vsearch (`character`) path to the `vsearch` executable
#' @param ... currently unused; reserved for future extensions.
#' @return `tibble::tibble` with columns `seq_id` (or `seq_idx` if `id_is_int`
#' is TRUE), `rank`, `parent_taxonomy`, `taxon`, and `prob`, where `seq_id`
#' (`seq_idx`) is the ID of a sequence from `query`, `rank` is the taxonomic
#' rank of the taxon, `parent_taxonomy` is the parent taxonomic unit of the
#' taxon, `taxon` is the taxon name, and `prob` is the probability of the
#' taxon being the correct classification of the sequence. Empty query inputs
#' return an empty tibble with the same schema.
#' @export
sintax <- function(
  query,
  ref,
  ncpu = NULL,
  id_is_int = FALSE,
  hash = NULL,
  files = NULL,
  seq_idx = NULL,
  vsearch = find_vsearch(),
  ...
) {
  # avoid R CMD check NOTE about global variables due to NSE
  seq_id <- taxonomy <- taxon <- NULL
  empty_out <- tibble::tibble(
    seq_id = character(),
    rank = rank2factor(character()),
    parent_taxonomy = character(),
    taxon = character(),
    prob = numeric()
  )
  checkmate::assert_file_exists(ref, access = "r")
  checkmate::assert_count(ncpu, null.ok = TRUE)
  checkmate::assert_file_exists(vsearch, "x")
  if (is.list(seq_idx) && length(seq_idx) > 1L) {
    stop(
      "`seq_idx` must not be a list with more than one partition.",
      call. = FALSE
    )
  }
  indexed_like <- inherits(query, "fastqindexr_index") ||
    seq_batch_is_fqi_path_set(query)
  if (!is.null(files) && !indexed_like) {
    stop(
      "`files` is only valid when `query` is a fastqindexr_index or .fqi/.qs2 paths.",
      call. = FALSE
    )
  }
  tmp_parent <- environment()
  qfiles <- seq_batch_make_chunk_files(
    seqs = query,
    files = files,
    seq_idx = seq_idx,
    ncpu = 1L,
    local_envir = tmp_parent
  )
  if (length(qfiles) == 0L) {
    if (id_is_int) {
      return(tibble::add_column(empty_out, seq_idx = integer(), .before = 1))
    }
    return(empty_out)
  }
  tout <- qfiles[[1L]]
  if (length(Biostrings::fasta.seqlengths(tout)) == 0L) {
    if (id_is_int) {
      return(tibble::add_column(empty_out, seq_idx = integer(), .before = 1))
    }
    return(empty_out)
  }
  version <- processx::run(vsearch, "--version")$stderr |>
    sub("vsearch v([0-9.]+).+", "\\1", x = _) |>
    strsplit(split = ".", fixed = TRUE) |>
    unlist() |>
    as.integer()
  has_random <- version[1] > 2L || (version[1] == 2L && version[2] >= 28L)
  result <- processx::run(
    vsearch,
    # fmt: skip
    c(
      "--sintax", tout,
      if (has_random) "--sintax_random",
      "--db", ref,
      "--tabbedout", "-",
      if (!is.null(ncpu)) c("--threads", ncpu)
    )
  )
  stopifnot(result$status == 0)
  out <- (if (length(result$stdout) > 0) {
    suppressWarnings(
      readr::read_delim(
        I(result$stdout),
        col_names = c("seq_id", "taxonomy"),
        delim = "\t",
        col_types = "cc-"
      ),
      "vroom_parse_issue"
    )
  } else {
    tibble::tibble(seq_id = character(), taxonomy = character())
  }) |>
    tidyr::separate_longer_delim(taxonomy, delim = ",") |>
    tidyr::separate_wider_regex(
      taxonomy,
      patterns = c(
        "[dkpcofgst]:",
        taxon = ".+",
        "\\(",
        prob = "[0-9.]+",
        "\\)"
      )
    ) |>
    dplyr::mutate(
      rank = rank2factor(tax_ranks()[seq_len(dplyr::n())]),
      parent_taxonomy = purrr::accumulate(
        taxon,
        \(...) paste(..., sep = ",")
      ) |>
        dplyr::lag(default = ""),
      .by = seq_id,
      .before = taxon
    )
  if (id_is_int) {
    out <- dplyr::mutate(
      out,
      seq_idx = as.integer(seq_id),
      .keep = "unused",
      .before = 1
    )
  }
  out
}

#' Class to store merged-fastq filter options for VSEARCH/USEARCH
#'
#' In all cases, merged read pairs which do not meet the filter criteria are
#' excluded.
#'
#' Called with no arguments, this returns the values parsed from
#' `pipeline_options.yaml` when available, otherwise the constructor defaults.
#'
#' It appears that vsearch does not count "N"s towards expected errors, even
#' though their quality score looks like each one represents 0.63 expected
#' errors. Thus it is possible to allow sequences with N's while still having 0
#' maximum expected errors.
#'
#' @param maxEE (`numeric` scalar) maximum absolute number of expected errors
#'   allowed in merged read pairs.
#' @param maxEE_rate (`numeric` scalar) maximum number of expected errors
#'   allowed in merged read pairs, as a fraction of the sequence length.
#' @param maxNs (`integer` scalar) maximum number of "N"s allowed in merged
#'   read pairs.
#' @param maxLen (`integer` scalar) maximum length allowed for merged read
#'   pairs.
#' @param minLen (`integer` scalar) minimum length allowed for merged read
#'   pairs.
#' @return An object of class `merged_filter_options`.
#' @export
merged_filter_options <- function(
  maxEE = 1,
  maxEE_rate = NULL,
  maxNs = NULL,
  maxLen = NULL,
  minLen = NULL
) {
  if (
    missing(maxEE) &&
      missing(maxEE_rate) &&
      missing(maxNs) &&
      missing(maxLen) &&
      missing(minLen)
  ) {
    stored <- getOption("optimotu.pipeline.merged_filter_options")
    if (!is.null(stored)) {
      return(stored)
    }
  }
  checkmate::assert_number(maxEE, lower = 0, null.ok = TRUE)
  checkmate::assert_number(maxEE_rate, lower = 0, upper = 1, null.ok = TRUE)
  checkmate::assert_count(maxNs, null.ok = TRUE)
  checkmate::assert_count(maxLen, positive = TRUE, null.ok = TRUE)
  checkmate::assert_count(minLen, null.ok = TRUE)
  structure(
    list(
      maxEE = maxEE,
      maxEE_rate = maxEE_rate,
      maxNs = maxNs,
      maxLen = maxLen,
      minLen = minLen
    ),
    class = "merged_filter_options"
  )
}

#' Names of options in merged_filter_options
#' @export
merged_filter_option_names <- c(
  "maxEE",
  "maxEE_rate",
  "maxNs",
  "maxLen",
  "minLen"
)

#' Update method for merged filter options
#' @param object ([`merged_filter_options`][merged_filter_options()]) existing
#' merged filter options object to modify
#' @param new_options (named `list`, single-row `data.frame`, named
#'  `character` vector, or named `numeric` vector) new values for the options.
#'  If a `data.frame`, then the values in each column should be all the same.
#' @param ... Additional arguments (ignored)
#' @exportS3Method stats::update
update.merged_filter_options <- function(object, new_options, ...) {
  checkmate::assert(
    checkmate::check_list(new_options, null.ok = TRUE),
    checkmate::check_data_frame(new_options, null.ok = TRUE),
    checkmate::check_character(new_options, null.ok = TRUE),
    checkmate::check_numeric(new_options, null.ok = TRUE)
  )
  if (is.null(new_options) || length(new_options) == 0L) {
    return(object)
  }
  for (nm in intersect(names(new_options), merged_filter_option_names)) {
    object[[nm]] <- new_options[[nm]]
  }
  do.call(merged_filter_options, object)
}

# Write an empty FASTA/FASTQ file, gzipped when the path ends in .gz.
write_empty_seq_file <- function(path) {
  dirs <- unique(dirname(path))
  dirs <- dirs[dirs != "."]
  if (length(dirs) > 0L) {
    dir.create(dirs, recursive = TRUE, showWarnings = FALSE)
  }
  for (p in path) {
    if (endsWith(p, ".gz")) {
      con <- gzfile(p, "wb")
      close(con)
    } else {
      file.create(p)
    }
  }
  path
}

empty_uc_cluster <- function() {
  structure(
    list(
      clusters = tibble::tibble(
        clust_idx = integer(),
        size = integer(),
        seq = character()
      ),
      map = tibble::tibble(
        clust_idx = integer(),
        seq_id = character()
      )
    ),
    class = "uc_cluster"
  )
}

# Start a process/pipeline for a single fastq_mergepairs call.
# Returns the process to wait on (gzip when compressing, otherwise vsearch).
fastq_merge_pairs_process <- function(vsearch, args, seq_out, compress) {
  if (isTRUE(compress)) {
    processx::pipeline$new(
      cmds = list(
        c(vsearch, unlist(as.character(args))),
        c("gzip", "-cf")
      ),
      stdout = seq_out,
      stderr = "|"
    )$get_processes()[[2]]
  } else {
    processx::process$new(
      command = vsearch,
      args = as.character(args),
      stderr = "|",
      poll_connection = TRUE
    )
  }
}

merge_pairs_finish <- function(proc) {
  # Drain stderr so a full pipe cannot stall the process.
  if (proc$has_error_connection()) {
    try(proc$read_error(), silent = TRUE)
  }
  if (proc$is_alive()) {
    return(FALSE)
  }
  status <- proc$get_exit_status()
  if (!identical(status, 0L)) {
    err <- tryCatch(proc$read_all_error(), error = function(e) "")
    stop(
      "vsearch/USEARCH fastq_mergepairs failed with exit status ",
      status,
      if (nzchar(err)) paste0(":\n", err) else ".",
      call. = FALSE
    )
  }
  TRUE
}

#' Assemble Illumina read pairs using USEARCH or VSEARCH
#' @param seq_R1 (`character` vector) path(s) of one or more FASTQ files,
#'   possibly gzipped, to be assembled in the forward orientation.
#' @param seq_R2 (`character` vector) path(s) of one or more FASTQ files,
#'   possibly gzipped, to be assembled in the reverse complement orientation.
#'   Must be the same number of files as `seq_R1`, and each file must have the
#'   same number of reads.
#' @param min_overlap (`integer` scalar) minimum length of the overlapping
#'   region in order to assemble a read pair.
#' @param seq_out (`character` vector) file names
#' @param max_mismatch (`numeric` scalar) if strictly < 1, the fraction of bases
#'   in the overlapping region which are allowed to be mismatches (
#'   `--fastq_maxdiffpct` argument, but should be expressed as a fraction rather
#'   than percent). If >= 1, the number of bases in the overlapping region which
#'   are allowed to be mismatches (`--fastq_maxdiffs` argument). In the latter
#'   case it should be integer-valued. Default: 10
#' @param threads (`integer` scalar) number of threads to use for *each*
#'   usearch/vsearch process. Individual processes eventually become I/O bound
#'   and more threads may provide little improvement, or even slow the total
#'   processing time. The point of diminishing returns is probably dependent on
#'   the system and the input size.
#' @param shards (`integer` scalar) number of parallel usearch/vsearch processes
#'   to run at once.  The total number of CPU threads used is something like
#'   `threads * shards` (although the R process itself will also consume some).
#' @param filter_options (result object from [`merged_filter_options()`]) additional
#'   options to filter the results; passed to USEARCH or VSEARCH.
#' @param compress (`logical` flag) if `TRUE`, then vsearch/usearch output is
#'   piped through `gzip` to compress it. The default autodetects based on
#'   whether the filenames in `seq_out` end with `".gz"`.
#' @param fastq (`logical` flag) if `TRUE`, then output is in FASTQ format;
#'   otherwise it is FASTA. By default this is autodetected from the filenames
#'   in `seq_out`.
#' @param vsearch (`character` string) path to USEARCH or VSEARCH executable.
#'
#' @return `character` vector giving the output file names (as given in
#' `seq_out`)
#' @export
vsearch_fastq_merge_pairs <- function(
  seq_R1,
  seq_R2,
  seq_out,
  min_overlap = 5,
  max_mismatch = 10,
  threads = 1,
  shards = min(local_cpus() %/% threads, length(seq_out)),
  filter_options = merged_filter_options(),
  compress = all(endsWith(seq_out, ".gz")),
  fastq = all(grepl("[.]fa?s?t?q", basename(seq_out))),
  vsearch = find_vsearch()
) {
  checkmate::assert_character(seq_R1)
  checkmate::assert_character(seq_R2, len = length(seq_R1))
  checkmate::assert_character(seq_out, len = length(seq_R1))
  if (length(seq_R1) == 0L) {
    return(character())
  }
  checkmate::assert_file_exists(seq_R1, access = "r")
  checkmate::assert_file_exists(seq_R2, access = "r")
  checkmate::assert_path_for_output(seq_out, overwrite = TRUE)
  checkmate::assert_count(min_overlap)
  checkmate::assert_integerish(min_overlap, lower = 5L)
  checkmate::assert_number(max_mismatch, lower = 0)
  checkmate::assert_count(threads, positive = TRUE)
  checkmate::assert_count(shards, positive = TRUE)
  checkmate::assert_class(filter_options, "merged_filter_options")
  checkmate::assert_flag(compress)
  checkmate::assert_flag(fastq)
  checkmate::assert_file_exists(vsearch, "x")
  write_empty_seq_file(seq_out)
  nonempty <- which(sequence_size(seq_R1) > 0L & sequence_size(seq_R2) > 0L)
  if (length(nonempty) == 0L) {
    return(seq_out)
  }
  seq_R1 <- seq_R1[nonempty]
  seq_R2 <- seq_R2[nonempty]
  seq_todo <- seq_out[nonempty]
  shards <- max(1L, min(as.integer(shards), length(seq_todo)))

  args <- c(
    # fmt: skip
    list(
      "--fastq_mergepairs", seq_R1,
      "--reverse", seq_R2,
      "--fastq_minovlen", min_overlap,
      "--threads", threads
    ),
    if (isTRUE(fastq)) {
      list("--fastqout")
    } else {
      list("--fastaout")
    },
    if (isTRUE(compress)) {
      # vsearch treats "-" as stdout; "--" is a literal output filename.
      list("-")
    } else {
      list(seq_todo)
    },
    if (max_mismatch >= 1) {
      list("--fastq_maxdiffs", max_mismatch)
    } else {
      list("--fastq_maxdiffpct", max_mismatch * 100)
    },
    if (!is.null(filter_options$maxEE)) {
      list("--fastq_maxee", filter_options$maxEE)
    },
    if (!is.null(filter_options$maxEE_rate)) {
      list("--fastq_maxee_rate", filter_options$maxEE_rate)
    },
    if (!is.null(filter_options$maxNs)) {
      list("--fastq_maxns", filter_options$maxNs)
    },
    if (!is.null(filter_options$maxLen)) {
      list("--fastq_maxmergelen", filter_options$maxLen)
    },
    if (!is.null(filter_options$minLen)) {
      list("--fastq_minmergelen", filter_options$minLen)
    }
  )

  args <- as.data.frame(args, stringsAsFactors = FALSE)
  processes <- vector("list", shards)
  i <- 0
  n_finished <- 0
  n_todo <- length(seq_todo)
  while (i < shards && i < n_todo) {
    i <- i + 1
    processes[[i]] <-
      fastq_merge_pairs_process(vsearch, args[i, ], seq_todo[i], compress)
  }
  while (n_finished < n_todo) {
    processx::poll(processes, -1)
    for (j in rev(seq_along(processes))) {
      if (merge_pairs_finish(processes[[j]])) {
        n_finished <- n_finished + 1
        if (i < n_todo) {
          i <- i + 1
          processes[[j]] <-
            fastq_merge_pairs_process(
              vsearch,
              args[i, ],
              seq_todo[i],
              compress
            )
        } else {
          processes <- processes[-j]
        }
      }
    }
  }
  seq_out
}

# start the process(es) for dereplicating and denoising a single file
unoise_process <- function(vsearch, seq_in, args) {
  processx::pipeline$new(
    cmds = list(
      c(
        vsearch,
        "--fastx_uniques",
        seq_in,
        "--sizeout",
        "--fastaout",
        "-"
      ),
      c(vsearch, unlist(as.character(args)))
    ),
    stdout = "|",
    stderr = "|"
  )$get_processes()[[2]]
}

#' Find denoised sequences using UNOISE3 (`vsearch --cluster_unoise`)
#'
#' This command actually executes both `--fastx_uniques` (or equivalent
#' dereplication) and `--cluster_unoise`, since the former is required for the
#' latter. The function name retains "unoise2" for historical reasons.
#'
#' @param seq (`character` vector) path(s) of one or more FASTA or FASTQ files,
#'   possibly gzipped
#' @param min_size (`integer` scalar) minimum abundance of a sequence to be
#'   considered a centroid (`--minsize`)
#' @param alpha (`numeric` scalar) alpha parameter for UNOISE (`--unoise_alpha`)
#' @param threads (`integer` scalar) number of threads to use per vsearch
#' process
#' @param shards (`integer` scalar) number of shards (parallel vsearch calls) to
#' use
#' @param vsearch (`character` scalar) path to vsearch executable
#' @return a named `list` of objects of class `uc_cluster`, where names match
#'   the values of `seq`. Each `uc_cluster` has two members, each a
#'   `data.frame`:
#'   - `clusters` with columns:
#'     - `clust_idx` - cluster index; integer from 0 to number of clusters - 1
#'     - `size` - the number of reads in the cluster
#'     - `seq` - the representative sequence (i.e. centroid) of the cluster
#'   - `map` with columns:
#'     - `clust_idx` - cluster index; integer from 0 to number of clusters - 1
#'     - `seq_id` - the sequence ID of a sequence which is a member of the
#'       cluster. Although input sequences are typically dereplicated and
#'       include a `";size="` annotation, it is stripped off so that the ID
#'       should match a read in the original file.
#'
#'   Empty inputs return an empty `uc_cluster` without invoking vsearch.
#' @export
vsearch_cluster_unoise2 <- function(
  seq,
  min_size = 8,
  alpha = 2.0,
  threads = 1,
  shards = min(local_cpus() %/% threads, length(seq)),
  vsearch = find_vsearch()
) {
  # avoid R CMD check NOTE about global variables due to NSE
  type <- clust_idx <- size <- seq_id <- NULL
  checkmate::assert_character(seq)
  if (length(seq) == 0L) {
    return(stats::setNames(list(), character()))
  }
  checkmate::assert_file_exists(seq, "r")
  checkmate::assert_count(min_size, positive = TRUE)
  checkmate::assert_number(alpha, lower = 0, finite = TRUE)
  checkmate::assert_count(threads, positive = TRUE)
  checkmate::assert_file_exists(vsearch, "x")
  # fmt: skip
  args <- c(
    "--cluster_unoise", "-",
    "--sizein",
    "--minsize", as.character(min_size),
    "--unoise_alpha", as.character(alpha),
    "--threads", as.character(threads),
    "--uc", "-",
    "--strand", "plus"
  )

  result <- vector("list", length(seq))
  names(result) <- seq
  nonempty <- which(sequence_size(seq) > 0L)
  empty <- setdiff(seq_along(seq), nonempty)
  result[empty] <- replicate(
    length(empty),
    empty_uc_cluster(),
    simplify = FALSE
  )
  if (length(nonempty) == 0L) {
    return(result)
  }
  seq_todo <- seq[nonempty]
  shards <- max(1L, min(shards, length(seq_todo)))
  checkmate::assert_count(shards, positive = TRUE)

  processes <- vector("list", shards)
  index <- integer(shards)
  i <- 0
  n_finished <- 0
  n_todo <- length(seq_todo)
  while (i < shards && i < n_todo) {
    i <- i + 1
    processes[[i]] <-
      unoise_process(vsearch, seq_todo[i], args)
    index[i] <- nonempty[i]
  }
  while (n_finished < n_todo) {
    poll_result <- processx::poll(processes, -1)
    for (j in rev(seq_along(poll_result))) {
      if (poll_result[[j]]["output"] == "ready") {
        if (processes[[j]]$is_incomplete_output()) {
          o <- processes[[j]]$read_output_lines()
          if (length(o) > 0) {
            o <- readr::read_tsv(
              I(o),
              # fmt: skip
              col_names = c("type", "clust_idx", "size", "sim", "strand", "na1",
                "na2", "na3", "seq_id", "hit_id"),
              col_types = "cii-----c-",
              na = c("NA", "*")
            ) |>
              dplyr::mutate(seq_id = sub(";size=[0-9]+", "", seq_id))
            result[[index[j]]] <- dplyr::bind_rows(result[[index[j]]], o)
          }
        } else {
          n_finished <- n_finished + 1
          if (processes[[j]]$get_exit_status() != 0) {
            stop(
              "vsearch/USEARCH process failed with exit status ",
              processes[[j]]$get_exit_status(),
              call. = FALSE
            )
          }
          if (is.null(result[[index[j]]])) {
            result[[index[j]]] <- empty_uc_cluster()
          } else {
            r <- result[[index[j]]] # just an alias for brevity
            r <- structure(
              list(
                clusters = dplyr::filter(r, type == "C") |>
                  dplyr::select(clust_idx, size, seq = seq_id),
                map = dplyr::filter(r, type != "C") |>
                  dplyr::select(clust_idx, seq_id)
              ),
              class = "uc_cluster"
            )
            # Replace sequence identifiers with sequences themselves.
            # There is no convenient vsearch/usearch output that provides both
            # the mapping and the centroid sequences, so this requires an extra
            # file read.  Seems easier to extract from the original file than
            # to have vsearch write them to a temp file.
            orig_seqs <- if (grepl(fastq_regex, seq[index[j]])) {
              Biostrings::readDNAStringSet(seq[index[j]], format = "fastq")
            } else {
              Biostrings::readDNAStringSet(seq[index[j]], format = "fasta")
            }
            # Vsearch strips off everything after the first whitespace in the
            # sequence ID, so we need to do the same here to match.
            names(orig_seqs) <- sub("\\s.*", "", names(orig_seqs))
            r$clusters$seq <- as.character(orig_seqs[r$clusters$seq])
            result[[index[j]]] <- r
          }
          if (i < n_todo) {
            i <- i + 1
            processes[[j]] <-
              unoise_process(vsearch, seq_todo[i], args)
            index[j] <- nonempty[i]
          } else {
            processes <- processes[-j]
            index <- index[-j]
          }
        }
      }
    }
  }
  result
}
