
#' Drop columns from a DADA2-style sequence table
#'
#' @param seqtable (`integer` matrix) a DADA2-style sequence table
#' @param which (`integer` or `character` vector) column indices to drop
#' @return a new sequence table with the specified columns removed
#' @export
drop_from_seqtable <- function(seqtable, which) {
  if (is.character(which)) which <- as.integer(which)
  checkmate::assert_integerish(which, lower = 1, upper = ncol(seqtable))
  if (length(which) == 0) {
    seqtable
  } else {
    seqtable[,-which,drop = FALSE]
  }
}

#### long sequence table ####

#' Convert an object into a long (i.e. sparse) sequence occurrence table
#'
#' @param x (`data.frame` as returned by `dada2::mergePairs()`, or integer matrix
#'  as returned by `dada2::makeSequenceTable()`, or a list of one of these.)
#' @param rc (logical flag) if TRUE, sequences in `x` will be reverse complemented.
#'
#' @return a `data.frame` with columns `sample`, `seq`, and `nread`

make_long_sequence_table <- function(x, rc = FALSE) {
  UseMethod("make_long_sequence_table", x)
}

#' @rdname make_long_sequence_table
#' @exportS3Method
make_long_sequence_table.data.frame <- function(x, rc = FALSE) {
  checkmate::assert_data_frame(x, col.names = "named")
  checkmate::assert_names(names(x), must.include = c("sequence", "abundance"))
  checkmate::assert_flag(rc)
  if ("accept" %in% names(x)) {
    checkmate::assert_logical(x$accept)
    x <- x[x$accept,]
  }
  out <- x[c("sequence", "abundance")]
  names(out) <- c("seq", "nread")
  if (isTRUE(rc)) out$seq <- dada2::rc(out$seq)
  out
}

#' @rdname make_long_sequence_table
#' @exportS3Method
make_long_sequence_table.matrix <- function(x, rc = FALSE) {
  checkmate::assert_integerish(x)
  checkmate::assert_flag(rc)
  if (isTRUE(rc)) colnames(x) <- dada2::rc(colnames(x))
  if (typeof(x) != "integer") mode(x) <- "integer"
  x[x==0L] <- NA_integer_
  as.data.frame(x) |>
    tibble::rownames_to_column("sample") |>
    tidyr::pivot_longer(-1, names_to = "seq", values_to = "nread", values_drop_na = TRUE)
}

#' @rdname make_long_sequence_table
#' @exportS3Method
make_long_sequence_table.list <- function(x, rc = FALSE) {
  # avoid R CMD check note for undefined global variables due to NSE
  nread <- NULL

  out <- if (checkmate::test_list(x, types = "data.frame")) {
    checkmate::assert_named(x)
    purrr::map_dfr(x, make_long_sequence_table.data.frame, rc = rc, .id = "sample")
  } else if (checkmate::test_list(x, types = "matrix")) {
    purrr::map_dfr(x, make_long_sequence_table.matrix, rc = rc)
  } else {
    stop("cannot determine entry type in make_long_sequence_table.list")
  }
  dplyr::summarize(out, nread = sum(nread), .by = c(sample, seq))
}

#' Convert an object into a long (i.e. sparse) sequence occurrence table where
#' sequences are stored as integer indices to a master list
#'
#' @param x (`data.frame` as returned by `dada2::mergePairs()`, or integer
#' matrix as returned by `dada2::makeSequenceTable()`, or a list of one of these.)
#' @param seqs (`character` vector, file name, or `Biostrings::XStringSet`)
#' master list of sequences
#' @param rc (logical flag) if TRUE, sequences in `x` will be reverse complemented.
#'
#' @return a `data.frame` with columns `sample`, `seq_idx`, and `nread`
#' @export

make_mapped_sequence_table <- function(x, seqs, rc = FALSE) {
  UseMethod("make_mapped_sequence_table", x)
}

#' @rdname make_mapped_sequence_table
#' @exportS3Method
make_mapped_sequence_table.data.frame <- function(x, seqs, rc = FALSE) {
  checkmate::assert_data_frame(x, col.names = "named")
  checkmate::assert_names(names(x), must.include = c("sequence", "abundance"))
  checkmate::assert_flag(rc)
  if ("accept" %in% names(x)) {
    checkmate::assert_logical(x$accept)
    x <- x[x$accept,]
  }
  if (checkmate::test_file_exists(seqs, "r")) {
    seqs <- Biostrings::readDNAStringSet(seqs)
  }
  out <- x[c("sequence", "abundance")]
  names(out) <- c("seq_idx", "nread")
  if (isTRUE(rc)) {
    out$seq_idx <- BiocGenerics::match(dada2::rc(out$seq_idx), seqs)
  } else {
    out$seq_idx <- BiocGenerics::match(out$seq_idx, seqs)
  }
  out
}

#' @rdname make_mapped_sequence_table
#' @exportS3Method
make_mapped_sequence_table.matrix <- function(x, seqs, rc = FALSE) {
  checkmate::assert_integerish(x)
  checkmate::assert_flag(rc)
  if (isTRUE(rc)) colnames(x) <- dada2::rc(colnames(x))
  if (checkmate::test_file_exists(seqs, "r")) {
    seqs <- Biostrings::readDNAStringSet(seqs)
  }
  colnames(x) <- BiocGenerics::match(colnames(x), seqs)
  if (typeof(x) != "integer") mode(x) <- "integer"
  x[x==0L] <- NA_integer_
  as.data.frame(x) |>
    tibble::rownames_to_column("sample") |>
    tidyr::pivot_longer(
      -1,
      names_to = "seq_idx",
      names_transform = as.integer,
      values_to = "nread",
      values_drop_na = TRUE
    )
}

#' @rdname make_mapped_sequence_table
#' @exportS3Method
make_mapped_sequence_table.list <- function(x, seqs, rc = FALSE) {
  # avoid R CMD check note for undefined global variables due to NSE
  nread <- seq_idx <- NULL

  if (checkmate::test_file_exists(seqs, "r")) {
    seqs <- Biostrings::readDNAStringSet(seqs)
  }
  out <- if (checkmate::test_list(x, types = "data.frame")) {
    checkmate::assert_named(x)
    if (length(x) == 0) {
      tibble::tibble(sample = character(), seq_idx = integer(), nread = integer())
    } else {
      purrr::map_dfr(x, make_mapped_sequence_table.data.frame, seqs = seqs, rc = rc, .id = "sample")
    }
  } else if (checkmate::test_list(x, types = "matrix")) {
    purrr::map_dfr(x, make_mapped_sequence_table.matrix, seqs = seqs, rc = rc)
  } else {
    stop("cannot determine entry type in make_long_sequence_table.list")
  }
  dplyr::summarize(out, nread = sum(nread), .by = c(sample, seq_idx))
}

#### deduplication ####

#' Remove "no-mismatch" duplicates from a mapped long sequence table
#' @param seqtable (`data.frame`) a mapped long sequence table
#' @param hits (`data.frame`) a table of hits from, e.g.,
#' `vsearch_usearch_global()`. Should have columns `query` and `hit`, both of
#' which contain integer indices which match the `seq_idx` column of `seqtable`.
#' @param abund_col (`character`) column name containing abundance information
#' @param sample_cols (`character`) column name(s) containing sample identifiers
#' @param merge (`logical`) if TRUE, merge duplicate sequences into a single row
#' per sample
#' @return a deduplicated sequence table
#' @export
#' @importFrom dplyr all_of any_of
deduplicate_seqtable <- function(seqtable, hits, abund_col = "nread",
                                 sample_cols = "sample", merge = TRUE) {
  # avoid R CMD check note for undefined global variables due to NSE
  hit <- query <- seq_idx <- NULL
  checkmate::assert_character(abund_col)
  checkmate::assert_character(sample_cols)
  checkmate::assert_data_frame(seqtable)
  checkmate::assert_names(
    names(seqtable),
    must.include = c("seq_idx", abund_col, sample_cols)
  )
  checkmate::assert_data_frame(hits)
  checkmate::assert_names(
    names(hits),
    must.include = c("query", "hit")
  )
  checkmate::assert_integer(hits$query, lower = 0L, any.missing = FALSE)
  checkmate::assert_integer(hits$hit, lower = 0L, any.missing = FALSE)
  checkmate::assert_integer(seqtable$seq_idx, lower = 0L, any.missing = FALSE)
  checkmate::assert_flag(merge)

  hits <- dplyr::arrange(hits, query)
  seqtable <- dplyr::left_join(
    seqtable,
    hits,
    by = c("seq_idx" = "query")
  ) |>
    dplyr::mutate(
      seq_idx = dplyr::coalesce(hit, seq_idx),
      .keep = "unused"
    )
  if (isTRUE(merge)) {
    seqtable <- dplyr::summarize(
      seqtable,
      dplyr::across(all_of(abund_col), sum),
      .by = any_of(c("seq_idx", sample_cols))
    )
  }
  seqtable$seq_idx <- seqtable$seq_idx - findInterval(seqtable$seq_idx, hits$query)
  seqtable
}

# TODO: allow sequence list which is not a file, not named with integers, etc.
# TODO: do this without reading the whole file into memory
#' Remove "no-mismatch" duplicates from a sequence file
#' @param seqs (`character`) file name of a sequence file
#' @param hits (`data.frame`) a table of hits from, e.g.,
#' `vsearch_usearch_global()`. Should have columns `query` and `hit`, both of
#' contain integer indices for sequences in `seqs`.
#' @param outfile (`character`) file name to write deduplicated sequences to
#' @return the file name of the deduplicated sequences
#' @export
deduplicate_seqs <- function(seqs, hits, outfile) {
  checkmate::assert_file_exists(seqs, "r")
  checkmate::assert_data_frame(hits)
  checkmate::assert_names(
    names(hits),
    must.include = c("query", "hit")
  )
  checkmate::assert_integer(hits$query, lower = 0L, any.missing = FALSE)
  checkmate::assert_integer(hits$hit, lower = 0L, any.missing = FALSE)
  if (nrow(hits) > 0) {
    out <- Biostrings::readBStringSet(seqs)[-hits$query]
    names(out) <- as.character(seq_along(out))
    write_sequence(out, outfile, compress = endsWith(outfile, ".gz"), compression_level = 9)
  } else {
    file.copy(seqs, outfile, overwrite = TRUE)
    outfile
  }
}

#' Deduplicate a sequence index vector
#' @param seq_idx (`integer`) a vector of sequence indices
#' @param hits (`data.frame`) a table of hits from, e.g.,
#' `vsearch_usearch_global()`. Should have columns `query` and `hit`, both
#' containing integers appearing in `seq_idx`.
#' @param merge (`logical`) if TRUE, merge duplicate sequence indices into a
#' single value
#' @return a deduplicated sequence index vector; the same length as `seq_idx`
#' if `merge` is FALSE, otherwise (possibly) shorter.
#' @export
deduplicate_seq_idx <- function(seq_idx, hits, merge = TRUE) {
  checkmate::assert_integerish(seq_idx, lower = 1)
  deduplicate_seqtable(
    seqtable = tibble::tibble(
      seq_idx = seq_idx
    ),
    hits = hits,
    sample_cols = character(),
    abund_col = character(),
    merge = merge
  )$seq_idx
}


#' Sort a sequence table
#' @param seqtable (`matrix` or `data.frame`) a sequence table as returned by
#' `dada2::makeSequenceTable()`, `make_long_sequence_table()`, or
#' `make_mapped_sequence_table()`
#' @param ... additional arguments passed to methods, currently ignored
#' @return for the `matrix` method, a sorted sequence table; for the `data.frame`
#' method, appropriate identifiers for the sequences, sorted.  These will be
#' either the sequences themselves, `seq_id`, or `seq_idx`, depending on the
#' input.
#' @export
sort_seq_table <- function(seqtable, ...) {
  UseMethod("sort_seq_table", seqtable)
}

#' @rdname sort_seq_table
#' @exportS3Method sort_seq_table matrix
sort_seq_table.matrix <- function(seqtable, ...) {
  # avoid R CMD check note for undefined global variables due to NSE
  var <- seq_id_out <- NULL

  colorder <- order(
    -colSums(seqtable > 0), # prevalence, highest to lowest
    -colSums(seqtable), # abundance, highest to lowest
    -apply(seqtable, 2, var), # variance, highest to lowest
    seqhash(colnames(seqtable)) # hash of sequence (pseudorandom but stable)
  )
  if (is.null(attr(seqtable, "map"))) {
    seqtable[order(rownames(seqtable)), colorder]
  } else {
    structure(
      seqtable[order(rownames(seqtable)), colorder],
      map = dplyr::mutate(attr(seqtable, "map"), seq_id_out = order(colorder)[seq_id_out])
    )
  }
}

#' @rdname sort_seq_table
#' @exportS3Method sort_seq_table data.frame
#' @param seqs (`character`, [`DNAStringSet`][`Biostrings::XStringSet-class`],
#' file name, or `data.frame`) sequences, for mapped sequence tables (currently
#' unused)
#' @param abund_col (`character`) column name containing abundance information
sort_seq_table.data.frame <- function(seqtable, seqs = NULL, abund_col = "nread", ...) {
  # avoid R CMD check note for undefined global variables due to NSE
  prevalence <- var <- NULL

  # TODO: add some verification here
  abund <- as.symbol(abund_col)
  seqorder <- dplyr::summarize(
    seqtable,
    prevalence = dplyr::n(),
    abundance = sum(!!abund),
    variance = if (prevalence == 1) 0 else var(!!abund),
    .by = any_of(c("seq", "seq_idx", "seq_id"))
  )
  seqorder$hash <-
    if ("seq" %in% names(seqorder)) {
      seqhash(seqorder$seq)
    } else if ("seq_idx" %in% names(seqorder)) {
      hash_sequences(seqs, use_names = FALSE)[as.integer(seqorder$seq_idx)]
    } else if ("seq_id" %in% names(seqorder)) {
      unname(hash_sequences(seqs, use_names = TRUE)[seqorder$seq_id])
    }
  out <- order(
    -seqorder$prevalence,
    -seqorder$abundance,
    -seqorder$variance,
    seqorder$hash
  )
  if ("seq" %in% names(seqorder)) {
    seqorder$seq[out]
  } else if ("seq_idx" %in% names(seqorder)) {
    out
  } else if ("seq_id" %in% names(seqorder)) {
    seqorder$seq_id[out]
  }
}
