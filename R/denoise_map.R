# SPDX-FileCopyrightText: 2026, Brendan Furneaux
# SPDX-License-Identifier: MIT

#' Map denoiser-local ASV ids to a master sequence list
#'
#' One row per denoiser ASV (per sample when `x` is a list). Matching against a
#' large `seq_all` is done once per call; pass the result to
#' [denoise_map_to_seqtable()] and to [dada2_read_map()] /
#' [unoise_read_map()].
#'
#' Rows with unmatched `seq_idx` (`NA`) are retained so fate maps can
#' distinguish "denoised but missing from `seq_all`" from "never denoised".
#'
#' @param x (`data.frame` as returned by `dada2::mergePairs()`, integer
#'   matrix as returned by `dada2::makeSequenceTable()`, `uc_cluster` as
#'   returned by [vsearch_cluster_unoise2()], or a named list of one of these)
#' @param seqs (`character` vector, file name, or `Biostrings::XStringSet`)
#'   master list of sequences
#' @param rc (`logical`) if `TRUE`, sequences in `x` are reverse-complemented
#'   before matching
#'
#' @return `data.frame` with columns:
#'   - `sample` (character; present for list and matrix methods)
#'   - `denoise_idx` (integer) denoiser-local ASV id: row index in
#'     `mergePairs` output, column index in a sequence table, or
#'     `uc_cluster$clusters$clust_idx` (0-based)
#'   - `seq_idx` (integer) index into `seqs`; `NA` if unmatched
#'   - `nread` (integer) abundance
#' @seealso [denoise_map_to_seqtable()], [make_mapped_sequence_table()]
#' @export
make_denoise_map <- function(x, seqs, rc = FALSE) {
  UseMethod("make_denoise_map", x)
}

#' @rdname make_denoise_map
#' @exportS3Method
make_denoise_map.data.frame <- function(x, seqs, rc = FALSE) {
  checkmate::assert_data_frame(x, col.names = "named")
  checkmate::assert_names(names(x), must.include = c("sequence", "abundance"))
  checkmate::assert_flag(rc)
  # Row index before any accept filter so it remains a valid index into
  # the original mergePairs table used by dada2_read_map().
  denoise_idx <- seq_len(nrow(x))
  if ("accept" %in% names(x)) {
    checkmate::assert_logical(x$accept)
    keep <- which(x$accept)
    x <- x[keep, , drop = FALSE]
    denoise_idx <- denoise_idx[keep]
  }
  tibble::tibble(
    denoise_idx = denoise_idx,
    seq_idx = match_to_seq_all(x$sequence, seqs, rc = rc),
    nread = as.integer(x$abundance)
  )
}

#' @rdname make_denoise_map
#' @exportS3Method
make_denoise_map.matrix <- function(x, seqs, rc = FALSE) {
  checkmate::assert_integerish(x)
  checkmate::assert_flag(rc)
  seq_idx <- match_to_seq_all(colnames(x), seqs, rc = rc)
  if (typeof(x) != "integer") {
    mode(x) <- "integer"
  }
  samples <- rownames(x)
  if (is.null(samples)) {
    samples <- as.character(seq_len(nrow(x)))
  }
  out <- vector("list", nrow(x))
  for (i in seq_len(nrow(x))) {
    nz <- which(x[i, ] > 0L)
    if (length(nz) == 0L) {
      out[[i]] <- tibble::tibble(
        sample = character(),
        denoise_idx = integer(),
        seq_idx = integer(),
        nread = integer()
      )
    } else {
      out[[i]] <- tibble::tibble(
        sample = samples[[i]],
        denoise_idx = as.integer(nz),
        seq_idx = seq_idx[nz],
        nread = as.integer(x[i, nz])
      )
    }
  }
  dplyr::bind_rows(out)
}

#' @rdname make_denoise_map
#' @exportS3Method
make_denoise_map.uc_cluster <- function(x, seqs, rc = FALSE) {
  checkmate::assert_class(x, "uc_cluster")
  checkmate::assert_flag(rc)
  if (nrow(x$clusters) == 0L) {
    return(tibble::tibble(
      denoise_idx = integer(),
      seq_idx = integer(),
      nread = integer()
    ))
  }
  tibble::tibble(
    denoise_idx = as.integer(x$clusters$clust_idx),
    seq_idx = match_to_seq_all(x$clusters$seq, seqs, rc = rc),
    nread = as.integer(x$clusters$size)
  )
}

#' @rdname make_denoise_map
#' @exportS3Method
make_denoise_map.list <- function(x, seqs, rc = FALSE) {
  checkmate::assert_flag(rc)

  if (length(x) == 0L) {
    return(tibble::tibble(
      sample = character(),
      denoise_idx = integer(),
      seq_idx = integer(),
      nread = integer()
    ))
  }

  if (checkmate::test_list(x, types = "data.frame")) {
    checkmate::assert_named(x)
    queries <- unlist(
      lapply(x, \(df) {
        if ("accept" %in% names(df)) {
          df <- df[df$accept, , drop = FALSE]
        }
        as.character(df$sequence)
      }),
      use.names = FALSE
    )
    lookup <- seq_idx_lookup(queries, seqs, rc = rc)
    purrr::map_dfr(
      x,
      make_denoise_map.data.frame,
      seqs = lookup,
      rc = rc,
      .id = "sample"
    )
  } else if (checkmate::test_list(x, types = "matrix")) {
    queries <- unlist(lapply(x, colnames), use.names = FALSE)
    lookup <- seq_idx_lookup(queries, seqs, rc = rc)
    purrr::map_dfr(
      x,
      make_denoise_map.matrix,
      seqs = lookup,
      rc = rc
    )
  } else if (all(vapply(x, inherits, logical(1), "uc_cluster"))) {
    checkmate::assert_named(x)
    queries <- unlist(
      lapply(x, \(u) as.character(u$clusters$seq)),
      use.names = FALSE
    )
    lookup <- seq_idx_lookup(queries, seqs, rc = rc)
    purrr::map_dfr(
      x,
      make_denoise_map.uc_cluster,
      seqs = lookup,
      rc = rc,
      .id = "sample"
    )
  } else {
    stop("cannot determine entry type in make_denoise_map.list")
  }
}

#' Project a denoise map to a mapped sequence occurrence table
#'
#' Drops unmatched and zero-abundance rows, then summarizes abundance by
#' `sample` and `seq_idx` when a `sample` column is present.
#'
#' @param denoise_map (`data.frame`) as returned by [make_denoise_map()]
#' @return `data.frame` with columns `sample` (if present in the input),
#'   `seq_idx`, and `nread`
#' @seealso [make_denoise_map()], [make_mapped_sequence_table()]
#' @export
denoise_map_to_seqtable <- function(denoise_map) {
  # avoid R CMD check NOTE for undefined globals due to NSE
  seq_idx <- nread <- sample <- NULL

  checkmate::assert_data_frame(denoise_map)
  checkmate::assert_names(
    names(denoise_map),
    must.include = c("seq_idx", "nread")
  )
  out <- dplyr::filter(denoise_map, !is.na(seq_idx), nread > 0L)
  if ("sample" %in% names(out)) {
    dplyr::summarize(out, nread = sum(nread), .by = c(sample, seq_idx)) |>
      dplyr::select(sample, seq_idx, nread)
  } else {
    dplyr::summarize(out, nread = sum(nread), .by = seq_idx) |>
      dplyr::select(seq_idx, nread)
  }
}
