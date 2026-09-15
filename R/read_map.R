# SPDX-FileCopyrightText: 2026, Brendan Furneaux
# SPDX-License-Identifier: MIT

#' Precompute a sequence-to-index lookup for a master ASV list
#'
#' Matching against a large `seq_all` (millions of ASVs) is dominated by
#' hashing the table, not by the number of queries. Call this once per chunk
#' with all unique query sequences, then pass the result to
#' [match_to_seq_all()] for per-sample lookups.
#'
#' @param queries (`character`) query sequences; `NA` values are dropped
#' @param seq_all (`character` vector of sequences,
#'   [`XStringSet`][Biostrings::XStringSet-class], or a readable FASTA path)
#'   unique ASV sequences
#' @param rc (`logical`) if `TRUE`, reverse-complement `queries` before
#'   matching. Stored on the result and checked by [match_to_seq_all()].
#' @return object of class `seq_idx_lookup` with elements `keys`, `idx`, and
#'   `rc`
#' @keywords internal
seq_idx_lookup <- function(queries, seq_all, rc = FALSE) {
  checkmate::assert_flag(rc)
  queries <- as.character(queries)
  keys <- unique(queries[!is.na(queries)])
  structure(
    list(
      keys = keys,
      idx = match_to_seq_all(keys, seq_all, rc = rc),
      rc = rc
    ),
    class = "seq_idx_lookup"
  )
}

#' Match sequences to a master ASV list
#'
#' @param sequences (`character`) query sequences, possibly with `NA`
#' @param seq_all (`character` vector of sequences,
#'   [`XStringSet`][Biostrings::XStringSet-class], a readable FASTA path, or a
#'   [seq_idx_lookup()] built for the same `rc`)
#'   unique ASV sequences
#' @param rc (`logical`) if `TRUE`, reverse-complement `sequences` before
#'   matching. Use when queries are reverse-oriented relative to `seq_all`.
#'   Must match the `rc` baked into a `seq_idx_lookup`.
#' @return (`integer`) 1-based indices into `seq_all`; `NA` where there is no
#'   match or the query is `NA`
#' @keywords internal
match_to_seq_all <- function(sequences, seq_all, rc = FALSE) {
  checkmate::assert_flag(rc)
  sequences <- as.character(sequences)
  out <- rep(NA_integer_, length(sequences))
  not_na <- !is.na(sequences)
  if (!any(not_na)) {
    return(out)
  }
  if (inherits(seq_all, "seq_idx_lookup")) {
    if (!identical(rc, seq_all$rc)) {
      stop(
        "match_to_seq_all() rc=",
        rc,
        " does not match seq_idx_lookup rc=",
        seq_all$rc,
        call. = FALSE
      )
    }
    out[not_na] <- seq_all$idx[match(sequences[not_na], seq_all$keys)]
    return(out)
  }
  if (checkmate::test_file_exists(seq_all, "r")) {
    seq_all <- Biostrings::readDNAStringSet(seq_all)
  }
  query <- sequences[not_na]
  if (isTRUE(rc)) {
    query <- as.character(Biostrings::reverseComplement(
      Biostrings::DNAStringSet(query)
    ))
  }
  out[not_na] <- as.integer(BiocGenerics::match(query, seq_all))
  out
}

empty_read_map <- function() {
  tibble::tibble(
    sample = character(),
    raw_idx = integer(),
    seq_idx = integer(),
    flags = raw()
  )
}

# TRUE when the caller passed a 0-sample batch. Path vectors must match
# `sample` in length; files are not checked because there are none.
read_map_empty_batch <- function(sample, ...) {
  checkmate::assert_character(sample, any.missing = FALSE)
  paths <- list(...)
  nms <- names(paths)
  for (nm in nms) {
    checkmate::assert_character(
      paths[[nm]],
      len = length(sample),
      any.missing = FALSE,
      .var.name = nm
    )
  }
  length(sample) == 0L
}

#' Subset a denoise map to one sample
#'
#' @param denoise_map (`data.frame`) as from [make_denoise_map()]
#' @param sample_name (`character` scalar)
#' @param n_samples (`integer` scalar) number of samples in the caller
#' @return rows of `denoise_map` for `sample_name`
#' @noRd
denoise_map_for_sample <- function(denoise_map, sample_name, n_samples) {
  checkmate::assert_data_frame(denoise_map)
  checkmate::assert_names(
    names(denoise_map),
    must.include = c("denoise_idx", "seq_idx")
  )
  if ("sample" %in% names(denoise_map)) {
    return(denoise_map[denoise_map$sample == sample_name, , drop = FALSE])
  }
  checkmate::assert_true(n_samples == 1L)
  denoise_map
}

#' Merge forward and reverse read maps
#'
#' Note that this is only needed for workflows where the sequences are not all
#' in the same orientation, not for ordinary Illumina R1 and R2. (Which are
#' merged earlier using `dada2::mergePairs()` or vsearch mergepairs.)
#'
#' @param read_map_fwd (`data.frame`) forward read map, as returned by
#'   [dada2_read_map()] or [unoise_read_map()]
#' @param read_map_rev (`data.frame`) reverse read map, as returned by
#'   [dada2_read_map()] or [unoise_read_map()]
#' @export
merge_read_maps <- function(read_map_fwd, read_map_rev) {
  # avoid R CMD check NOTE: no visible binding for global variable
  raw_idx <- seq_idx_fwd <- seq_idx_rev <- flags_fwd <- flags_rev <- NULL

  dplyr::full_join(
    read_map_fwd,
    read_map_rev,
    by = c("sample", "raw_idx"),
    suffix = c("_fwd", "_rev")
  ) |>
    dplyr::transmute(
      sample,
      raw_idx,
      seq_idx = dplyr::coalesce(seq_idx_fwd, seq_idx_rev),
      flags = flags_fwd | flags_rev
    )
}

#' Pipe LULU and/or UNCROSS annotations onto a read-map expression
#'
#' Plan-time helper for `targets` commands. Wraps `read_map_expr` in
#' [add_lulu_to_read_map()] and/or [add_uncross_to_read_map()] according to
#' [do_lulu()] and [do_tag_jump()].
#'
#' The result is spliced into a target command with `!!`:
#' ```
#' tar_target(read_map, !!with_read_map_annotate(quote(dada2_read_map(...))))
#' ```
#' Built with [base::substitute()] rather than `rlang::expr(!!x |> f())`.
#' In a pipe, `!!` is parsed as two `!` operators, so that form would run
#' the annotation at plan time instead of capturing an expression.
#'
#' @param read_map_expr (`language`) quoted expression that produces a fate
#'   map as from [dada2_read_map()] or [unoise_read_map()].
#' @param lulu_map (`language`) expression for the LULU map target.
#'   Default `lulu_seq_map`.
#' @param seqtable (`language` or `NULL`) community table passed to
#'   [remove_tag_jumps()]. Default `seqtable_lulu` when LULU is enabled,
#'   otherwise `seqtable_raw`.
#' @param uncross (`language`) expression for the UNCROSS target. Default
#'   `uncross`.
#' @return (`language`) `read_map_expr`, possibly piped through LULU and
#'   UNCROSS annotations.
#' @seealso [add_lulu_to_read_map()], [add_uncross_to_read_map()]
#' @export
with_read_map_annotate <- function(
  read_map_expr,
  lulu_map = quote(lulu_seq_map),
  seqtable = NULL,
  uncross = quote(uncross)
) {
  checkmate::assert(
    checkmate::check_class(read_map_expr, "call"),
    checkmate::check_class(read_map_expr, "name"),
    .var.name = "read_map_expr"
  )
  checkmate::assert(
    checkmate::check_class(lulu_map, "call"),
    checkmate::check_class(lulu_map, "name"),
    .var.name = "lulu_map"
  )
  checkmate::assert(
    checkmate::check_class(uncross, "call"),
    checkmate::check_class(uncross, "name"),
    .var.name = "uncross"
  )
  expr <- read_map_expr
  if (do_lulu()) {
    expr <- substitute(
      EXPR |> optimotu.pipeline::add_lulu_to_read_map(LULU),
      list(EXPR = expr, LULU = lulu_map)
    )
  }
  if (isTRUE(do_tag_jump())) {
    if (is.null(seqtable)) {
      seqtable <- if (do_lulu()) {
        quote(seqtable_lulu)
      } else {
        quote(seqtable_raw)
      }
    } else {
      checkmate::assert(
        checkmate::check_class(seqtable, "call"),
        checkmate::check_class(seqtable, "name"),
        .var.name = "seqtable"
      )
    }
    expr <- substitute(
      EXPR |> optimotu.pipeline::add_uncross_to_read_map(PRE, UNCROSS),
      list(EXPR = expr, PRE = seqtable, UNCROSS = uncross)
    )
  }
  expr
}
