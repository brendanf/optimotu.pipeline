#' Map the fate of individual reads through dada2 dereplication, denoising, and merge.
#'
#' @param dadaF ([`dada2::dada-class`] object or list of such objects) denoised
#' forward reads
#' @param derepF ([`dada2::derep-class`] object or list of such objects)
#' dereplicated forward reads
#' @param dadaR ([`dada2::dada-class`] object or list of such objects) denoised
#' reverse reads
#' @param derepR ([`dada2::derep-class`] object or list of such objects)
#' dereplicated reverse reads
#' @param merged (`data.frame` returned by `dada2::mergePairs()` or list of
#' such objects) results of merginf the denoised reads in dadaF and dadaR
#'
#' @return a `data.frame` with three columns:
#'   - `fwd_idx` (integer) index of forward ASV in `dadaF`
#'   - `rev_idx` (integer) index of reverse ASV in `dadaR`
#'   - `merge_idx` (integer) row index of merged ASV `merged`
#' Each row of this `data.frame` represents a single read in the fastq files
#' originally passed to `dada2::derepFastq()`, and the rows are in the same
#' order as the reads.
#' If the inputs were lists, then the output is a list of `data.frame`s as
#' described above.
dada_merge_map <- function(dadaF, derepF, dadaR, derepR, merged) {
  if (
    all(
      methods::is(dadaF, "dada"),
      methods::is(dadaR, "dada"),
      methods::is(derepF, "derep"),
      methods::is(derepR, "derep"),
      methods::is(merged, "data.frame")
    )
  ) {
    tibble::tibble(
      fwd_idx = dadaF$map[derepF$map],
      rev_idx = dadaR$map[derepR$map]
    ) |>
      dplyr::left_join(
        tibble::rowid_to_column(merged[c("forward", "reverse")], "merge_idx"),
        by = c("fwd_idx" = "forward", "rev_idx" = "reverse")
      )
  } else if (
    all(
      rlang::is_bare_list(dadaF),
      rlang::is_bare_list(dadaR),
      rlang::is_bare_list(derepF),
      rlang::is_bare_list(derepR),
      rlang::is_bare_list(merged)
    )
  ) {
    purrr::pmap(list(dadaF, derepF, dadaR, derepR, merged), dada_merge_map)
  }
}

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

empty_seq_map <- function() {
  tibble::tibble(
    sample = character(),
    raw_idx = integer(),
    seq_idx = integer(),
    flags = raw()
  )
}

#' Map the fate of individual reads through merging to find unique reads
#'
#' Accepts a single sample or a chunk of samples. For a chunk, all unique
#' sequences in `merged` are matched to `seq_all` once, then each sample is
#' processed against that shared lookup.
#'
#' @param sample (`character`) sample name(s)
#' @param fq_raw (`character`) raw FASTQ R1 file path(s)
#' @param fq_trim (`character`) trimmed FASTQ R1 file path(s)
#' @param fq_filt (`character`) filtered FASTQ R1 file path(s)
#' @param dadaF ([`dada2::dada-class`] or list of such) denoised R1
#' @param derepF ([`dada2::derep-class`] or list of such) dereplicated R1
#' @param dadaR ([`dada2::dada-class`] or list of such) denoised R2
#' @param derepR ([`dada2::derep-class`] or list of such) dereplicated R2
#' @param merged (`data.frame` as returned by `dada2::mergePairs()`, or a named
#'   list of such) result of merging `dadaF` and `dadaR`
#' @param seq_all (`character` vector of sequences,
#'   [`XStringSet`][Biostrings::XStringSet-class], or a readable FASTA path,
#'   e.g. a `tar_file` target) unique ASV sequences
#' @param rc (`logical`) if `TRUE`, sequences in `merged` are reverse-complemented
#'  relative to `seq_all`.
#'
#' @return `data.frame` with columns:
#'  - `sample` (character) the sample name
#'  - `raw_idx` (integer) the index of the sequence in the raw file
#'  - `seq_idx` (integer) the index of the sequence in `seq_all`
#'  - `flags` (raw) bitset indicating the presence of the sequence at different
#'    stages:
#'    0x01 = trimmed
#'    0x02 = filtered
#'    0x04 = denoised & merged
#'    0x08 = survived UNCROSS (set later by [add_uncross_to_seq_map()] when
#'      `is_tag_jump` is `FALSE`; not set here)
#'
#' Bits `0x10`--`0x80` are reserved for ASV-level filter results
#' (`asv_map$result`: chimera/spike/model), not per-read fate flags.
#' @export
seq_map <- function(
  sample,
  fq_raw,
  fq_trim,
  fq_filt,
  dadaF,
  derepF,
  dadaR,
  derepR,
  merged,
  seq_all,
  rc = FALSE
) {
  # avoid R CMD check NOTE: no visible binding for global variable
  raw_idx <- seq_idx <- trim_idx <- filt_idx <- dada_idx <- NULL

  checkmate::assert_character(sample, min.len = 1L, any.missing = FALSE)
  checkmate::assert_character(fq_raw, len = length(sample), any.missing = FALSE)
  checkmate::assert_character(
    fq_trim,
    len = length(sample),
    any.missing = FALSE
  )
  checkmate::assert_character(
    fq_filt,
    len = length(sample),
    any.missing = FALSE
  )
  checkmate::assert_file_exists(fq_raw, "r")
  checkmate::assert_file_exists(fq_trim, "r")
  checkmate::assert_file_exists(fq_filt, "r")
  checkmate::assert_flag(rc)

  if (methods::is(merged, "data.frame")) {
    checkmate::assert_true(length(sample) == 1L)
    merged <- stats::setNames(list(merged), sample)
    dadaF <- list(dadaF)
    derepF <- list(derepF)
    dadaR <- list(dadaR)
    derepR <- list(derepR)
  } else {
    checkmate::assert_list(merged, len = length(sample))
    checkmate::assert_list(dadaF, len = length(sample))
    checkmate::assert_list(derepF, len = length(sample))
    checkmate::assert_list(dadaR, len = length(sample))
    checkmate::assert_list(derepR, len = length(sample))
  }

  queries <- unlist(
    lapply(merged, \(m) as.character(m$sequence)),
    use.names = FALSE
  )
  lookup <- seq_idx_lookup(queries, seq_all, rc = rc)

  out <- vector("list", length(sample))
  for (i in seq_along(sample)) {
    smap <- fastq_seq_map(fq_raw[[i]], fq_trim[[i]], fq_filt[[i]])
    dmap <- dada_merge_map(
      dadaF[[i]],
      derepF[[i]],
      dadaR[[i]],
      derepR[[i]],
      merged[[i]]
    )
    merge_seq_idx <- match_to_seq_all(
      merged[[i]]$sequence,
      lookup,
      rc = rc
    )
    smap$dada_idx <- smap$seq_idx <-
      merge_seq_idx[dmap$merge_idx[smap$filt_idx]]
    out[[i]] <- dplyr::transmute(
      smap,
      sample = sample[[i]],
      raw_idx,
      seq_idx,
      flags = as.raw(
        ifelse(is.na(trim_idx), 0, 0x01) +
          ifelse(is.na(filt_idx), 0, 0x02) +
          ifelse(is.na(dada_idx), 0, 0x04)
      )
    )
  }
  if (length(out) == 0L) {
    return(empty_seq_map())
  }
  dplyr::bind_rows(out)
}

#' Merge forward and reverse sequence maps
#'
#' Note that this is only needed for workflows where the sequences are not all in
#' the same orientation, not for ordinary Illumina R1 and R2. (Which are merged
#' earlier using `dada2::mergePairs()`.)
#'
#' @param seqmap_fwd (`data.frame`) forward sequence map, as returned by `seq_map()`
#' @param seqmap_rev (`data.frame`) reverse sequence map, as returned by `seq_map()`
#' @export
merge_seq_maps <- function(seqmap_fwd, seqmap_rev) {
  # avoid R CMD check NOTE: no visible binding for global variable
  raw_idx <- seq_idx_fwd <- seq_idx_rev <- flags_fwd <- flags_rev <- NULL

  dplyr::full_join(
    seqmap_fwd,
    seqmap_rev,
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

#' Pipe LULU and/or UNCROSS annotations onto a fate-map expression
#'
#' Plan-time helper for `targets` commands. Wraps `seqmap_expr` in
#' [add_lulu_to_seq_map()] and/or [add_uncross_to_seq_map()] according to
#' [do_lulu()] and [do_tag_jump()].
#'
#' The result is spliced into a target command with `!!`:
#' ```
#' tar_target(dada_map, !!with_seqmap_annotate(quote(seq_map(...))))
#' ```
#' Built with [base::substitute()] rather than `rlang::expr(!!x |> f())`.
#' In a pipe, `!!` is parsed as two `!` operators, so that form would run
#' the annotation at plan time instead of capturing an expression.
#'
#' @param seqmap_expr (`language`) quoted expression that produces a fate
#'   map as from [seq_map()] or [unoise_seq_map()].
#' @param lulu_map (`language`) expression for the LULU map target.
#'   Default `lulu_asv_map`.
#' @param seqtable (`language` or `NULL`) community table passed to
#'   [remove_tag_jumps()]. Default `seqtable_lulu` when LULU is enabled,
#'   otherwise `seqtable_raw`.
#' @param uncross (`language`) expression for the UNCROSS target. Default
#'   `uncross`.
#' @return (`language`) `seqmap_expr`, possibly piped through LULU and
#'   UNCROSS annotations.
#' @seealso [add_lulu_to_seq_map()], [add_uncross_to_seq_map()]
#' @export
with_seqmap_annotate <- function(
  seqmap_expr,
  lulu_map = quote(lulu_asv_map),
  seqtable = NULL,
  uncross = quote(uncross)
) {
  checkmate::assert(
    checkmate::check_class(seqmap_expr, "call"),
    checkmate::check_class(seqmap_expr, "name"),
    .var.name = "seqmap_expr"
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
  expr <- seqmap_expr
  if (do_lulu()) {
    expr <- substitute(
      EXPR |> optimotu.pipeline::add_lulu_to_seq_map(LULU),
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
      EXPR |> optimotu.pipeline::add_uncross_to_seq_map(PRE, UNCROSS),
      list(EXPR = expr, PRE = seqtable, UNCROSS = uncross)
    )
  }
  expr
}
