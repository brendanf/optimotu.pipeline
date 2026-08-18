#' Remove potential tag-jump from sequence table
#' @param seqtable (`data.frame`) sequence table, as returned by
#' `make_long_sequence_table()` or `make_mapped_sequence_table()`
#' @param f (`numeric`) expected cross-talk rate, e.g. 0.01
#' @param p (`numeric`) power to rise the exponent (default, 1; use 1/2 or 1/3
#' to make curves more steep)
#' @param id_col (`character`) name of the column uniquely identifying the
#' sequence
#' @return `data.frame` with the input columns (including `id_col`) plus
#' `total`, `uncross`, and `is_tag_jump`. Row order matches `seqtable`.
#' @export
#   core by Vladimir Mikryukov,
#   edited for 'targets' by Sten Anslan
#   modified to match OptimOTU style by Brendan Furneaux
remove_tag_jumps <- function(seqtable, f, p, id_col = "seq") {
  # avoid R CMD check NOTE for undeclared globals
  nread <- NULL

  checkmate::assert_data_frame(seqtable)
  checkmate::assert_names(
    names(seqtable),
    must.include = c(id_col, "sample", "nread")
  )
  ## Load ASV table
  cat("...Number of ASVs: ", dplyr::n_distinct(seqtable[[id_col]]), "\n")
  n <- dplyr::n_distinct(seqtable$sample)
  cat("...Number of samples: ", n, "\n")

  ## UNCROSS score (with original parameter - take a root from the exp in denominator, to make curves more steep)
  uncross_score <- function(x, N, n, f = 0.01, tmin = 0.1, p = 1) {
    # x = ASV abundance in a sample
    # N = total ASV abundance
    # n = number of samples
    # f = expected cross-talk rate, e.g. 0.01
    # tmin = min score to be considered as cross-talk
    # p = power to rise the exponent (default, 1; use 1/2 or 1/3 to make cureves more stepp)

    z <- f * N / n # Expected treshold
    sc <- 2 / (1 + exp(x / z)^p) # t-score
    data.frame(uncross = sc, is_tag_jump = sc >= tmin)
  }

  ## Estimate total abundance of sequence per plate
  out <- seqtable |>
    dplyr::mutate(
      total = sum(nread, na.rm = TRUE),
      .by = dplyr::all_of(id_col)
    )

  ## Esimate UNCROSS score
  out <- cbind(
    out,
    uncross_score(
      x = out$nread,
      N = out$total,
      n = n,
      f = as.numeric(f),
      p = as.numeric(p)
    )
  )
  cat("...Number of tag-jumps: ", sum(out$is_tag_jump, na.rm = TRUE), "\n")
  # fwrite(x = TJ, file = "TagJump_stats.txt", sep = "\t")

  ## Remove detected tag-jumps from the ASV table
  out
}

#' Add uncrossing information to a sequence map
#'
#' Sets bit `0x08` when the read's current `seq_idx` is present in `uncross`
#' for that sample and `is_tag_jump` is `FALSE`. Join-miss and tag-jumps both
#' leave `0x08` unset. Extra columns on `seqmap` (e.g. `denoise_idx`) are
#' preserved.
#'
#' When `uncross` includes `seq_idx` (the default after
#' [remove_tag_jumps()]), the join uses that column. Otherwise keys are rebuilt
#' from `seqtable_raw` by row position.
#'
#' Call this after [add_lulu_to_seq_map()] so `seq_idx` is the LULU parent
#' when LULU is enabled.
#'
#' @param seqmap (`data.frame`) sequence map, as returned by `seq_map()`
#' @param seqtable_raw (`data.frame`) sequence table passed to
#'   [remove_tag_jumps()]. Used to rebuild join keys when `uncross` does not
#'   contain `seq_idx`.
#' @param uncross (`data.frame`) uncrossing information, as returned by
#' `remove_tag_jumps()`.
#' @return `data.frame` with the same columns as `seqmap`, but with the `flags`
#' column updated to include the `is_tag_jump` information from `uncross`.
#' @seealso [with_seqmap_annotate()], [add_lulu_to_seq_map()]
#' @export
add_uncross_to_seq_map <- function(seqmap, seqtable_raw, uncross) {
  # avoid R CMD check NOTE for undeclared globals
  flags <- is_tag_jump <- seq_idx <- NULL

  checkmate::assert_data_frame(seqmap)
  checkmate::assert_names(
    names(seqmap),
    must.include = c("sample", "seq_idx", "flags")
  )
  checkmate::assert_data_frame(uncross)
  checkmate::assert_names(
    names(uncross),
    must.include = c("sample", "is_tag_jump")
  )

  uncross_keys <- if ("seq_idx" %in% names(uncross)) {
    dplyr::select(uncross, dplyr::all_of(c("sample", "seq_idx", "is_tag_jump")))
  } else {
    checkmate::assert_data_frame(seqtable_raw)
    checkmate::assert_names(
      names(seqtable_raw),
      must.include = c("sample", "seq_idx")
    )
    tibble::tibble(
      sample = seqtable_raw$sample,
      seq_idx = seqtable_raw$seq_idx,
      is_tag_jump = uncross$is_tag_jump
    )
  }

  dplyr::left_join(
    seqmap,
    uncross_keys,
    by = c("sample", "seq_idx")
  ) |>
    dplyr::mutate(
      flags = flags | as.raw(ifelse(is.na(is_tag_jump) | is_tag_jump, 0, 0x08))
    ) |>
    dplyr::select(-is_tag_jump)
}

#' Summarize uncrossing information
#' @param uncross (`data.frame`) uncrossing information, as returned by
#' `remove_tag_jumps()`.
#' @return `data.frame` with columns `sample`, `Total_reads`,
#' `Number_of_TagJump_Events`, `TagJump_reads`, and `ReadPercent_removed`.
#' @export
summarize_uncross <- function(uncross) {
  # avoid R CMD check NOTE for undeclared globals
  nread <- is_tag_jump <- TagJump_reads <- Total_reads <- NULL
  uncross |>
    dplyr::summarize(
      Total_reads = sum(nread),
      Number_of_TagJump_Events = sum(is_tag_jump),
      TagJump_reads = sum(nread[is_tag_jump], na.rm = TRUE),
      ReadPercent_removed <- TagJump_reads / Total_reads * 100,
      .by = sample
    )
}
