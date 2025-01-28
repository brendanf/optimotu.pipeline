#' Remove potential tag-jump from sequence table
#' @param seqtable (`data.frame`) sequence table, as returned by
#' `make_long_sequence_table()` or `make_mapped_sequence_table()`
#' @param f (`numeric`) expected cross-talk rate, e.g. 0.01
#' @param p (`numeric`) power to rise the exponent (default, 1; use 1/2 or 1/3
#' to make curves more steep)
#' @param id_col (`character`) name of the column uniquely identifying the
#' sequence
#' @return `data.frame` with columns `sample`, `nread`, `total`, `uncross`, and
#' `is_tag_jump`
#   core by Vladimir Mikryukov,
#   edited for 'targets' by Sten Anslan
#   modified to match OptimOTU style by Brendan Furneaux
remove_tag_jumps <- function(seqtable, f, p, id_col = "seq") {
  # avoid R CMD check NOTE for undeclared globals
  nread <- NULL

  checkmate::assert_data_frame(seqtable)
  checkmate::assert_names(names(seqtable), must.include = c(id_col, "sample", "nread"))
  ## Load ASV table
  cat("...Number of ASVs: ", dplyr::n_distinct(seqtable[[id_col]]), "\n")
  n <- dplyr::n_distinct(seqtable$sample)
  cat("...Number of samples: ", n, "\n")

  ## UNCROSS score (with original parameter - take a root from the exp in denominator, to make curves more steep)
  uncross_score <- function(x, N, n, f = 0.01, tmin = 0.1, p = 1){
    # x = ASV abundance in a sample
    # N = total ASV abundance
    # n = number of samples
    # f = expected cross-talk rate, e.g. 0.01
    # tmin = min score to be considered as cross-talk
    # p = power to rise the exponent (default, 1; use 1/2 or 1/3 to make cureves more stepp)

    z <- f * N / n               # Expected treshold
    sc <- 2 / (1 + exp(x/z)^p)   # t-score
    data.frame(uncross = sc, is_tag_jump = sc >= tmin)
  }

  ## Estimate total abundance of sequence per plate
  out <- seqtable |>
    dplyr::mutate(total = sum(nread, na.rm = TRUE), .by = dplyr::all_of(id_col)) |>
    dplyr::select(-dplyr::all_of(id_col))



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
#' @param seqmap (`data.frame`) sequence map, as returned by `seq_map()`
#' @param seqtable_raw (`data.frame`) raw sequence table, as returned by
#' `make_mapped_sequence_table()`.
#' @param uncross (`data.frame`) uncrossing information, as returned by
#' `remove_tag_jumps()`.
#' @return `data.frame` with the same columns as `seqmap`, but with the `flags`
#' column updated to include the `is_tag_jump` information from `uncross`.
#'
add_uncross_to_seq_map <- function(seqmap, seqtable_raw, uncross) {
  # avoid R CMD check NOTE for undeclared globals
  raw_idx <- seq_idx <- flags <- is_tag_jump <- NULL

  dplyr::left_join(
    seqmap,
    tibble::tibble(
      sample = seqtable_raw$sample,
      seq_idx = seqtable_raw$seq_idx,
      is_tag_jump = uncross$is_tag_jump
    ),
    by = c("sample", "seq_idx")
  ) |>
    dplyr::transmute(
      sample,
      raw_idx,
      seq_idx,
      flags = flags | as.raw(ifelse(is.na(is_tag_jump) | is_tag_jump, 0, 0x08))
    )
}

#' Summarize uncrossing information
#' @param uncross (`data.frame`) uncrossing information, as returned by
#' `remove_tag_jumps()`.
#' @return `data.frame` with columns `sample`, `Total_reads`,
#' `Number_of_TagJump_Events`, `TagJump_reads`, and `ReadPercent_removed`.
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
