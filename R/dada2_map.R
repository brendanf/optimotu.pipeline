# SPDX-FileCopyrightText: 2026, Brendan Furneaux
# SPDX-License-Identifier: MIT

#' Map the fate of individual reads through dada2 dereplication, denoising,
#' and merge.
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
#' such objects) results of merging the denoised reads in dadaF and dadaR
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
#' @keywords internal
dada2_merge_map <- function(dadaF, derepF, dadaR, derepR, merged) {
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
    purrr::pmap(list(dadaF, derepF, dadaR, derepR, merged), dada2_merge_map)
  }
}

#' Map the fate of individual reads through DADA2 denoising and merge
#'
#' Accepts a single sample or a chunk of samples. Sequence-to-index matching
#' against `seq_all` must already have been done by [make_denoise_map()];
#' pass that result as `denoise_map`. A length-0 `sample` (empty seq batch)
#' returns a 0-row fate map with the usual columns.
#'
#' @param sample (`character`) sample name(s); length 0 is allowed and
#'   returns an empty fate map
#' @param fq_raw (`character`) raw FASTQ R1 file path(s)
#' @param fq_trim (`character`) trimmed FASTQ R1 file path(s)
#' @param fq_filt (`character`) filtered FASTQ R1 file path(s)
#' @param dadaF ([`dada2::dada-class`] or list of such) denoised R1
#' @param derepF ([`dada2::derep-class`] or list of such) dereplicated R1
#' @param dadaR ([`dada2::dada-class`] or list of such) denoised R2
#' @param derepR ([`dada2::derep-class`] or list of such) dereplicated R2
#' @param merged (`data.frame` as returned by `dada2::mergePairs()`, or a named
#'   list of such) result of merging `dadaF` and `dadaR`
#' @param denoise_map (`data.frame`) as returned by [make_denoise_map()] for
#'   the same `merged` object(s); must include `denoise_idx` and `seq_idx`,
#'   and `sample` when mapping more than one sample
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
#'    0x08 = survived UNCROSS (set later by [add_uncross_to_read_map()] when
#'      `is_tag_jump` is `FALSE`; not set here)
#'
#' Bits `0x10`--`0x80` are reserved for ASV-level filter results
#' (`asv_map$result`: chimera/spike/model), not per-read fate flags.
#' @export
dada2_read_map <- function(
  sample,
  fq_raw,
  fq_trim,
  fq_filt,
  dadaF,
  derepF,
  dadaR,
  derepR,
  merged,
  denoise_map
) {
  # avoid R CMD check NOTE: no visible binding for global variable
  raw_idx <- seq_idx <- trim_idx <- filt_idx <- denoise_local <- NULL

  if (
    read_map_empty_batch(
      sample,
      fq_raw = fq_raw,
      fq_trim = fq_trim,
      fq_filt = fq_filt
    )
  ) {
    return(empty_read_map())
  }
  checkmate::assert_file_exists(fq_raw, "r")
  checkmate::assert_file_exists(fq_trim, "r")
  checkmate::assert_file_exists(fq_filt, "r")
  checkmate::assert_data_frame(denoise_map)

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

  out <- vector("list", length(sample))
  for (i in seq_along(sample)) {
    smap <- fastq_seq_map(fq_raw[[i]], fq_trim[[i]], fq_filt[[i]])
    dmap <- dada2_merge_map(
      dadaF[[i]],
      derepF[[i]],
      dadaR[[i]],
      derepR[[i]],
      merged[[i]]
    )
    dm_i <- denoise_map_for_sample(denoise_map, sample[[i]], length(sample))
    smap$denoise_local <- dmap$merge_idx[smap$filt_idx]
    smap$seq_idx <- dm_i$seq_idx[match(smap$denoise_local, dm_i$denoise_idx)]
    out[[i]] <- dplyr::transmute(
      smap,
      sample = sample[[i]],
      raw_idx,
      seq_idx,
      flags = as.raw(
        ifelse(is.na(trim_idx), 0, 0x01) +
          ifelse(is.na(filt_idx), 0, 0x02) +
          ifelse(is.na(denoise_local), 0, 0x04)
      )
    )
  }
  if (length(out) == 0L) {
    return(empty_read_map())
  }
  dplyr::bind_rows(out)
}
