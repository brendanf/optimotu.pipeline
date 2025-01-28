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
  if (all(
    methods::is(dadaF, "dada"),
    methods::is(dadaR, "dada"),
    methods::is(derepF, "derep"),
    methods::is(derepR, "derep"),
    methods::is(merged, "data.frame")
  )) {
    tibble::tibble(
      fwd_idx = dadaF$map[derepF$map],
      rev_idx = dadaR$map[derepR$map]
    ) |>
      dplyr::left_join(
        tibble::rowid_to_column(merged[c("forward", "reverse")], "merge_idx"),
        by = c("fwd_idx" = "forward", "rev_idx" = "reverse")
      )
  } else if (all(
    rlang::is_bare_list(dadaF),
    rlang::is_bare_list(dadaR),
    rlang::is_bare_list(derepF),
    rlang::is_bare_list(derepR),
    rlang::is_bare_list(merged)
  )) {
    purrr::pmap(list(dadaF, derepF, dadaR, derepR, merged), dada_merge_map)
  }
}

#' Map the fate of individual reads through merging to find unique reads
#'
#' @param sample (`character`) name of the sample
#' @param fq_raw (`character`) name of the raw fastq R1 file
#' @param fq_trim (`character`) name of the trimmed fastq R1 file
#' @param fq_filt (`character`) name of the filtered fastq R1 file
#' @param dadaF ([`dada2::dada-class`]) denoised R1
#' @param derepF ([`dada2::derep-class`]) dereplicated R1
#' @param dadaR ([`dada2::dada-class`]) denoised R2
#' @param derepR ([`dada2::dada-class`]) dereplicated R2
#' @param merged (`data.frame` as returned by `dada2::mergePairs()`) result of
#' merging `dadaF` and `dadaR`
#' @param seq_all (`character`) unique ASV sequences
#' @param rc (`logical`) if `TRUE`, sequences in `merged` are reverse-complemented
#'  relative to `seq_all`.
#'
#' @return `data.frame` with columns:
#'  -`sample` (character) the sample name
#'  -`raw_idx` (integer) the index of the sequence in the raw file; see `seq_map()`
#'  -`seq_idx` (integer) the index of the sequence in seq_all
#'  -`flags` (raw) bitset indicating the presence of the sequence at different stages:
#'    0x01 = trimmed
#'    0x02 = filtered
#'    0x04 = denoised & merged
seq_map <- function(sample, fq_raw, fq_trim, fq_filt, dadaF, derepF, dadaR,
                    derepR, merged, seq_all, rc = FALSE) {
  # avoid R CMD check NOTE: no visible binding for global variable
  raw_idx <- seq_idx <- trim_idx <- filt_idx <- dada_idx <- NULL

  seq_map <- fastq_seq_map(fq_raw, fq_trim, fq_filt)
  dada_map <- dada_merge_map(dadaF, derepF, dadaR, derepR, merged)
  seq_map$dada_idx <-
    seq_map$seq_idx <- match(merged$sequence, seq_all)[dada_map$merge_idx[seq_map$filt_idx]]
  dplyr::transmute(
    seq_map,
    sample = sample,
    raw_idx,
    seq_idx,
    flags = as.raw(
      ifelse(is.na(trim_idx), 0, 0x01) +
        ifelse(is.na(filt_idx), 0, 0x02) +
        ifelse(is.na(dada_idx), 0, 0x04)
    )
  )
}

#' Merge forward and reverse sequence maps
#'
#' Note that this is only needed for workflows where the sequences are not all in
#' the same orientation, not for ordinary Illumina R1 and R2. (Which are merged
#' earlier using `dada2::mergePairs()`.)
#'
#' @param seqmap_fwd (`data.frame`) forward sequence map, as returned by `seq_map()`
#' @param seqmap_rev (`data.frame`) reverse sequence map, as returned by `seq_map()`
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
