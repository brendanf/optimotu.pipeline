# SPDX-FileCopyrightText: 2026, Brendan Furneaux
# SPDX-License-Identifier: MIT

# UNOISE helpers for per-read mapping. Clustering itself is
# vsearch_cluster_unoise2() in R/vsearch.R; targets call that directly.

# Map each input sequence ID to the dereplicated representative ID.
# Re-runs vsearch --fastx_uniques (intentionally not stored by
# vsearch_cluster_unoise2(), which pipes derep into UNOISE).
# When `merged_ids` is supplied, emptiness is taken from that vector so
# sequence_size() is not repeated after a prior FASTQ name read.
vsearch_derep_uc <- function(
  seq,
  vsearch = find_vsearch(),
  merged_ids = NULL
) {
  checkmate::assert_string(seq)
  checkmate::assert_file_exists(seq, "r")
  empty <- if (is.null(merged_ids)) {
    sequence_size(seq) == 0L
  } else {
    length(merged_ids) == 0L
  }
  if (isTRUE(empty)) {
    return(tibble::tibble(seq_id = character(), unique_id = character()))
  }
  derep_out <- withr::local_tempfile(pattern = "derep", fileext = ".fasta")
  uc <- processx::run(
    vsearch,
    c(
      "--fastx_uniques",
      seq,
      "--uc",
      "-",
      "--fastaout",
      derep_out
    ),
    error_on_status = TRUE
  )$stdout
  if (!nzchar(uc)) {
    return(tibble::tibble(seq_id = character(), unique_id = character()))
  }
  parsed <- readr::read_tsv(
    I(uc),
    col_names = c(
      "type",
      "clust_idx",
      "size",
      "sim",
      "strand",
      "na1",
      "na2",
      "na3",
      "seq_id",
      "hit_id"
    ),
    col_types = "c-------cc",
    na = c("NA", "*")
  )
  parsed <- parsed[parsed$type %in% c("S", "H"), ]
  parsed$seq_id <- sub(";size=[0-9]+", "", parsed$seq_id)
  parsed$hit_id <- sub(";size=[0-9]+", "", parsed$hit_id)
  tibble::tibble(
    seq_id = parsed$seq_id,
    unique_id = ifelse(parsed$type == "S", parsed$seq_id, parsed$hit_id)
  )
}

read_seq_ids <- function(seq) {
  # Prefer the C++ name reader for FASTQ; fall back to Biostrings for FASTA.
  if (grepl(fastq_regex, seq)) {
    return(as.character(fastq_names(seq)))
  }
  if (sequence_size(seq) == 0L) {
    return(character())
  }
  seqs <- Biostrings::readBStringSet(seq, format = "fasta")
  sub("\\s.*", "", names(seqs))
}

#' Map the fate of individual reads through UNOISE merge and denoising
#'
#' Analog of [dada2_read_map()] for the UNOISE path, where merging and quality
#' filtering happen before denoising.
#'
#' Accepts a single sample or a chunk of samples. Sequence-to-index matching
#' against `seq_all` must already have been done by [make_denoise_map()];
#' pass that result as `denoise_map`.
#'
#' Bit `0x02` (filter) is set when the raw read is present in the merged,
#' quality-filtered FASTQ. Bit `0x04` (denoise) is set when that merged
#' sequence maps to an ASV in `uc`.
#'
#' Dereplication is repeated here on `fq_merged`. That duplicates work already
#' done inside [vsearch_cluster_unoise2()], which is an intentional tradeoff
#' between computation (dereplication is fast) and storage (the dereplication
#' map is large).
#'
#' @param sample (`character`) sample name(s)
#' @param fq_raw (`character`) raw FASTQ R1 file path(s)
#' @param fq_trim (`character`) trimmed FASTQ R1 file path(s)
#' @param fq_merged (`character`) merged and filtered FASTQ file path(s)
#' @param uc (`uc_cluster` or named list of such) result of
#'   [vsearch_cluster_unoise2()] for the sample(s)
#' @param denoise_map (`data.frame`) as returned by [make_denoise_map()] for
#'   the same `uc` object(s); must include `denoise_idx` and `seq_idx`, and
#'   `sample` when mapping more than one sample
#'
#' @return `data.frame` with columns:
#'  - `sample` (character) the sample name
#'  - `raw_idx` (integer) the index of the sequence in the raw file
#'  - `seq_idx` (integer) the index of the sequence in `seq_all`
#'  - `flags` (raw) bitset indicating the presence of the sequence at different
#'    stages:
#'    0x01 = trimmed
#'    0x02 = merged and quality-filtered
#'    0x04 = denoised
#'    0x08 = survived UNCROSS (set later by [add_uncross_to_read_map()] when
#'      `is_tag_jump` is `FALSE`; not set here)
#'
#' Bits `0x10`--`0x80` are reserved for ASV-level filter results
#' (`asv_map$result`: chimera/spike/model), not per-read fate flags.
#' @export
unoise_read_map <- function(
  sample,
  fq_raw,
  fq_trim,
  fq_merged,
  uc,
  denoise_map
) {
  # avoid R CMD check NOTE: no visible binding for global variable
  raw_idx <- seq_idx <- trim_idx <- filt_idx <- denoise_local <- NULL

  checkmate::assert_character(sample, min.len = 1L, any.missing = FALSE)
  checkmate::assert_character(fq_raw, len = length(sample), any.missing = FALSE)
  checkmate::assert_character(
    fq_trim,
    len = length(sample),
    any.missing = FALSE
  )
  checkmate::assert_character(
    fq_merged,
    len = length(sample),
    any.missing = FALSE
  )
  checkmate::assert_file_exists(fq_raw, "r")
  checkmate::assert_file_exists(fq_trim, "r")
  checkmate::assert_file_exists(fq_merged, "r")
  checkmate::assert_data_frame(denoise_map)

  if (inherits(uc, "uc_cluster")) {
    checkmate::assert_true(length(sample) == 1L)
    uc <- stats::setNames(list(uc), sample)
  } else {
    checkmate::assert_list(uc, len = length(sample))
    for (u in uc) {
      checkmate::assert_class(u, "uc_cluster")
    }
  }

  out <- vector("list", length(sample))
  for (i in seq_along(sample)) {
    smap <- fastq_seq_map(fq_raw[[i]], fq_trim[[i]], fq_merged[[i]])
    merged_ids <- read_seq_ids(fq_merged[[i]])
    derep_map <- vsearch_derep_uc(fq_merged[[i]], merged_ids = merged_ids)
    u <- uc[[i]]
    dm_i <- denoise_map_for_sample(denoise_map, sample[[i]], length(sample))
    merged_seq_id <- merged_ids[smap$filt_idx]
    unique_id <- derep_map$unique_id[match(merged_seq_id, derep_map$seq_id)]
    smap$denoise_local <- u$map$clust_idx[match(unique_id, u$map$seq_id)]
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
