# SPDX-FileCopyrightText: 2026, Brendan Furneaux
# SPDX-License-Identifier: MIT

# UNOISE helpers for per-read mapping. Clustering itself is
# vsearch_cluster_unoise2() in R/vsearch.R; targets call that directly.

# Map each input sequence ID to the dereplicated representative ID.
# Re-runs vsearch --fastx_uniques (intentionally not stored by
# vsearch_cluster_unoise2(), which pipes derep into UNOISE).
vsearch_derep_uc <- function(seq, vsearch = find_vsearch()) {
  checkmate::assert_string(seq)
  checkmate::assert_file_exists(seq, "r")
  if (sequence_size(seq) == 0L) {
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
  if (sequence_size(seq) == 0L) {
    return(character())
  }
  seqs <- if (grepl(fastq_regex, seq)) {
    Biostrings::readBStringSet(seq, format = "fastq")
  } else {
    Biostrings::readBStringSet(seq, format = "fasta")
  }
  sub("\\s.*", "", names(seqs))
}

#' Map the fate of individual reads through UNOISE merge and denoising
#'
#' Analog of [`seq_map()`] for the UNOISE path, where merging and quality
#' filtering happen before denoising.
#'
#' Bit `0x02` (filter) is set when the raw read is present in the merged,
#' quality-filtered FASTQ. Bit `0x04` (denoise & merge) is set when that
#' merged sequence maps to an ASV in `uc`.
#'
#' Dereplication is repeated here on `fq_merged`. That duplicates work already
#' done inside [`vsearch_cluster_unoise2()`], which is an intentional tradeoff
#' between computation (dereplication is fast) and storage (the dereplication
#' map is large).
#'
#' @param sample (`character`) name of the sample
#' @param fq_raw (`character`) name of the raw FASTQ R1 file
#' @param fq_trim (`character`) name of the trimmed FASTQ R1 file
#' @param fq_merged (`character`) name of the merged and filtered FASTQ file
#' @param uc (`uc_cluster`) result of [`vsearch_cluster_unoise2()`] for this
#'   sample
#' @param seq_all (`character` or `XStringSet`) unique ASV sequences
#' @param rc (`logical`) if `TRUE`, centroid sequences in `uc` are
#'   reverse-complemented relative to `seq_all`
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
#' @export
unoise_seq_map <- function(
  sample,
  fq_raw,
  fq_trim,
  fq_merged,
  uc,
  seq_all,
  rc = FALSE
) {
  # avoid R CMD check NOTE: no visible binding for global variable
  raw_idx <- seq_idx <- trim_idx <- filt_idx <- dada_idx <- NULL

  checkmate::assert_string(sample)
  checkmate::assert_file_exists(fq_raw, "r")
  checkmate::assert_file_exists(fq_trim, "r")
  checkmate::assert_file_exists(fq_merged, "r")
  checkmate::assert_class(uc, "uc_cluster")
  checkmate::assert_flag(rc)

  seq_map <- fastq_seq_map(fq_raw, fq_trim, fq_merged)
  merged_ids <- read_seq_ids(fq_merged)
  derep_map <- vsearch_derep_uc(fq_merged)
  merged_seq_id <- merged_ids[seq_map$filt_idx]
  unique_id <- derep_map$unique_id[match(merged_seq_id, derep_map$seq_id)]
  clust_idx <- uc$map$clust_idx[match(unique_id, uc$map$seq_id)]
  centroid_seq <- uc$clusters$seq[match(clust_idx, uc$clusters$clust_idx)]
  if (isTRUE(rc) && length(centroid_seq) > 0L) {
    not_na <- !is.na(centroid_seq)
    centroid_seq[not_na] <- as.character(Biostrings::reverseComplement(
      Biostrings::DNAStringSet(centroid_seq[not_na])
    ))
  }
  if (checkmate::test_file_exists(seq_all, "r")) {
    seq_all <- Biostrings::readDNAStringSet(seq_all)
  }
  seq_all_chr <- as.character(seq_all)
  seq_map$dada_idx <- seq_map$seq_idx <- match(centroid_seq, seq_all_chr)
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
