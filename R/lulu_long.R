# Implementation of LULU algorithm for a long OTU table.
# Written without reference to the LULU code, based only on algorithm
# description in the manuscript and documentation.

#' LULU secondary denoising
#'
#' The LULU algorithm is described by Frøslev et al. (2017) and implemented in
#' lulu::lulu()`.  This is an independent implementation of the algorithm as
#' described in that paper, adapted to work with long-format OTU tables. It also
#' assumes that the match table was generated on a per-sample basis.  This means
#' that there should be multiple entries for each pair of co-occurring OTUs if
#' they co-occur in multiple samples.  Although it may seem that this means
#' that the match list and the computation required to generate it are larger
#' than necessary, it is in fact smaller and more efficient in most cases,
#' because the majority of OTU pairs do not co-occur in any samples, and so
#' there is no need to compute or store their distances. Additionally computing
#' the distances sample-wise allows easy parallelization of that step.
#'
#' There is a [longstanding issue](https://github.com/tobiasgf/lulu/issues/8)
#' in the original LULU implementation which prevents it from actually
#' identifying "daughter" sequences when they occur less than 100% of the time
#' with their "parent" sequence, even though the default value for
#' `min_cooccurrence` is 0.95. This implementation, like
#' [mumu](https://github.com/frederic-mahe/mumu), implements the algorithm as
#' described rather than as implemented in `lulu::lulu()`, but other than this
#' it is tested to give identical results to `lulu::lulu()` when run on the same
#' data.
#'
#' This function implements only the core LULU algorithm, meaning it determines
#' which OTUs are to be merged.  The actual merging of the OTU table is then
#' performed by `lulu_table()`.  Like the original implementation, calculating
#' pairwise distances is left to the user, e.g. using BLAST, USEARCH, VSEARCH,
#' etc. [optimotu::seq_distmx] is also a convenient way.
#'
#' @references Frøslev, T. G., Kjøller, R., Bruun, H. H., Ejrnæs, R., Brunbjerg,
#' A. K., Pietroni, C., & Hansen, A. J. (2017). Algorithm for post-clustering
#' curation of DNA amplicon data yields reliable biodiversity estimates.
#' Nature Communications, 8(1), 1188.
#' <https://www.nature.com/articles/s41467-017-01312-x>
#'
#' @param otu_table (`data.frame`) the long-format table to denoise. Must have
#' columns `seq_id` (`character`) xor `seq_idx` (`integer`), as well as
#' `nread` (`integer`). Additional columns such as the sample id or metadata are
#' allowed but ignored.
#' @param match_table (`data.frame`) table of pairwise matches between
#' co-occuring OTUs. Must have columns `seq_id1`and `seq_id2` (both `character`)
#' xor `seq_idx1` and `seq_idx2` (both `integer`), as well as `nread1` and
#' `nread2` (`integer`) and `dist` (`numeric`). Additional columns are ignored.
#' @param max_dist (`numeric` scalar) maximum pairwise distance for two OTUs to
#' be considered for merging. This may be given as a percentage between 0 and
#' 100, as a fraction between 0.0 and 1.0, or any other non-negative scale, but
#' should be on the same scale as the values of `dist` in `match_table`.
#' Default: `0.1` (intended as a fraction).
#' @param min_abundance_ratio (`numeric` scalar) minimum ratio of
#' "parent":"daughter" abundances for two OTUs to be considered for merging.
#' The exact meaning is dependent on the value of `use_mean_abundance_ratio`.
#' Default: `1`.
#' @param min_cooccurrence_ratio (`numeric` scalar) minimum ratio of
#' co-occurrences to total occurrences of the "daughter" for two OTUs to be
#' considered for merging. Default: `1`.
#' @param use_mean_abundance_ratio (`logical` flag) if `TRUE`,
#' `min_abundance_ratio` refers to the mean of the relative abundances between
#' the "parent" and "daughter" sequences in all samples where they co-occur. If
#' `FALSE`, the minimum is enforced for every sample. Default: `FALSE`.
#' @param id_is_int (`logical` flag) if `TRUE`, the OTU identifiers are
#' integers; in this case they are expected to be in columns `seq_idx*`. If
#' `FALSE`, the OTU identifiers are character strings, and are expected to be
#' in columns `seq_id*`. Default: `TRUE` if `otu_table` has a column `seq_idx`,
#' `FALSE` otherwise.
#' @param id_is_sorted (`logical` flag) if `TRUE`, the values of `seq_id*` or
#' `seq_idx*` are sorted in decreasing order of prevalence (number of samples),
#' with ties broken by decreasing total abundance (`nread`). This skips the
#' computational cost of sorting the OTUs, but will cause an error if they are
#' not sorted as described. Default: `TRUE` if `id_is_int` is `TRUE`, otherwise
#' `FALSE`.
#'
#' @returns a two-column `data.frame` with columns `seq_idx` and `lulu_idx`
#' or `seq_id` and `lulu_id`. `seq_id*` includes all values which occur in
#' `otu_table`, and `lulu_id*` gives the corresponding value in `lulu_table`;
#' these differ only for sequences which were determined to be minor variants
#' of another sequence.
#' @export
lulu_map <- function(
  otu_table,
  match_table,
  max_dist = lulu_max_dist(),
  min_abundance_ratio = lulu_min_abundance_ratio(),
  min_cooccurrence_ratio = lulu_min_cooccurrence_ratio(),
  use_mean_abundance_ratio = lulu_use_mean_abundance_ratio(),
  id_is_int = "seq_idx" %in% names(otu_table),
  id_is_sorted = id_is_int
) {
  checkmate::assert_flag(id_is_int)
  id_col <- if (id_is_int) "seq_idx" else "seq_id"
  id_col_sym <- rlang::sym(id_col)
  id1_col <- paste0(id_col, "1")
  id1_col_sym <- rlang::sym(id1_col)
  id2_col <- paste0(id_col, "2")
  id2_col_sym <- rlang::sym(id2_col)

  checkmate::assert_data_frame(otu_table)
  checkmate::assert_names(
    names(otu_table),
    must.include = c("nread", id_col)
  )
  if (id_is_int) {
    checkmate::assert_integerish(otu_table$seq_idx, lower = 1L)
  } else {
    checkmate::assert_character(otu_table$seq_id)
  }
  if ("seq_id" %in% names(otu_table) && "seq_idx" %in% names(otu_table)) {
    stop("`otu_table` may have column `seq_id` or `seq_idx`, but not both")
  }
  checkmate::assert_integerish(otu_table$nread, lower = 1L)

  checkmate::assert_data_frame(match_table)
  checkmate::assert_names(
    names(match_table),
    must.include = c(
      c("nread1", "nread2", "dist", id1_col, id2_col)
    )
  )
  checkmate::assert_integerish(match_table$nread1, lower = 1)
  checkmate::assert_integerish(match_table$nread2, lower = 1)
  checkmate::assert_numeric(match_table$dist, lower = 0)
  if (id_is_int) {
    checkmate::assert_names(
      names(match_table),
      disjunct.from = c("seq_id1", "seq_id2")
    )
    checkmate::assert_integerish(match_table$seq_idx1, lower = 1L)
    checkmate::assert_integerish(match_table$seq_idx2, lower = 1L)
  } else {
    checkmate::assert_names(
      names(match_table),
      disjunct.from = c("seq_idx1", "seq_idx2")
    )
    checkmate::assert_character(match_table$seq_id1)
    checkmate::assert_character(match_table$seq_id2)
  }
  checkmate::assert_numeric(max_dist, lower = 0)
  checkmate::assert_numeric(min_abundance_ratio, lower = 0, upper = 1)

  match_table <- dplyr::filter(
    match_table,
    dist <= max_dist
  )

  # Sort the otu table if it is not already sorted
  if (isFALSE(id_is_sorted)) {
    seq_counts <- dplyr::count(otu_table, !!id_col_sym, name = "nsample")
    otu_table <- otu_table |>
      dplyr::left_join(seq_counts, by = id_col) |>
      dplyr::mutate(nread_total = sum(nread), .by = !!id_col_sym) |>
      dplyr::arrange(
        dplyr::desc(nsample),
        dplyr::desc(nread_total),
        !!id_col_sym
      )
    # Make an ordered factor with the levels sorted.
    otu_table <- dplyr::mutate(
      otu_table,
      !!id_col_sym := ordered(!!id_col_sym, levels = unique(!!id_col_sym))
    )
    match_table <- dplyr::mutate(
      match_table,
      !!id1_col_sym := ordered(!!id1_col_sym, levels = levels(otu_table[[id_col]])),
      !!id2_col_sym := ordered(!!id2_col_sym, levels = levels(otu_table[[id_col]]))
    )
  }

  # Ensure that the first sequence is always the one with the smaller id
  match_table <- match_table |>
    dplyr::filter(!!id1_col_sym != !!id2_col_sym) |>
    dplyr::mutate(
      seq_id3 = !!id1_col_sym,
      nread3 = nread1,
      swap = !!id1_col_sym > !!id2_col_sym,
      nread1 = ifelse(swap, nread2, nread1),
      nread2 = ifelse(swap, nread3, nread2),
      !!id1_col_sym := ifelse(swap, !!id2_col_sym, !!id1_col_sym),
      !!id2_col_sym := ifelse(swap, seq_id3, !!id2_col_sym)
    )

  # Actual LULU algorithm implemented in C++.
  # This is where it will fail if the ids are not sorted as required.
  out <- lulu_map_impl(
    match_id1 = match_table[[id1_col]],
    match_id2 = match_table[[id2_col]],
    match_nread1 = match_table$nread1,
    match_nread2 = match_table$nread2,
    match_dist = match_table$dist,
    seq_idx = otu_table[[id_col]],
    nread = otu_table$nread,
    max_dist = max_dist,
    min_abundance_ratio = min_abundance_ratio,
    min_cooccurrence_ratio = min_cooccurrence_ratio,
    use_mean_abundance_ratio = use_mean_abundance_ratio
  )

  # Map the result back to the original IDs
  if (is.ordered(otu_table[[id_col]])) {
    out <- dplyr::mutate(out, dplyr::across(everything(), as.character))
    if (id_is_int) {
      out <- dplyr::mutate(out, dplyr::across(everything(), as.integer))
    } else {
      out <- dplyr::rename(out, seq_id = seq_idx, lulu_id = lulu_idx)
    }
  }
  out
}

#' Apply LULU denoising to an OTU table
#'
#' @param lulu_map (`data.frame`) output of a call to `lulu_map()`
#' @param otu_table (`data.frame`) the long-format table to denoise. Must have
#' columns `seq_id` (`character`) xor `seq_idx` (`integer`), as well as
#' `nread` (`integer`). Additional columns should collectively define the
#' sample. Typically this includes the sample name, and sample-level metadata
#' are also fine to include, but sequence-level data are not. If no such values
#' are included, then the entire table is treated as a single sample.
#'
#' @returns a `data.frame` with the same columns as `otu_table`, but with rows
#' combined according to `lulu_map`.
#'
#' @export
lulu_table <- function(
    lulu_map,
    otu_table
) {
  checkmate::assert_data_frame(lulu_map)
  checkmate::assert_data_frame(otu_table)
  checkmate::assert(
    checkmate::check_integerish(otu_table$seq_idx, lower = 1L),
    checkmate::check_character(otu_table$seq_id)
  )
  if ("seq_id" %in% names(otu_table) && "seq_idx" %in% names(otu_table)) {
    stop("`otu_table` may have column `seq_id` or `seq_idx`, but not both")
  }
  id_as_int <- "seq_idx" %in% names(otu_table)
  id_col <- if (id_as_int) "seq_idx" else "seq_id"
  id_col_sym <- rlang::sym(id_col)
  lulu_col <- if (id_as_int) "lulu_idx" else "lulu_id"
  lulu_col_sym <- rlang::sym(lulu_col)
  checkmate::assert_names(
    names(otu_table),
    must.include = c("nread", id_col)
  )
  checkmate::assert_integerish(otu_table$nread, lower = 1L)
  checkmate::assert_names(
    names(lulu_map),
    must.include = c(id_col, lulu_col)
  )
  checkmate::assert_subset(otu_table[[id_col]], lulu_map[[id_col]])
  otu_table |>
    dplyr::left_join(lulu_map, by = id_col) |>
    dplyr::summarize(
      nread = sum(nread),
      .by = c(everything(), -nread, -!!id_col_sym)
    ) |>
    dplyr::rename(!!id_col_sym := !!rlang::sym(lulu_col))
}

#' Compute pairwise distances for LULU
#'
#' This is a pipeline step which calculates pairwise distances between
#' sequences in a single sample. It is intended to be used inside
#' `dplyr::summarize()` after grouping by sample.
#'
#' @param seqall_file (`character` string) name(s) of one or more gzipped FASTA
#' files containing the master list of sequences. Headers should be parsable as
#' integers, and should be the 1-based index of the sequence in the file.
#' @param seqall_index (`character` string) name of the index file for
#' `seqall_file`, as generated by [fastx_gz_index()]. Not supported for
#' multiple input files; in that case the files will be fully read into memory
#' one at a time.
#' @param seqtable (`data.frame`) the long-format table to denoise. Must have
#' columns `seq_idx` (`integer`) and `nread` (`integer`), where `seq_idx` refers
#' to the index of sequences in `seqall_file`. Additional columns
#' are ignored.
#' @param threshold (`numeric` scalar) maximum pairwise distance to report
#' between two sequences.
#' @param dist_config (`optimotu_dist_config`) configuration for distance
#' calculation, as returned by [optimotu::dist_config()] or its helpers.
#' @param parallel_config (`optimotu_parallel_config`) configuration for
#' parallel execution, as returned by [optimotu::parallel_config()] or its
#' helpers. Only the number of threads is used.
#' @param ... additional arguments for dependency tracking. These are ignored.
#' @export
#' @returns a `data.frame` with columns `seq_idx1`, `seq_idx2`, `dist`,
#' `align_length`, `n_gap`, `max_gap`, `nread1`, and `nread2`.
#' @keywords internal

lulu_distmx <- function(
  seqall_file,
  seqtable,
  seqall_index = NULL,
  threshold = lulu_max_dist(),
  dist_config = lulu_dist_config(),
  parallel_config = optimotu::parallel_concurrent(local_cpus()),
  ...
) {
  if (length(seqall_file) == 1L && !is.null(seqall_index)) {
    seqs <- fastx_gz_random_access_extract(
      infile = seqall_file,
      index = seqall_index,
      i = seqtable$seq_idx,
      outfile = withr::local_tempfile(fileext = ".fasta")
    )
  } else {
    seqs <- lapply(
      seq_along(seqall_file),
      function(i) {
        s <- Biostrings::readBStringSet(seqall_file[i])
        which_seq <- intersect(as.character(seqtable$seq_idx), names(s))
        s[which_seq]
      }
    )
    seqs <- do.call(c, seqs)
  }

  distmx <- optimotu::seq_distmx(
    seqs,
    threshold = optimotu::threshold_as_dist(threshold),
    dist_config = dist_config,
    parallel_config = parallel_config,
    details = "gapstats",
    span = "global",
    id_is_int = FALSE # this currently returns indices rather than parsed names
  )

  dplyr::transmute(
    distmx,
    seq_idx1 = as.integer(seq_id1),
    seq_idx2 = as.integer(seq_id2),
    dist = dist2,
    align_length,
    n_gap = n_insert + n_delete,
    max_gap = pmax(max_insert, max_delete)
  ) |>
    dplyr::left_join(
      dplyr::select(seqtable, seq_idx, nread1 = nread),
      by = c("seq_idx1" = "seq_idx")
    ) |>
    dplyr::left_join(
      dplyr::select(seqtable, seq_idx, nread2 = nread),
      by = c("seq_idx2" = "seq_idx")
    )
}
