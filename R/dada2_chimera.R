#' Generate a table of denovo bimeras
#'
#' This function performs the per-sample *de novo* chimera detection step of
#' `dada2::isBimeraDenovoTable()`, but does not do the consensus step. This
#' allows the function to be performed on smaller partitions of the sequence
#' table to reduce total memory usage. The results can then be combined using
#' `combine_bimera_denovo_tables()`. The implementation calls the internal API
#' of DADA2 and so it may be fragile to changes in the DADA2 package.
#'
#' @param seqtab (`matrix` or `data.frame`) input sequence table, as returned by
#' `dada2::makeSequenceTable()`, `make_long_sequence_table()`, or
#' `make_mapped_sequence_table()`
#' @param seqs (file name of a FASTA file, `character` vector,
#' [`DNAStringSet`][Biostrings::XStringSet-class], or `NULL`) sequences to use
#' in the table. If `NULL`, the sequences must be included in `seqtab`, either
#' as column names (for the `matrix` method) or as a column named `seq` (for the
#' `data.frame` method).
#' @param minFoldParentOverAbundance (`numeric(1)`) see `dada2::removeBimeraDenovo()`
#' @param minParentAbundance (`numeric(1)`) see `dada2::removeBimeraDenovo()`
#' @param allowOneOff (`logical(1)`) see `dada2::removeBimeraDenovo()`
#' @param minOneOffParentDistance (`numeric(1)`) see `dada2::removeBimeraDenovo()`
#' @param maxShift (`numeric(1)`) see `dada2::removeBimeraDenovo()`
#' @param multithread (`logical(1)`) if `TRUE`, use multiple threads for the
#' computation. If `FALSE`, use a single thread. If an integer, use that many
#' threads. If `NULL`, use the default number of threads.
#' @param ... additional arguments to be passed to the method
#' @return a `data.frame` with one of the columns `seq`, `seq_id`, or `seq_idx`,
#' as well as columns `nsam`, and `nflag`, where `seq`
#' is the full sequence, `seq_id` is a string identifier, and `seq_idx` is an
#' integer identifier, `nsam` is the number of samples in which the sequence
#' appears, and `nflag` is the number of samples in which the sequence is
#' flagged as a bimera.
#' @export
bimera_denovo_table <- function(
    seqtab,
    seqs = NULL,
    minFoldParentOverAbundance = 1.5,
    minParentAbundance = 2,
    allowOneOff = FALSE,
    minOneOffParentDistance = 4,
    maxShift = 16,
    multithread = FALSE,
    ...
) {
  UseMethod("bimera_denovo_table", seqtab)
}

#' @rdname bimera_denovo_table
#' @exportS3Method bimera_denovo_table matrix
bimera_denovo_table.matrix <- function(
    seqtab,
    seqs = colnames(seqtab),
    minFoldParentOverAbundance = 1.5,
    minParentAbundance = 2,
    allowOneOff = FALSE,
    minOneOffParentDistance = 4,
    maxShift = 16,
    multithread = FALSE,
    ...
) {
  if (is.null(seqs)) seqs <- colnames(seqtab)
  if (isTRUE(multithread)) {
    RcppParallel::setThreadOptions(numThreads = "auto")
  } else if (isFALSE(multithread)) {
    RcppParallel::setThreadOptions(numThreads = 1)
  } else {
    checkmate::assert_count(multithread, positive = TRUE)
    RcppParallel::setThreadOptions(numThreads = multithread)
  }
  dada2:::C_table_bimera2(
    mat = seqtab,
    seqs = seqs,
    min_fold = minFoldParentOverAbundance,
    min_abund = minParentAbundance,
    allow_one_off = allowOneOff,
    min_one_off_par_dist = minOneOffParentDistance,
    match = dada2::getDadaOpt("MATCH"),
    mismatch = dada2::getDadaOpt("MISMATCH"),
    gap_p = dada2::getDadaOpt("GAP_PENALTY"),
    max_shift = maxShift
  ) |>
    tibble::as_tibble() |>
    tibble::add_column(seq = seqs)
}

#' @rdname bimera_denovo_table
#' @exportS3Method bimera_denovo_table data.frame
bimera_denovo_table.data.frame <- function(
    seqtab,
    seqs = NULL,
    minFoldParentOverAbundance = 1.5,
    minParentAbundance = 2,
    allowOneOff = FALSE,
    minOneOffParentDistance = 4,
    maxShift = 16,
    multithread = FALSE,
    ...
) {
  # avoid R CMD check NOTE about undeclared global variables
  nread <- nflag <- nsam <- NULL

  seq_col <- intersect(c("seq", "seq_id", "seq_idx"), names(seqtab))
  if (length(seq_col) == 0) {
    stop("seqtab must have at least one of columns 'seq', 'seq_id', or 'seq_idx'")
  }
  seq_col <- seq_col[1]
  if (seq_col != "seq") {
    if (length(seqs) == 1 && file.exists(seqs)) {
      seqs <- Biostrings::readDNAStringSet(seqs)
    }
    if (methods::is(seqs, "XStringSet")) seqs <- as.character(seqs)
    if (!is.character(seqs)) {
      stop("'seqs' must be a filename of a fasta file, an XStringSet, or a character")
    }
    if (identical(seq_col, "seq_id") && !rlang::is_named(seqs)) {
      stop("'seqs' must be named if sequences are identified by 'seq_id' in 'seqtab'")
      checkmate::assert_subset(names(seqs), seqtab$seq_id)
    }
  }
  n_asv <- dplyr::n_distinct(seqtab[[seq_col]])
  n_sample <- dplyr::n_distinct(seqtab$sample)
  # max seqtable size for one partition is 1 Gb (== 2^30 bytes)
  # (not including sequences)
  # R integers are 32 bit (== 4 bytes)
  # If there are more than 250M ASVs in a single sample, then the matrix ends
  # up larger (but at that point the size of the sequences themselves is a bigger problem!)
  n_partition <- ceiling(n_asv*n_sample*4/2^30)
  sample_splits <-
    split(unique(seqtab$sample), rep(seq_len(n_partition), length.out = n_sample))
  out <- list()
  for (s in sample_splits) {
    m <- dplyr::filter(seqtab, sample %in% s) |>
      tidyr::pivot_wider(
        names_from = all_of(seq_col),
        values_from = nread,
        values_fill = list(nread = 0L)
      ) |>
      tibble::column_to_rownames("sample") |>
      as.matrix()

    switch(seq_col,
           seq_id = colnames(m) <- seqs[colnames(m)],
           seq_idx = colnames(m) <- seqs[as.integer(colnames(m))]
    )

    out_m <- bimera_denovo_table.matrix(
      seqtab = m,
      minFoldParentOverAbundance = minFoldParentOverAbundance,
      minParentAbundance = minParentAbundance,
      allowOneOff = allowOneOff,
      minOneOffParentDistance = minOneOffParentDistance,
      maxShift = maxShift,
      multithread = multithread
    )
    switch(
      seq_col,
      seq_idx = out_m$seq <- match(out_m$seq, seqs),
      seq_id = out_m$seq <- names(seqs)[match(out_m$seq, seqs)]
    )
    names(out_m)[3] <- "seq_idx"

    out <- c(
      out,
      list(out_m)
    )
  }

  dplyr::bind_rows(out) |>
    dplyr::summarize(nflag = sum(nflag), nsam = sum(nsam), .by = all_of(seq_col))
}

#' Combine bimera tables from `bimera_denovo_table()`
#' @param bimdf (`data.frame`) output of one or more calls to
#' `bimera_denovo_table()`
#' @param minSampleFraction (`numeric(1)`) see `dada2::isBimeraDenovoTable()`
#' @param ignoreNNegatives (`integer(1)`) see `dada2::isBimeraDenovoTable()`
#' @param verbose (`logical(1)`) see `dada2::isBimeraDenovoTable()`
#' @return if `bimdf` was generated from a
#' [dada2-style sequence table][dada2::makeSequenceTable()] or a
#' [long sequence table][make_long_sequence_table()] with an explicit `seq` column,
#' then a character vector of sequences that are bimeras. If `bimdf` was a
#' [mapped sequence table][make_mapped_sequence_table()] with a `seq_id` or
#' `seq_idx` column, then values corresponding to the `seq_id` or `seq_idx` of
#' the bimeras.
#' @export
combine_bimera_denovo_tables <- function(
    bimdf,
    minSampleFraction = 0.9,
    ignoreNNegatives = 1L,
    verbose = FALSE
) {
  seq_col <- intersect(c("seq", "seq_id", "seq_idx"), names(bimdf))
  if (length(seq_col) == 0) {
    stop("bimdf must have at least one of columns 'seq', 'seq_id', or 'seq_idx'")
  }
  seq_col <- seq_col[1]

  bimdf <- dplyr::summarize(
    bimdf,
    dplyr::across(everything(), sum),
    .by = any_of(seq_col)
  )
  ## This snippet modified from DADA2
  bims.out <- with(
    bimdf,
    nflag >= nsam | (nflag > 0 & nflag >= (nsam - ignoreNNegatives) * minSampleFraction)
  )
  if (verbose)
    message("Identified ", sum(bims.out), " bimeras out of ",
            length(bims.out), " input sequences.")
  ## end snippet from DADA2
  bimdf[[seq_col]][bims.out]
}
