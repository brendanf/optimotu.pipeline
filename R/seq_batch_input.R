#' Internal helpers for indexed / sequential sequence batching
#'
#' Shared building blocks for wrappers that partition sequence inputs across
#' parallel external-tool workers. Tempfiles are created only when
#' materializing per-chunk FASTA paths for subprocesses; in-memory conversion
#' does not use disk. Those tempfiles use `withr::local_tempfile(...,
#' .local_envir = local_envir)` so they persist until the caller exits (e.g.
#' [hmmalign()]).
#'
#' @name seq_batch_input
#' @keywords internal
NULL

#' Whether `x` is a character vector of existing `.fqi` index paths
#'
#' @noRd
seq_batch_is_fqi_path_set <- function(x) {
  if (!is.character(x) || anyNA(x)) {
    return(FALSE)
  }
  if (!all(file.exists(x))) {
    return(FALSE)
  }
  all(grepl("\\.fqi$", x, ignore.case = TRUE))
}

#' Total FASTA record count across paths (for `seq_idx` validation)
#'
#' @noRd
seq_batch_n_fasta_path_records <- function(infiles) {
  sum(vapply(
    infiles,
    function(p) length(Biostrings::fasta.seqlengths(p)),
    integer(1L)
  ))
}

#' Materialize per-chunk temp FASTA paths for up to `ncpu` parallel workers
#'
#' Returns `character(0)` when there are no sequences to materialize (no
#' records, or an empty resolved `seq_idx` stream). Otherwise returns one path
#' per non-empty chunk (length at most `min(ncpu, m)` for `m` selected records),
#' or the passthrough vector from the special case below.
#'
#' Includes a special case for `ncpu == length(seqs)` and seq_idx is `NULL`,
#' where the input files are used directly, avoiding the need for temporary
#' files. This will likely cause downstream errors with fastq inputs.
#' Gzipped FASTA files also use this special case, but do generate temporary
#' files.
#'
#' @param seqs (`character` vector, [`XStringSet`][Biostrings::XStringSet],
#'   `data.frame`, or `fastqindexr_index`)
#'   Source sequences. `character` vectors may be literal sequences (in which
#'   case they should typically be named) or file paths; files may be either
#'   FASTA (possibly gzipped) or `.fqi` indexes for
#'   [`fastqindexr::read_fqi_index()`].
#' @param files (`character` vector) optional per-file paths overriding those
#'   stored in the index, if `seqs` is a
#'   [`fastqindexr_index`][fastqindexr::create_index()] object or `.fqi`
#'   path(s). Useful after moving inputs or for `targets` dependency tracking.
#' @param seq_idx (`integer` vector) optional 1-based indices into the logical
#'   sequence stream (`NULL` means all sequences in order). Applies after
#'   concatenating multiple FASTA inputs, and supports duplicates and
#'   reordering.
#' @param ncpu (`integer`) maximum number of parallel processes to use. The
#'   selected sequences are split into up to `ncpu` contiguous chunks of nearly
#'   equal size (fewer when there are fewer sequences than `ncpu`).
#' @param local_envir (`environment`) the caller frame that must keep tempfiles
#'   alive until chunk paths are consumed (defaults to [parent.frame()]).
#' @return `character` vector of per-chunk FASTA paths, or `character(0)` if
#'   there is nothing to extract.
#' @noRd
seq_batch_make_chunk_files <- function(
  seqs,
  files,
  seq_idx = NULL,
  ncpu = local_cpus(),
  local_envir = parent.frame()
) {
  index_obj <- NULL
  if (seq_batch_is_fqi_path_set(seqs)) {
    seqs <- fastqindexr::read_fqi_index(
      fqi_path = seqs,
      files = files,
      type = "auto"
    )
    # fall through to next block
  }
  if (inherits(seqs, "fastqindexr_index")) {
    index_obj <- seqs
    seqs <- if (!is.null(files)) files else index_obj$files
  }
  if (checkmate::test_file_exists(seqs, "r")) {
    if (length(seqs) == ncpu && is.null(seq_idx)) {
      outfiles <- seqs
      gzip_files <- which(endsWith(outfiles, ".gz"))
      gzip_processes <- vector("list", length(gzip_files))
      for (i in gzip_files) {
        outfiles[i] <- withr::local_tempfile(
          fileext = ".fasta",
          .local_envir = local_envir
        )
        gzip_processes[[i]] <- processx::process$new(
          command = "zcat",
          args = seqs[i],
          stdout = outfiles[i]
        )
      }
      gzip_return <- 0L
      for (i in seq_along(gzip_processes)) {
        gzip_return <- union(
          gzip_processes[[i]]$wait()$get_exit_status(),
          gzip_return
        )
      }
      stopifnot(identical(gzip_return, 0L))
      is_empty <- vapply(outfiles, file.size, 0) == 0
      outfiles <- outfiles[!is_empty]
      return(outfiles)
    }
    n_rec <- seq_batch_n_fasta_path_records(seqs)
    if (n_rec == 0L) {
      return(character())
    }
    idx_full <- resolve_linear_seq_idx(n_rec, seq_idx)
    chunks <- partition_vector_equal_ncpu(idx_full, ncpu)
    nonempty <- vapply(chunks, length, integer(1L)) > 0L
    chunks_ne <- chunks[nonempty]
    if (length(chunks_ne) == 0L) {
      return(character())
    }
    out_ne <- replicate(
      length(chunks_ne),
      withr::local_tempfile(
        fileext = ".fasta",
        .local_envir = local_envir
      )
    )
    fastqindexr::extract_sequences_to_file(
      index = index_obj,
      seq_idx = chunks_ne,
      file = seqs,
      outfile = out_ne,
      type = "fasta",
      input_type = "auto",
      append = FALSE,
      compress = FALSE,
      collapse_sequence_lines = TRUE
    )
    return(out_ne)
  }
  n_rec <- sequence_size(seqs)
  if (n_rec == 0L) {
    return(character())
  }
  checkmate::assert_count(n_rec)
  idx_full <- resolve_linear_seq_idx(n_rec, seq_idx)
  chunks_i <- partition_vector_equal_ncpu(idx_full, ncpu)
  nonempty <- vapply(chunks_i, length, integer(1L)) > 0L
  chunks_i <- chunks_i[nonempty]
  if (length(chunks_i) == 0L) {
    return(character())
  }
  vapply(
    chunks_i,
    function(ii) {
      tseq <- withr::local_tempfile(
        fileext = ".fasta",
        .local_envir = local_envir
      )
      write_sequence(select_sequence(seqs, ii), tseq)
      tseq
    },
    character(1)
  )
}
