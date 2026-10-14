#' Align a batch of sequences to an amplicon CM or HMM
#'
#' Writes one aligned FASTA for a (possibly subset) sequence batch. For
#' `model_type = "HMM"`, delegates to [hmmalign()] with A2M output. For
#' `model_type = "CM"`, runs Infernal `cmalign`, keeps consensus (RF) columns
#' via [consensus_columns()], and writes the result. Optionally also writes a
#' persistent Infernal `--sfile` (CM only) and returns both paths.
#'
#' This is the shared worker for batched model-alignment targets (query ASVs,
#' outgroup references, supplemental ASVs). It does not wrap `nhmmer` or CM
#' filter-only jobs that discard the alignment.
#'
#' @param seqs ([`XStringSet`][Biostrings::XStringSet-class], `character`
#'   vector, `data.frame`, or
#'   [`fastqindexr_index`][fastqindexr::create_index()] object) query
#'   sequences. A `character` vector may contain literal sequences (named) or
#'   file paths to FASTA / `.fqi` / `.qs2` indexes (same rules as
#'   [hmmalign()]).
#' @param model (`character` string) path to a CM or HMM file
#' @param outfile (`character` string) path for the aligned FASTA (gzip if
#'   `.gz`)
#' @param model_type (`character` string) `"CM"` or `"HMM"`
#' @param files (`NULL` or `character`) when `seqs` is a
#'   `fastqindexr_index` or `.fqi` paths, optional per-file paths overriding
#'   those stored in the index (same length as the index file list)
#' @param seq_idx (`NULL` or `integer`) optional 1-based indices into the
#'   logical sequence stream (`NULL` means all sequences in order)
#' @param sfile (`NULL` or `character` string) if non-`NULL`, write Infernal
#'   `--sfile` scores here (CM only). Returns `c(outfile, sfile)`. Errors if
#'   set with `model_type = "HMM"`.
#' @param ncpu (`integer`) CM: threads for `cmalign`. HMM: max parallel
#'   `hmmalign` processes (see [hmmalign()])
#' @param hmmalign (`character`) path to the `hmmalign` executable
#' @param ... ignored; reserved for `targets` dependency tracking (e.g. file
#'   hashes) without changing behavior
#' @return `character`: `outfile`, or `c(outfile, sfile)` when `sfile` is set
#' @export
align_fasta_batch_to_amplicon_model <- function(
  seqs,
  model,
  outfile,
  model_type = c("CM", "HMM"),
  files = NULL,
  seq_idx = NULL,
  sfile = NULL,
  ncpu = local_cpus(),
  hmmalign = find_hmmalign(),
  ...
) {
  model_type <- match.arg(model_type)
  checkmate::assert_string(model)
  checkmate::assert_file_exists(model, access = "r")
  checkmate::assert_string(outfile)
  checkmate::assert_string(sfile, null.ok = TRUE)
  checkmate::assert_count(ncpu, positive = TRUE)
  if (!is.null(sfile) && identical(model_type, "HMM")) {
    stop(
      "`sfile` is only valid for model_type = \"CM\" ",
      "(HMM scores use nhmmer, not hmmalign).",
      call. = FALSE
    )
  }

  if (identical(model_type, "HMM")) {
    # Parameter name shadows the R function; call via namespace.
    return(optimotu.pipeline::hmmalign(
      seqs = seqs,
      hmm = model,
      outfile = outfile,
      outformat = "A2M",
      files = files,
      seq_idx = seq_idx,
      ncpu = ncpu,
      hmmalign = hmmalign,
      ...
    ))
  }

  # CM path: materialize one FASTA, cmalign, consensus columns, write.
  compress <- endsWith(outfile, ".gz")
  tmp_parent <- environment()
  chunk_files <- seq_batch_make_chunk_files(
    seqs = seqs,
    files = files,
    seq_idx = seq_idx,
    ncpu = 1L,
    local_envir = tmp_parent
  )
  if (length(chunk_files) == 0L) {
    ensure_directory(outfile)
    write_sequence(
      Biostrings::DNAStringSet(),
      outfile,
      compress = compress
    )
    if (!is.null(sfile)) {
      ensure_directory(sfile)
      writeLines(character(), sfile)
      return(c(outfile, sfile))
    }
    return(outfile)
  }

  if (!is.null(sfile)) {
    ensure_directory(sfile)
  }
  aln <- inferrnal::cmalign(
    model,
    chunk_files[[1L]],
    global = TRUE,
    notrunc = TRUE,
    dnaout = TRUE,
    cpu = ncpu,
    sfile = sfile
  )
  write_sequence(
    consensus_columns(aln),
    outfile,
    compress = compress
  )
  if (!is.null(sfile)) {
    c(outfile, sfile)
  } else {
    outfile
  }
}
