#' Identify sequences using BayesANT
#' @param query (`character` vector, `data.frame`, [`DNAStringSet`][Biostrings::XStringSet-class],
#' or `character` file name) the sequences to identify. If a file name, it
#' should point to a FASTA file (possibly gzipped).
#' @param model (`character` file name or `BayesANT` model) the BayesANT model
#' to use.  If a file name, it must be in .rds or .qs2 format.  Legacy `.qs`
#' files are no longer supported; convert to `.qs2` or `.rds`.
#' @param ncpu (`integer`) the number of CPU cores to use
#' @param id_is_int (`logical`) if `TRUE`, parse the sequence IDs as integers
#' @param n_top_taxa (`integer`) the number of top taxa to return
#' @param min_prob (`numeric`) the minimum probability to return
#' @param file (`character` vector) optional per-file paths overriding those
#' stored in the index, if `query` is a
#' [`fastqindexr_index`][fastqindexr::create_index()] object or `.fqi` /
#' `.qs2` path(s). Useful after moving inputs or for `targets` dependency
#' tracking.
#' @param seq_idx (`integer` vector) optional 1-based indices into the logical
#' sequence stream (`NULL` means all sequences in order). Applies after
#' concatenating multiple FASTA inputs, and supports duplicates and
#' reordering.
#' @param ... ignored; reserved for `targets` dependency tracking (e.g.
#'   `hash = seqbatch_hash`) without changing behavior.
#' @return a `data.frame` with columns `seq_id` (or `seq_idx` if `id_is_int` is
#' `TRUE`), `rank`, `parent_taxonomy`, `taxon`, and `prob`, where `seq_id`
#' (`seq_idx`) is the ID of a sequence from `query`, `rank` is the taxonomic
#' rank at which the prediction is being made, `parent_taxonomy` is the parent
#' taxonomic unit of the taxon at the previous rank, `taxon` is the name of the
#' predicted taxon, and `prob` is the probability of the predicted taxon being
#' the correct classification of the sequence from `query` at the current rank.
#' If `query` is empty, returns an empty table with the same column schema.
#' @export
bayesant <- function(
  query,
  model,
  ncpu = local_cpus(),
  id_is_int = FALSE,
  n_top_taxa = 20,
  min_prob = 0.01,
  file = NULL,
  seq_idx = NULL,
  ...
) {
  checkmate::assert_count(ncpu)
  checkmate::assert_flag(id_is_int)
  checkmate::assert_count(n_top_taxa)
  checkmate::assert_number(min_prob, lower = 0, upper = 1)
  checkmate::assert_character(file, null.ok = TRUE)
  checkmate::assert_integerish(seq_idx, null.ok = TRUE)

  indexed_like <- inherits(query, "fastqindexr_index") ||
    seq_batch_is_fqi_path_set(query)
  if (!is.null(file) && !indexed_like) {
    stop(
      "`file` is only valid when `query` is a fastqindexr_index or ",
      ".fqi/.qs2 paths.",
      call. = FALSE
    )
  }
  # avoid R CMD check NOTE
  taxon <- i <- leaf_prob <- prob <- seq_id <- NULL
  if (is.character(model) && file.exists(model)) {
    if (endsWith(model, ".rds")) {
      model <- readRDS(model)
    } else if (endsWith(model, ".qs")) {
      .stop_qs_deprecated(".qs BayesANT model files")
    } else if (endsWith(model, ".qs2")) {
      if (!requireNamespace("qs2")) {
        stop(
          "qs2 package is required but not installed. Please install it",
          " using `install.packages('qs2')`."
        )
      }
      model <- qs2::qs_read(model, nthreads = ncpu)
    } else if (endsWith(model, ".qdata")) {
      stop("Model file must be in .rds or .qs2 format.")
    }
  } else if (is.character(model) && !file.exists(model)) {
    stop(
      "Model file ",
      model,
      " does not exist.",
      " Please provide a valid file path."
    )
  } else if (!inherits(model, "BayesANT")) {
    stop("Model must be a BayesANT model or a valid file path.")
  }
  query <- seq_batch_character(query, file, seq_idx)
  checkmate::assert_character(
    query,
    min.chars = 1,
    names = "unique",
    any.missing = FALSE
  )
  out <- stats::predict(
    model,
    query,
    return_probs = TRUE,
    cores = ncpu,
    n_top_taxa = n_top_taxa
  )$top_n_probs

  # result is a list of data frames, one for each query sequence
  # each data frame has columns named after the taxonomic ranks, plus
  # `leaf_prob`
  out <- dplyr::bind_rows(out, .id = "seq_id") |>
    tibble::rowid_to_column("i") |>
    tidyr::pivot_longer(
      cols = any_of(tax_ranks()),
      names_to = "rank",
      values_to = "taxon",
      names_transform = rank2factor
    ) |>
    dplyr::mutate(
      parent_taxonomy = purrr::accumulate(taxon, paste, sep = ",") |>
        dplyr::lag(),
      .by = i
    ) |>
    dplyr::summarize(
      prob = sum(leaf_prob),
      .by = c("seq_id", "rank", "parent_taxonomy", "taxon")
    ) |>
    dplyr::filter(prob >= min_prob)

  if (id_is_int) {
    out <- out |>
      dplyr::mutate(seq_idx = as.integer(seq_id), .keep = "unused", .before = 1)
  }

  out
}

#' Inject code to read BayesANT file from disk or target
#' @noRd
read_bayesant_model <- function(model = bayesant_model()) {
  if (is.null(model)) {
    quote(bayesant_model)
  } else {
    ext <- tools::file_ext(model) |>
      tolower()
    if (ext == "rds") {
      substitute(readRDS(model), list(model = model))
    } else if (ext == "qs") {
      .stop_qs_deprecated(".qs bayesant_model files")
    } else if (ext == "qs2") {
      substitute(qs2::qs_read(model), list(model = model))
    } else {
      stop("Unsupported file type '", ext, "' for bayesant_model")
    }
  }
}
