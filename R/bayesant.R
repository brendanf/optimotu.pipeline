#' Identify sequences using BayesANT
#' @param query (`character` vector) the sequences to identify
#' @param model (`character` file name or `BayesANT` model) the BayesANT model
#' to use.  If a file name, it must be in .rds, .qs, or .qs2 format.
#' @param ncpu (`integer`) the number of CPU cores to use
#' @param id_is_int (`logical`) if `TRUE`, parse the sequence IDs as integers
#' @param n_top_taxa (`integer`) the number of top taxa to return
#' @param min_prob (`numeric`) the minimum probability to return
#' @return a `data.frame` with columns `seq_id` (or `seq_idx` if `id_is_int` is
#' `TRUE`), `rank`, `parent_taxonomy`, `taxon`, and `prob`, where `seq_id`
#' (`seq_idx`) is the ID of a sequence from `query`, `rank` is the taxonomic
#' rank at which the prediction is being made, `parent_taxonomy` is the parent
#' taxonomic unit of the taxon at the previous rank, `taxon` is the name of the
#' predicted taxon, and `prob` is the probability of the predicted taxon being
#' the correct classification of the sequence from `query` at the current rank.
#' @export
bayesant <- function(
  query,
  model,
  ncpu = local_cpus(),
  id_is_int = FALSE,
  n_top_taxa = 20,
  min_prob = 0.01
) {
  # avoid R CMD check NOTE
  taxon <- i <- leaf_prob <- prob <- seq_id <- NULL
  if (is.character(model) && file.exists(model)) {
    if (endsWith(model, ".rds")) {
      model <- readRDS(model)
    } else if (endsWith(model, ".qs")) {
      if (!requireNamespace("qs")) {
        stop(
          "qs package is required but not installed. Please install it",
          " using `install.packages('qs')`."
        )
      }
      model <- qs::qread(model, nthreads = ncpu)
    } else if (endsWith(model, ".qs2")) {
      if (!requireNamespace("qs2")) {
        stop(
          "qs2 package is required but not installed. Please install it",
          " using `install.packages('qs2')`."
        )
      }
      model <- qs2::qs_read(model, nthreads = ncpu)
    } else if (endsWith(model, ".qdata")) {
      stop("Model file must be in .rds, .qs, or .qs2 format.")
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
  if (checkmate::test_file_exists(query, access = "r")) {
    query <- BayesANT::read.BayesANT.testDNA(query)
  }
  if (methods::is(query, "DNAStringSet")) {
    query <- as.character(query)
  }
  if (is.data.frame(query)) {
    query <- `names<-`(
      query[[find_seq_col(query)]],
      query[[find_name_col(query)]]
    )
  }
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
