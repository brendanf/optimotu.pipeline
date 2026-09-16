#### taxonomic ranks ####

#' Get or set the taxonomic ranks to use in the pipeline
#' @return a `character` vector of taxonomic ranks
#' @export
tax_ranks <- function() {
  getOption(
    "optimotu.pipeline.tax_ranks",
    c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  )
}

#' @rdname tax_ranks
#' @return a `list` of `symbol`s representing the taxonomic ranks
#' @export
tax_rank_vars <- function() {
  rlang::syms(tax_ranks())
}

#' @rdname tax_ranks
#' @param value (`character` vector) the taxonomic ranks to use in the pipeline,
#' in order from most inclusive (e.g., kingdom) to least inclusive (e.g., species)
#' @return the previously set taxonomic ranks
#' @export
set_tax_ranks <- function(value) {
  options("optimotu.pipeline.tax_ranks" = value)[[
    "optimotu.pipeline.tax_ranks"
  ]]
}

#' Get or set the offset to use when converting between integers and taxonomic
#' ranks
#' @return an integer offset
#' @export
rank_offset <- function() {
  getOption("optimotu.pipeline.rank_offset", 0)
}

#' @rdname rank_offset
#' @param value (`integer`) the offset to use when converting between integers
#' and taxonomic ranks
#' @return the previously set offset
#' @export
set_rank_offset <- function(value) {
  options("optimotu.pipeline.rank_offset" = value)[[
    "optimotu.pipeline.rank_offset"
  ]]
}

#' Get or set the taxonomic ranks which are assumed to be "known"
#' @return a `character` vector of taxonomic ranks
#' @export
known_ranks <- function() {
  getOption("optimotu.pipeline.known_ranks", tax_ranks()[1])
}

#' @rdname known_ranks
#' @return a `list` of `symbol`s representing the known
#' @export
known_rank_vars <- function() {
  rlang::syms(known_ranks())
}

#' @rdname known_ranks
#' @param value (`character` vector) the taxonomic ranks to assume are "known"
#' @return the previously set known ranks
#' @export
set_known_ranks <- function(value) {
  options("optimotu.pipeline.known_ranks" = value)[[
    "optimotu.pipeline.known_ranks"
  ]]
}

#' Get or set the taxa which are assumed to be "known"
#' @return a `character` vector of taxa, with the same length as `known_ranks()`
#' @export
known_taxa <- function() {
  getOption("optimotu.pipeline.known_taxa", "Fungi")
}

#' @rdname known_taxa
#' @param value (`character` vector) the taxa to assume are "known"
#' @return the previously set known taxa
#' @export
set_known_taxa <- function(value) {
  options("optimotu.pipeline.known_taxa" = value)[[
    "optimotu.pipeline.known_taxa"
  ]]
}

#' @rdname known_ranks
#' @return a `character` vector of taxonomic ranks
#' @export
unknown_ranks <- function() {
  setdiff(tax_ranks(), known_ranks())
}

#' @rdname known_ranks
#' @return a `list` of `symbol`s representing the unknown ranks
#' @export
unknown_rank_vars <- function() {
  rlang::syms(unknown_ranks())
}

#' Shortcuts for accessing taxonomy options
#' @name tax_opts
#' @param ranks (`character` vector) the taxonomic ranks to use
NULL

#' @describeIn tax_opts Get the root taxonomic rank as a string
#' @export
root_rank <- function() {
  tax_ranks()[1]
}

#' @describeIn tax_opts Get the root taxonomic rank as a symbol
#' @export
root_rank_var <- function() {
  rlang::sym(root_rank())
}

#' @describeIn tax_opts Get the second taxonomic rank as a string
#' @export
second_rank <- function() {
  tax_ranks()[2]
}

#' @describeIn tax_opts Get the second taxonomic rank as a symbol
#' @export
second_rank_var <- function() {
  rlang::sym(second_rank())
}

#' @describeIn tax_opts Get the ingroup taxonomic rank as a string
#' @export
ingroup_rank <- function() {
  tax_ranks()[length(known_ranks())]
}

#' @describeIn tax_opts Get the ingroup taxonomic rank as a symbol
#' @export
ingroup_rank_var <- function() {
  rlang::sym(ingroup_rank())
}

#' @describeIn tax_opts Get the ingroup taxon as a string
#' @export
ingroup_taxon <- function() {
  known_taxa()[length(known_taxa())]
}

#' @describeIn tax_opts Get the tip taxonomic rank as a string
#' @export
tip_rank <- function(ranks = tax_ranks()) {
  ranks[length(ranks)]
}

#' @describeIn tax_opts Get the tip taxonomic rank as a symbol
#' @export
tip_rank_var <- function(ranks = tax_ranks()) {
  rlang::sym(tip_rank(ranks))
}

#' Convert taxonomic ranks from strings and integers to ordered factors
#'
#' The least inclusive rank is given the smallest value in the ordering,
#' so for instance `"species" < "kingdom"`.
#'
#' @param x the taxonomic ranks to convert
#' @param ranks (`character` vector) the full list of taxonomic ranks, in order
#' from most inclusive to least inclusive
#' @return an ordered factor of taxonomic ranks
#' @export
rank2factor <- function(x, ranks = tax_ranks()) {
  factor(x, levels = rev(ranks), ordered = TRUE)
}

#' @rdname rank2factor
#' @export
int2rankfactor <- function(x, ranks = tax_ranks()) {
  rank2factor(ranks[x + rank_offset()], ranks)
}

#' Get superordinate or subordinate taxonomic ranks
#'
#' @param x (`character` string) the taxonomic rank to find superordinate or
#' subordinate ranks for
#' @param ranks (`character` vector) the full list of taxonomic ranks, in order
#' from most inclusive to least inclusive
#' @return a `character` vector of superordinate or subordinate taxonomic ranks
#' @export
superranks <- function(x = tip_rank(), ranks = tax_ranks()) {
  ranks[rank2factor(ranks, ranks) > x]
}

#' @rdname superranks
#' @return a `list` of `symbol`s representing the superordinate ranks
#' @export
superrank_vars <- function(x = tip_rank(), ranks = tax_ranks()) {
  rlang::syms(superranks(x, ranks))
}

#' @rdname superranks
#' @export
subranks <- function(x = ingroup_rank(), ranks = tax_ranks()) {
  ranks[rank2factor(ranks, ranks) < x]
}

#' @rdname superranks
#' @return a `list` of `symbol`s representing the subordinate ranks
#' @export
subrank_vars <- function(x = ingroup_rank(), ranks = tax_ranks()) {
  rlang::syms(subranks(x, ranks))
}

#' Convert wide taxonomy/probability tables to long format
#'
#' @param taxonomy (`data.frame`) wide table with one taxon column per rank and
#' a `seq_id` column.
#' @param probability (`data.frame`) wide table with one probability column per
#' rank and a `seq_id` column.
#' @param ranks (`character` vector) taxonomic ranks to pivot.
#' @return `tibble` with columns `seq_id`, `rank`, `parent_taxonomy`, `taxon`,
#' and `prob`.
#' @export
tax_table_wide_to_long <- function(taxonomy, probability, ranks = tax_ranks()) {
  # avoid R CMD check note for undefined global variables due to NSE
  seq_id <- parent_taxonomy <- taxon <- prob <- NULL
  checkmate::assert_data_frame(taxonomy)
  checkmate::assert_data_frame(probability)
  checkmate::assert_character(ranks, min.len = 1)
  dplyr::full_join(
    tidyr::pivot_longer(
      taxonomy,
      dplyr::all_of(ranks),
      names_to = "rank",
      values_to = "taxon"
    ),
    tidyr::pivot_longer(
      probability,
      dplyr::all_of(ranks),
      names_to = "rank",
      values_to = "prob"
    ),
    by = c("seq_id", "rank")
  ) |>
    dplyr::mutate(
      rank = factor(rank, levels = ranks, ordered = TRUE),
      parent_taxonomy = NA_character_
    ) |>
    dplyr::select(seq_id, rank, parent_taxonomy, taxon, prob)
}

#' Combine tip classifications to build a full PROTAX taxonomy
#' @param ... (`character` vectors) the classifications to combine; these should
#' be comma-delimited classifications from most inclusive rank to least
#' inclusive.
#' @return a `data.frame` with columns `taxon_id`, `parent_id`, `rank`,
#' `classification`, and `prior`, as required to write a Protax "`taxonomy`"
#' file
#' @export
build_taxonomy <- function(...) {
  # avoid R CMD check NOTE for undeclared global variables
  classification <- parent <- prior <- taxon_id <- parent_id <- NULL

  tax <- tibble::tibble(
    classification = union(...),
    rank = ifelse(
      classification == "root",
      0L,
      stringr::str_count(classification, stringr::fixed(",")) + 1L
    ),
    parent = ifelse(
      rank <= 1,
      "root",
      sub(",[^,]+$", "", classification)
    )
  )
  tax <- split(tax, tax$rank)
  tax[[1]]$taxon_id <- 0L
  tax[[1]]$prior <- 1
  tax[[1]]$parent_id <- 0L
  tax[[2]]$taxon_id <- 1L:2L
  tax[[2]]$prior <- c(0.99, 0.01)

  for (r in length(tax):3L) {
    if (r == length(tax)) {
      tax[[r]]$prior <- 0.99 / nrow(tax[[r]])
    } else {
      tax[[r]]$prior <- NULL
      tax[[r]] <- dplyr::left_join(
        tax[[r]],
        dplyr::group_by(tax[[r + 1]], parent) |>
          dplyr::summarise(prior = sum(prior)) |>
          dplyr::rename(classification = parent),
        by = "classification"
      )
    }
  }
  for (r in 2L:length(tax)) {
    tax[[r]]$taxon_id <- seq_len(nrow(tax[[r]])) + max(tax[[r - 1L]]$taxon_id)
    tax[[r]]$parent_id <- NULL
    tax[[r]] <- dplyr::left_join(
      tax[[r]],
      dplyr::select(
        tax[[r - 1]],
        parent = classification,
        parent_id = taxon_id
      ),
      by = "parent"
    )
  }
  dplyr::bind_rows(tax) |>
    dplyr::select(taxon_id, parent_id, rank, classification, prior)
}

#' Format a classification for use as a Sintax reference DB
#' @param s (`character` vector) the classification to format; should be
#' comma-delimited from most inclusive rank to least inclusive
sintax_format <- function(s) {
  s <- sub(",", ";p:", s, fixed = TRUE)
  s <- sub(",", ";c:", s, fixed = TRUE)
  s <- sub(",", ";o:", s, fixed = TRUE)
  s <- sub(",", ";f:", s, fixed = TRUE)
  s <- sub(",", ";g:", s, fixed = TRUE)
  s <- sub(",", ";s:", s, fixed = TRUE)
  s <- chartr(";", ",", s)
  paste0("tax=d:", s, ";")
}

#' @rdname build_taxonomy
#' @export
#' @importFrom dplyr one_of
#' @importFrom rlang `:=`
build_taxonomy_new <- function(...) {
  # avoid R CMD check NOTE for undeclared global variables
  classification <- parent_classification <- prior <- taxon_id <- parent_id <-
    n <- NULL
  tax <- tibble::tibble(classification = unique(c(...))) |>
    tidyr::separate_wider_delim(
      1,
      delim = ",",
      names = tax_ranks(),
      too_few = "align_start"
    ) |>
    dplyr::mutate(
      !!root_rank() := ifelse(
        !!root_rank() == "root",
        NA_character_,
        !!root_rank()
      )
    )

  # remove all taxa with children
  for (rank in rev(seq_along(tax_ranks()))) {
    rankname <- tax_rank_vars()[rank + 1]
    tax <- dplyr::filter(
      tax,
      if (any(!is.na(!!rankname))) !is.na(!!rankname) else TRUE,
      .by = all_of(tax_ranks()[1:rank])
    )
  }
  tax$n <- 1L
  tax_out <- list()
  tax_out[[1]] <- tibble::tibble(
    taxon_id = 0L,
    parent_id = 0L,
    rank = 0L,
    classification = "root",
    prior = 1.0
  )
  for (i in seq_along(tax_ranks())) {
    rankname <- tax_rank_vars()[i]
    tax_i <- dplyr::select(tax, one_of(tax_ranks()[1:i]), n)
    tax_i <- tax_i[!is.na(tax_i[[tax_ranks()[i]]]), ]
    tax_i <- dplyr::summarize(tax_i, n = sum(n), .by = one_of(tax_ranks()[1:i]))
    if (i == 1) {
      tax_i$parent_classification <- "root"
    } else {
      tax_i$parent_classification <-
        apply(tax_i[seq_len(i - 1)], 1, paste, collapse = ",", simplify = TRUE)
    }
    tax_i <- dplyr::left_join(
      tax_i,
      dplyr::select(
        tax_out[[i]],
        parent_classification = classification,
        parent_id = taxon_id,
        prior
      ),
      by = "parent_classification"
    )
    # TODO: allow to specify priors
    # TODO: unknown species?
    tax_i <- dplyr::mutate(tax_i, prior = prior * n / sum(n), .by = parent_id)
    if (i > 1L) {
      tax_i <- dplyr::mutate(
        tax_i,
        classification = paste(parent_classification, {{ rankname }}, sep = ",")
      )
    } else {
      tax_i <- dplyr::mutate(
        tax_i,
        classification = {{ rankname }}
      ) |>
        dplyr::arrange(classification == "nonFungi")
    }
    tax_i$taxon_id <- seq_len(nrow(tax_i)) + max(tax_out[[i]]$taxon_id)
    tax_i$rank <- i
    tax_out[[i + 1]] <-
      dplyr::select(tax_i, taxon_id, parent_id, rank, classification, prior)
  }
  dplyr::bind_rows(tax_out)
}

#' Truncate classification(s) at a given rank
#' @param s (`character` vector`) the classification(s) to truncate
#' @param rank (`integer`) the rank to truncate at
#' @return a `character` vector of truncated classifications
#' @export
truncate_taxonomy <- function(s, rank) {
  regex <- paste0("(^([^,]+,){", rank - 1L, "}[^,]+).*")
  out <- gsub(regex, "\\1", s)
  out[!grepl(regex, s)] <- NA_character_
  out
}

#' Remove mycobank numbers from genus and species names
#' @param taxon (`character`) the taxon name to remove the mycobank number from
#' @return the taxon name with the mycobank number removed
#' @export
remove_mycobank_number <- function(taxon) {
  dplyr::if_else(startsWith(taxon, "pseudo"), taxon, sub("_[0-9]+$", "", taxon))
}

#### classifier cascade helpers ####

#' Taxon labels that do not count as outgroup for best_hit screening
#'
#' Shared between the reclassify screen and
#' `final_asv_known_outgroup` construction.
#'
#' @return `character` vector of ignored taxon labels
#' @export
best_hit_ingroup_exceptions <- function() {
  c("unspecified", "Eukaryota_kgd_Incertae_sedis", "None")
}

#' Ranks a classifier actually assigns
#'
#' @param scope_ranks (`character`) ranks in the classifier's closed scope
#' @param ranks (`character`) full rank list
#' @return `character` vector of assigned ranks
#' @export
assigned_ranks <- function(scope_ranks, ranks = tax_ranks()) {
  checkmate::assert_character(scope_ranks, any.missing = FALSE)
  checkmate::assert_character(ranks, any.missing = FALSE, min.len = 1)
  setdiff(ranks, scope_ranks)
}

#' Sequence IDs that fail the primary classifier confidence / scope screen
#'
#' An ID is included when either:
#' 1. There is no known (non-`NA`, non-`"unk"`) taxon at `rank` with
#'    `prob >= min_prob`, or
#' 2. `screen` is supplied and the best-hit taxon at any `scope_ranks` is
#'    present and not equal to the corresponding `scope_taxa` (ignoring
#'    [best_hit_ingroup_exceptions()]).
#'
#' @param tax_prob (`data.frame`) long taxonomy/probability table from the
#'   primary classifier
#' @param rank (`character`) rank used for the probability gate
#' @param min_prob (`numeric`) minimum known-taxon probability
#' @param all_ids vector of all sequence IDs in the batch
#' @param id_col (`character`) ID column name (`"seq_idx"` or `"seq_id"`)
#' @param screen (`data.frame` or `NULL`) best-hit table with `id_col` and
#'   one column per `scope_ranks`
#' @param scope_ranks (`character`) closed-scope ranks for the primary
#' @param scope_taxa (`character`) expected taxa at `scope_ranks`
#' @return vector of IDs needing reclassification (same type as `all_ids`)
#' @export
ids_needing_reclassify <- function(
  tax_prob,
  rank,
  min_prob,
  all_ids,
  id_col = "seq_idx",
  screen = NULL,
  scope_ranks = character(),
  scope_taxa = character()
) {
  checkmate::assert_data_frame(tax_prob)
  checkmate::assert_string(rank)
  checkmate::assert_number(min_prob, lower = 0, upper = 1)
  checkmate::assert_string(id_col)
  checkmate::assert_character(scope_ranks, any.missing = FALSE)
  checkmate::assert_character(scope_taxa, any.missing = FALSE)
  if (length(scope_ranks) != length(scope_taxa)) {
    stop("`scope_ranks` and `scope_taxa` must have the same length.")
  }
  checkmate::assert_names(
    names(tax_prob),
    must.include = c(id_col, "rank", "taxon", "prob")
  )

  ids <- unique(all_ids)
  if (length(ids) == 0L) {
    return(ids)
  }

  known_ok <- tax_prob |>
    dplyr::filter(
      .data[["rank"]] == .env$rank,
      !is.na(.data$taxon),
      .data$taxon != "unk",
      .data$prob >= .env$min_prob
    ) |>
    dplyr::distinct(.data[[id_col]]) |>
    dplyr::pull(1)

  need <- setdiff(ids, known_ok)

  if (!is.null(screen) && length(scope_ranks) > 0L) {
    checkmate::assert_data_frame(screen)
    checkmate::assert_names(
      names(screen),
      must.include = c(id_col, scope_ranks)
    )
    exceptions <- best_hit_ingroup_exceptions()
    screen_need <- rep(FALSE, nrow(screen))
    for (i in seq_along(scope_ranks)) {
      col <- scope_ranks[[i]]
      taxon_i <- scope_taxa[[i]]
      vals <- screen[[col]]
      screen_need <- screen_need |
        (!is.na(vals) & !vals %in% c(taxon_i, exceptions))
    }
    screened_ids <- screen[[id_col]][screen_need]
    need <- union(need, screened_ids)
  }

  # Preserve type/order relative to all_ids
  ids[ids %in% need]
}

#' Combine primary and fallback taxonomy tables for a cascade
#'
#' Drops primary rows for `reclassify_ids`, binds fallback rows, and tags a
#' `classifier` column.
#'
#' @param primary (`data.frame`) primary classifier output
#' @param fallback (`data.frame`) fallback classifier output
#' @param reclassify_ids IDs that were sent to the fallback
#' @param id_col (`character`) ID column name
#' @param primary_name (`character`) label for primary rows
#' @param fallback_name (`character`) label for fallback rows
#' @return combined `data.frame`
#' @export
combine_tax_prob <- function(
  primary,
  fallback,
  reclassify_ids,
  id_col = "seq_idx",
  primary_name = "primary",
  fallback_name = "fallback"
) {
  checkmate::assert_data_frame(primary)
  checkmate::assert_data_frame(fallback)
  checkmate::assert_string(id_col)
  checkmate::assert_string(primary_name)
  checkmate::assert_string(fallback_name)
  checkmate::assert_names(names(primary), must.include = id_col)
  checkmate::assert_names(names(fallback), must.include = id_col)

  kept_primary <- primary[
    !primary[[id_col]] %in% reclassify_ids,
    ,
    drop = FALSE
  ]
  if (nrow(kept_primary) > 0L) {
    kept_primary$classifier <- primary_name
  } else {
    kept_primary$classifier <- character()
  }
  if (nrow(fallback) > 0L) {
    fallback$classifier <- fallback_name
  } else if (!"classifier" %in% names(fallback)) {
    fallback$classifier <- character()
  }

  # Align columns for bind_rows
  all_cols <- union(names(kept_primary), names(fallback))
  for (col in setdiff(all_cols, names(kept_primary))) {
    kept_primary[[col]] <- rep(NA, nrow(kept_primary))
  }
  for (col in setdiff(all_cols, names(fallback))) {
    fallback[[col]] <- rep(NA, nrow(fallback))
  }
  dplyr::bind_rows(kept_primary[all_cols], fallback[all_cols])
}

#' Build prior taxonomy rows for a classifier's closed scope
#'
#' @param ids vector of sequence IDs
#' @param id_col (`character`) name for the ID column
#' @param scope_ranks (`character`) closed-scope ranks
#' @param scope_taxa (`character`) taxa at those ranks
#' @param ranks (`character`) full ordered rank list for factor conversion
#' @return `data.frame` with id, rank, taxon, prob=1, parent_taxonomy
#' @export
scope_prior_taxonomy <- function(
  ids,
  id_col = "seq_idx",
  scope_ranks,
  scope_taxa,
  ranks = tax_ranks()
) {
  checkmate::assert_string(id_col)
  checkmate::assert_character(scope_ranks, any.missing = FALSE)
  checkmate::assert_character(scope_taxa, any.missing = FALSE)
  if (length(scope_ranks) != length(scope_taxa)) {
    stop("`scope_ranks` and `scope_taxa` must have the same length.")
  }
  if (length(scope_ranks) == 0L || length(ids) == 0L) {
    out <- tibble::tibble(
      rank = rank2factor(character(), ranks),
      taxon = character(),
      prob = numeric(),
      parent_taxonomy = character()
    )
    out[[id_col]] <- if (is.integer(ids)) integer() else character()
    return(dplyr::relocate(out, dplyr::all_of(id_col)))
  }
  priors <- tibble::tibble(
    rank = rank2factor(scope_ranks, ranks),
    taxon = scope_taxa,
    prob = 1.0
  ) |>
    dplyr::arrange(dplyr::desc(rank)) |>
    dplyr::mutate(
      parent_taxonomy = purrr::accumulate(
        .x = dplyr::lag(taxon, default = NA_character_),
        .f = \(x, y) if (is.na(x)) y else paste(x, y, sep = ",")
      )
    )
  tidyr::crossing(tibble::tibble(!!id_col := ids), priors)
}

#' Build study-ingroup / outgroup / unknown labels per sequence
#'
#' Assigned taxonomy from open classifiers wins at `known_ranks`. Closed
#' classifier winners (and sequences with no assignment at those ranks) fall
#' back to `best_hit` when provided. Injected scope priors are not treated as
#' assigned evidence.
#'
#' @param tax_prob long taxonomy table (`seq_id`/`seq_idx`, `rank`, `taxon`,
#'   `prob`, optional `classifier`)
#' @param best_hit wide best-hit table with `id_col` and `known_ranks` columns,
#'   or `NULL`
#' @param id_col ID column name
#' @param known_ranks study-ingroup ranks
#' @param known_taxa study-ingroup taxa
#' @param open_classifiers classifier labels whose assignments are authoritative
#'   at `known_ranks`; empty means no classifier assignments count (use
#'   `best_hit` only)
#' @param has_classifier_column whether `tax_prob` has a `classifier` column
#' @return `list` of tibbles: `known_ingroup`, `known_outgroup`, `unknown_outin`
#' @export
study_ingroup_membership <- function(
  tax_prob,
  best_hit = NULL,
  id_col = "seq_id",
  study_ranks = known_ranks(),
  study_taxa = known_taxa(),
  open_classifiers = character(),
  has_classifier_column = "classifier" %in% names(tax_prob)
) {
  checkmate::assert_data_frame(tax_prob)
  checkmate::assert_string(id_col)
  checkmate::assert_character(study_ranks, any.missing = FALSE)
  checkmate::assert_character(study_taxa, any.missing = FALSE)
  checkmate::assert_character(open_classifiers, any.missing = FALSE)
  exceptions <- best_hit_ingroup_exceptions()
  known_ranks <- study_ranks
  known_taxa <- study_taxa

  all_ids <- unique(tax_prob[[id_col]])
  if (!is.null(best_hit)) {
    all_ids <- union(all_ids, unique(best_hit[[id_col]]))
  }

  # Assigned rows at known ranks from open classifiers only
  assigned <- tax_prob
  if (length(known_ranks) == 0L) {
    assigned <- assigned[0, , drop = FALSE]
  } else {
    assigned <- dplyr::filter(
      assigned,
      .data$rank %in% .env$known_ranks,
      !is.na(.data$taxon),
      .data$taxon != "unk"
    )
    if (has_classifier_column) {
      if (length(open_classifiers) == 0L) {
        assigned <- assigned[0, , drop = FALSE]
      } else {
        assigned <- dplyr::filter(
          assigned,
          .data$classifier %in% .env$open_classifiers
        )
      }
    } else if (length(open_classifiers) == 0L) {
      # No classifier column and no open classifiers: ignore tax_prob for
      # study-ingroup ranks (closed sole classifier → best_hit only).
      assigned <- assigned[0, , drop = FALSE]
    }
  }

  assigned_wide <- if (nrow(assigned) == 0L) {
    out <- tibble::tibble(!!id_col := all_ids[0])
    for (r in known_ranks) {
      out[[r]] <- character()
    }
    out
  } else {
    assigned |>
      dplyr::arrange(dplyr::desc(prob)) |>
      dplyr::distinct(.data[[id_col]], rank, .keep_all = TRUE) |>
      dplyr::select(dplyr::all_of(c(id_col, "rank", "taxon"))) |>
      tidyr::pivot_wider(names_from = "rank", values_from = "taxon")
  }

  # Fill gaps from best_hit
  if (!is.null(best_hit) && length(known_ranks) > 0L) {
    checkmate::assert_names(
      names(best_hit),
      must.include = c(id_col, known_ranks)
    )
    merged <- dplyr::full_join(
      assigned_wide,
      dplyr::select(best_hit, dplyr::all_of(c(id_col, known_ranks))),
      by = id_col,
      suffix = c("", "_bh")
    )
    for (r in known_ranks) {
      bh <- paste0(r, "_bh")
      if (bh %in% names(merged)) {
        merged[[r]] <- dplyr::coalesce(merged[[r]], merged[[bh]])
        merged[[bh]] <- NULL
      }
    }
  } else {
    merged <- assigned_wide
  }

  # Ensure all ids present
  if (length(all_ids) > 0L) {
    merged <- dplyr::full_join(
      tibble::tibble(!!id_col := all_ids),
      merged,
      by = id_col
    )
  }

  if (length(known_ranks) == 0L || nrow(merged) == 0L) {
    empty <- merged[0, , drop = FALSE]
    return(list(
      known_ingroup = empty,
      known_outgroup = empty,
      unknown_outin = empty
    ))
  }

  outgroup_cols <- character(length(known_ranks))
  ingroup_ok <- rep(TRUE, nrow(merged))
  outgroup_any <- rep(FALSE, nrow(merged))
  for (i in seq_along(known_ranks)) {
    rank_i <- known_ranks[[i]]
    taxon_i <- known_taxa[[i]]
    vals <- merged[[rank_i]]
    is_out <- !is.na(vals) & !vals %in% c(taxon_i, exceptions)
    is_in <- !is.na(vals) & vals == taxon_i
    outgroup_cols[i] <- paste0(rank_i, "_outgroup")
    merged[[outgroup_cols[i]]] <- is_out
    outgroup_any <- outgroup_any | is_out
    ingroup_ok <- ingroup_ok & (is_in | is.na(vals))
  }
  # known_ingroup: all known ranks that are present match study taxa, and the
  # tip known rank is present and matches
  tip <- known_ranks[[length(known_ranks)]]
  tip_taxon <- known_taxa[[length(known_taxa)]]
  known_ingroup <- merged[
    !outgroup_any & !is.na(merged[[tip]]) & merged[[tip]] == tip_taxon,
    ,
    drop = FALSE
  ]
  known_outgroup <- merged[outgroup_any, , drop = FALSE]
  unknown_outin <- merged[
    !outgroup_any &
      (is.na(merged[[tip]]) | merged[[tip]] != tip_taxon),
    ,
    drop = FALSE
  ]

  keep_cols <- c(id_col, known_ranks)
  list(
    known_ingroup = dplyr::select(known_ingroup, dplyr::all_of(keep_cols)),
    known_outgroup = dplyr::select(known_outgroup, dplyr::all_of(keep_cols)),
    unknown_outin = dplyr::select(unknown_outin, dplyr::all_of(keep_cols))
  )
}
