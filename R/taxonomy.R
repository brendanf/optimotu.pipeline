#### taxonomic ranks ####

#' Get or set the taxonomic ranks to use in the pipeline
#' @return a `character` vector of taxonomic ranks
#' @export
tax_ranks <- function() {
  getOption("optimotu.pipeline.tax_ranks",
            c("kingdom", "phylum", "class", "order", "family", "genus", "species"))
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
  options("optimotu.pipeline.tax_ranks" = value)[["optimotu.pipeline.tax_ranks"]]
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
  options("optimotu.pipeline.rank_offset" = value)[["optimotu.pipeline.rank_offset"]]
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
  options("optimotu.pipeline.known_ranks" = value)[["optimotu.pipeline.known_ranks"]]
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
  options("optimotu.pipeline.known_taxa" = value)[["optimotu.pipeline.known_taxa"]]
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
tip_rank <- function() {
  tax_ranks()[length(tax_ranks())]
}

#' @describeIn tax_opts Get the tip taxonomic rank as a symbol
#' @export
tip_rank_var <- function() {
  rlang::sym(tip_rank())
}

#' Define the taxonomic ranks to use in the pipeline
#'
#' @param ranks (`list` or `character` vector) the taxonomic ranks to use in the
#' pipeline, in order from most inclusive (e.g., kingdom) to least inclusive
#' (e.g., species). Values may be either named or unnamed. When named, the name is
#' taken to be the rank, and the value is the "in-group" taxon at that rank, i.e.
#' the taxon for which results are desired. When unnamed, the value is taken to
#' be the rank. Example: `list(kingdom = "Fungi", "phylum", "class", "order", "family", "genus", "species")`
#' @return `NULL`.  This function is called for its side effect, which is to
#' configure global options.
#' @export
define_taxonomy <- function(ranks) {
  checkmate::assert(
    checkmate::check_list(
      ranks,
      types = c("character", "list"),
      min.len = 1
    ),
    checkmate::check_character(
      ranks,
      unique = TRUE,
      min.len = 1
    )
  )
  KNOWN_TAXA <- purrr::keep(ranks, ~ dplyr::cumall(checkmate::test_list(.x))) |>
    unlist()
  UNKNOWN_RANKS <- purrr::discard(ranks, ~ dplyr::cumall(checkmate::test_list(.x))) |>
    unlist()
  if (length(UNKNOWN_RANKS) == 0 || !is.null(names(UNKNOWN_RANKS))) {
    stop(
      "Option 'ranks' should start from the most inclusive rank (e.g. kingdom)\n",
      "  and continue to the least inclusive rank (e.g. species).  Optionally the first\n",
      "  rank(s) may be defined (e.g. '- kingdom: Fungi') but subsequent ranks must be \n",
      "  undefined (e.g. '- class')."
    )
  }
  set_known_ranks(names(KNOWN_TAXA))
  set_known_taxa(unname(KNOWN_TAXA))
  set_tax_ranks(c(known_ranks(), UNKNOWN_RANKS))
}

#' Convert taxonomic ranks from character vectors and integers to ordered factors
#'
#' The least inclusive rank is given the smallest value in the ordering,
#' so for instance `"species" < "kingdom"`.
#'
#' @param x the taxonomic ranks to convert
#' @return an ordered factor of taxonomic ranks
#' @export
rank2factor <- function(x) {
  factor(x, levels = rev(tax_ranks()), ordered = TRUE)
}

#' @rdname rank2factor
#' @export
int2rankfactor <- function(x) {
  rank2factor(tax_ranks()[x + rank_offset()])
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
  ranks[rank2factor(ranks) > x]
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
  ranks[rank2factor(ranks) < x]
}

#' @rdname superranks
#' @return a `list` of `symbol`s representing the subordinate ranks
#' @export
subrank_vars <- function(x = ingroup_rank(), ranks = tax_ranks()) {
  rlang::syms(subranks(x, ranks))
}

#' Combine tip classifications to build a full PROTAX taxonomy
#' @param ... (`character` vectors) the classifications to combine; these should
#' be comma-delimited classifications from most inclusive rank to least
#' inclusive.
#' @return a `data.frame` with columns `taxon_id`, `parent_id`, `rank`,
#' `classification`, and `prior`, as required to write a Protax "`taxonomy`" file
#' @export
build_taxonomy <- function(...) {
  # avoid R CMD check NOTE for undeclared global variables
  classification <- parent <- prior <- taxon_id <- parent_id <- NULL

  tax <- tibble::tibble(
    classification = union(...),
    rank = ifelse(
      classification == "root",
      0L,
      stringr::str_count(classification, stringr::fixed(",")) + 1L),
    parent = ifelse(
      rank <= 1,
      "root",
      sub(",[^,]+$", "", classification)
    )
  )
  tax <- split(tax, tax$rank)
  tax[[1]]$taxon_id = 0L
  tax[[1]]$prior = 1
  tax[[1]]$parent_id = 0L
  tax[[2]]$taxon_id = 1L:2L
  tax[[2]]$prior = c(0.99, 0.01)

  for (r in length(tax):3L) {
    if (r == length(tax)) {
      tax[[r]]$prior = 0.99/nrow(tax[[r]])
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
    dplyr::mutate(!!root_rank() := ifelse(!!root_rank() == "root", NA_character_, !!root_rank()))

  # remove all taxa with children
  for (rank in 6:1) {
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
    tax_i <- tax_i[!is.na(tax_i[[tax_ranks()[i]]]),]
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
        classification = paste(parent_classification, {{rankname}}, sep = ",")
      )
    } else {
      tax_i <- dplyr::mutate(
        tax_i,
        classification = {{rankname}}
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
  regex <- paste0("(^([^,]+,){", rank-1L, "}[^,]+).*")
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
