#' Calculate clustering thresholds for each taxon, falling back to its ancestor
#' taxa as necessary
#'
#' @param rank (`character` string) the rank within which clustering will be
#' performed
#' @param taxon_table (`data.frame`) "wide" taxonomy table; column "seq_id"
#' gives the sequence ID, and columns {`root_rank()`} to {`tip_rank()`} (e.g.,
#' "kingdom" to "species") give the taxonomy at each rank.
#' @param optima (`data.frame`) optimum clustering thresholds within
#' various taxa;
#' column "rank" gives the rank which is approximated by
#' clustering;
#' "superrank" gives the rank of the taxon within which the clustering threshold
#' was optimized;
#' "supertaxon" gives that taxon name;
#' "threshold" gives the optimum clustering threshold;
#' "value" gives the value of the optimization target at the optimum threshold;
#' optionally, "conf_level" gives a string description of the confidence level;
#' optionally, "metric" gives a string description of the optimization target.
#' Typically, "conf_level" and "measure" are only needed when optimization has
#' been performed for multiple confidence levels or measures, and so the
#' arguments by the same name should also be included.
#' @param conf_level (`character` string or NULL) if given, then only rows of
#' `optima` with this value in the "conf_level" column will be used
#' @param metric (`character` string or NULL) if given, then only rows of
#' `optima` with this value in the "metric" column will be used
#' @param default (`character string`) default taxon to define threshold to use
#' when taxonomy is unknown. default: `ingroup_taxon()`
#'
#'
#' @return (named `numeric`, where names are taxa and values are clustering
#' thresholds)
#' unknown sequences, grouped by taxonomy at the parent rank. Sequences where
#' the parent rank is also unknown are in the item named `"_NA_"`.
#' @export
calc_taxon_thresholds <- function(rank, taxon_table, optima,
                                  conf_level = NULL, metric = NULL, default = ingroup_taxon()) {
  # avoid R CMD check NOTE: no visible binding for global variable
  superrank <- supertaxon <- threshold <- NULL

  # check arguments
  checkmate::assert_character(rank)
  checkmate::assert_subset(rank, tax_ranks())
  checkmate::assert_true(rank2factor(rank) > rank2factor(tip_rank()))
  checkmate::assert_data_frame(taxon_table)
  checkmate::assert_names(
    names(taxon_table),
    must.include = c("seq_id", superranks(rank), rank, subranks(rank)[1])
  )
  checkmate::assert_data_frame(optima)
  checkmate::assert_names(
    names(optima),
    must.include = c("rank", "superrank", "supertaxon", "threshold")
  )
  checkmate::assert_string(conf_level, null.ok = TRUE)
  checkmate::assert_string(metric, null.ok = TRUE)
  checkmate::assert_string(default)

  # symbol version of the rank, for NSE
  rank_name <- rlang::sym(rank)

dplyr::select(taxon_table, !!root_rank_var():!!rank_name) |>
    dplyr::filter(!is.na(!!rank_name)) |>
    unique() |>
    purrr::reduce(
      c(superranks(rank), rank),
      function(thresholds, r) {
        dplyr::left_join(
          thresholds,
          dplyr::filter(
            optima,
            rank == subranks(!!rank)[1],
            superrank == r,
            if (is.null(!!conf_level)) TRUE else conf_level == !!conf_level,
            if (is.null(!!metric)) TRUE else metric == !!metric
          ) |>
            dplyr::select(
              !!r := supertaxon,
              !!(paste0("threshold_", r)) := threshold
            ),
          by = r
        )
      },
      .init = _
    ) |>
    (\(x) dplyr::transmute(
      x,
      !!rank_name := !!rank_name,
      threshold = x |>
        dplyr::select(dplyr::starts_with("threshold")) |>
        rev() |>
        do.call(dplyr::coalesce, args = _)
    )
    )() |>
    tibble::deframe() |>
    c("_NA_" = dplyr::filter(
      optima,
      rank == subranks(!!rank)[1],
      supertaxon == default,
      if (is.null(!!conf_level)) TRUE else conf_level == !!conf_level,
      if (is.null(!!metric)) TRUE else metric == !!metric
    )$threshold)
}

#' Convert thresholds which may be distances or similarities to distances
#' @param thresholds (`numeric`) thresholds to convert. These may be distances
#' (0 is identity) with values between 0.0 and 0.5, or between 0.0 and 50.0,
#' in which case they are interpreted as percentages;
#' or similarities (1 or 100 is identity) with values between 0.5 and 1.0, or
#' between 50.0 and 100.0. Note that if *all* percentage distances are < 0.5,
#' they will be misinterpreted as fractional distances.
#' @return (`numeric`) thresholds as distances, in the range from 0.0 to 1.0
#' @export
threshold_as_dist <- function(thresholds) {
  if (any(thresholds > 50)) {
    1 - 0.01 * thresholds
  } else if (any(thresholds > 1)) {
    0.01 * thresholds
  } else if (any(thresholds > 0.5) ) {
    1 - thresholds
  } else {
    thresholds
  }
}
