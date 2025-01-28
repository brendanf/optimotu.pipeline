#' Calculate clustering thresholds for each taxon, falling back to its ancestor
#' taxa as necessary
#'
#' @param rank (`character` string) the rank within which clustering will be
#' performed
#' @param conf_level (`character` string) as in `fmeasure_optima`
#' @param taxon_table (`data.frame`) taxonomy table; column "seq_id" gives the
#' sequence ID, and columns {`root_rank()`} to {`tip_rank()`} (e.g., "kingdom"
#' to "species") give the taxonomy at each rank.
#' @param fmeasure_optima (`data.frame`) optimum clustering thresholds within
#' various taxa; column "rank" gives the rank which is approximated by
#' clustering; "superrank" gives the rank of the taxon within which the
#' clustering threshold was optimized; "supertaxon" gives that taxon name;
#' "conf_level" gives a string description of the confidence level threshold for
#' taxonomic assignments; "threshold" gives the optimum clustering threshold;
#' "f_measure" gives the F measure at the optimum threshold.
#' @param default (`character string`) default taxon to define threshold to use
#' when taxonomy is unknown. default: `ingroup_taxon()`
#'
#'
#' @return (named `numeric`, where names are taxa and values are clustering
#' thresholds)
#' unknown sequences, grouped by taxonomy at the parent rank. Sequences where
#' the parent rank is also unknown are in the item named `"_NA_"`.
#' @export
calc_taxon_thresholds <- function(rank, conf_level, taxon_table,
                                  fmeasure_optima, default = ingroup_taxon()) {
  # avoid R CMD check NOTE: no visible binding for global variable
  superrank <- supertaxon <- threshold <- NULL

  rank_name <- rlang::sym(rank)
  dplyr::select(taxon_table, !!root_rank():!!rank_name) |>
    dplyr::filter(!is.na(!!rank_name)) |>
    unique() |>
    purrr::reduce(
      c(superranks(rank), rank),
      function(thresholds, r) {
        dplyr::left_join(
          thresholds,
          dplyr::filter(
            fmeasure_optima,
            rank == subranks(!!rank)[1],
            superrank == r,
            conf_level == !!conf_level
          ) |>
            dplyr::select(
              !!r := supertaxon,
              !!paste0("threshold_", r) := threshold
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
      fmeasure_optima,
      rank == subranks(!!rank)[1],
      supertaxon == default,
      conf_level == !!conf_level
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

#' Calculate clustering thresholds for each taxon, falling back to its ancestor
#' taxa as necessary
#'
#' @param rank (`character` string) the rank within which clustering will be
#' performed
#' @param conf_level (`character` string) as in `fmeasure_optima`
#' @param taxon_table (`data.frame`) taxonomy table; column "seq_id" gives the
#' sequence ID, and columns {`root_rank()`} to {`tip_rank()`} (e.g., "kingdom" to
#' "species") give the taxonomy at each rank.
#' @param fmeasure_optima (`data.frame`) optimum clustering thresholds within
#' various taxa; column "rank" gives the rank which is approximated by
#' clustering; "superrank" gives the rank of the taxon within which the
#' clustering threshold was optimized; "supertaxon" gives that taxon name;
#' "conf_level" gives a string description of the confidence level threshold for
#' taxonomic assignments; "threshold" gives the optimum clustering threshold;
#' "f_measure" gives the F measure at the optimum threshold.
#' @param default (`character string`) default taxon to define threshold to use
#' when taxonomy is unknown. default: `ingroup_taxon()`
#'
#'
#' @return (named `list` of `double` vectors)
calc_subtaxon_thresholds <- function(rank, conf_level, taxon_table,
                                     fmeasure_optima, default = ingroup_taxon()) {
  # avoid R CMD check NOTE: no visible binding for global variable
  subrank <- superrank <- supertaxon <- threshold <- NULL

  rank_name <- rlang::sym(rank)
  dplyr::select(taxon_table, {{root_rank()}}:{{rank_name}}) |>
    tidyr::crossing(subrank = subranks(rank)) |>
    dplyr::filter(!is.na(!!rank_name)) |>
    unique() |>
    purrr::reduce(
      c(superranks(rank), rank),
      function(thresholds, r) {
        dplyr::left_join(
          thresholds,
          dplyr::filter(
            fmeasure_optima,
            rank %in% subranks(!!rank),
            superrank == r,
            conf_level == !!conf_level
          ) |>
            dplyr::select(
              subrank = rank,
              !!r := supertaxon,
              !!paste0("threshold_", r) := threshold
            ),
          by = c("subrank", r)
        )
      },
      .init = _
    ) |>
    (\(x) dplyr::transmute(
      x,
      subrank = rank2factor(subrank),
      !!rank_name := !!rank_name,
      threshold = x |>
        dplyr::select(dplyr::starts_with("threshold")) |>
        rev() |>
        do.call(dplyr::coalesce, args = _) |>
        threshold_as_dist()
    )
    )() |>
    dplyr::arrange(subrank) |>
    (\(x) split(x, x[[rank]]))() |>
    lapply(dplyr::select, !any_of(rank)) |>
    lapply(tibble::deframe) |>
    lapply(cummax) |>
    c(
      "_NA_" = dplyr::filter(
        fmeasure_optima,
        rank %in% subranks(!!rank),
        supertaxon == default,
        conf_level == !!conf_level
      ) |>
        dplyr::transmute(rank = rank2factor(rank), threshold_as_dist(threshold)) |>
        dplyr::arrange(rank) |>
        tibble::deframe() |>
        cummax() |>
        list()
    )
}
