#' Calculate a full preclosed taxon table
#'
#' This calculates the full table of inputs to be passed to
#' [do_closed_ref_cluster()], with a single group for each parent taxon.
#' It is also used internally by [large_preclosed_taxon_table()] and
#' [small_preclosed_taxon_table()].
#' @param known_taxon_table (`data.frame`) a table of known taxon assignments
#' @param asv_taxsort (`data.frame`) a table of ASV taxonomic assignments.
#'   Must include `seq_idx` (integer, 1-based index after taxonomic sorting)
#'   and exactly one of `seq_idx_in` (1-based index before sorting) or
#'   `seq_id` (character sequence identifier). When using `seq_id`, the taxon
#'   table must include `seq_id`; a taxon table with only `seq_idx` is not
#'   allowed.
#' @param rank (`character`) the taxonomic rank to calculate
#' @param parent_rank (`character`) the taxonomic rank to use as a parent
#' @param tax_ranks (`character`) the taxonomic ranks to use
#'
#' @return (`data.frame`) a table of sequences for which closed-reference
#' clustering is required. This includes both sequences which are not assigned
#' to the target rank which will act as queries, and identified sequences
#' which belong to the same parent rank as one or more queries. The results
#' includes a `tar_group` column for dynamic branching over each parent
#' taxon. If there are no queries, a single row is returned to branch over.
#' @export
full_preclosed_taxon_table <- function(
  known_taxon_table,
  asv_taxsort,
  rank,
  parent_rank,
  tax_ranks
) {
  # avoid R CMD check NOTEs
  seq_idx <- seq_id <- seq_idx_in <- NULL
  # check input
  checkmate::assert_string(rank)
  checkmate::assert_string(parent_rank)
  checkmate::assert_character(tax_ranks)
  checkmate::assert_subset(rank, tax_ranks)
  checkmate::assert_subset(parent_rank, tax_ranks)

  checkmate::assert_data_frame(known_taxon_table)
  checkmate::assert_names(
    names(known_taxon_table),
    must.include = c(rank, parent_rank)
  )
  checkmate::assert(
    checkmate::check_names(names(known_taxon_table), must.include = "seq_id"),
    checkmate::check_names(names(known_taxon_table), must.include = "seq_idx"),
    combine = "or"
  )
  checkmate::assert_character(known_taxon_table[[rank]])
  checkmate::assert_character(known_taxon_table[[parent_rank]])
  checkmate::assert_integerish(known_taxon_table[["seq_idx"]], null.ok = TRUE)
  checkmate::assert_character(known_taxon_table[["seq_id"]], null.ok = TRUE)

  rank_sym <- rlang::sym(rank)
  super_ranks <- superranks(rank, tax_ranks)

  out <- join_table_with_taxsort(known_taxon_table, asv_taxsort) |>
    dplyr::arrange(dplyr::pick(all_of(c(super_ranks, "seq_idx")))) |>
    dplyr::group_by(dplyr::pick(all_of(parent_rank))) |>
    dplyr::filter(
      any(is.na(!!rank_sym)),
      any(!is.na(!!rank_sym))
    ) |>
    targets::tar_group()
  if ("seq_id" %in% names(out)) {
    out <- dplyr::select(out, -any_of("seq_idx_in"))
  }
  ensure_table(out)
}

#' Calculate a pre-closed-reference table for taxa with many sequences
#'
#' This calculates a pre-closed-reference table for taxa with many sequences,
#' which may be efficiently parallelized within a single parent taxon.
#'
#' @inheritParams full_preclosed_taxon_table
#' @param min_ops (`integer` scalar) the minimum number of operations required
#' for a single parent taxon to be considered "large"
#'
#' @return (`data.frame`) a table of sequences for which closed-reference
#' clustering is required. This includes both sequences which are not assigned
#' to the target rank which will act as queries, and identified sequences
#' which belong to the same parent rank as one or more queries. The results
#' includes a `tar_group` column for dynamic branching over each parent
#' taxon. If there are no queries, a single row is returned to branch over.
#' @export
large_preclosed_taxon_table <- function(
  known_taxon_table,
  asv_taxsort,
  rank,
  parent_rank,
  tax_ranks,
  min_ops = 1e6
) {
  checkmate::assert_number(min_ops, lower = 1)

  rank_sym <- rlang::sym(rank)

  fptt <- full_preclosed_taxon_table(
    known_taxon_table,
    asv_taxsort,
    rank,
    parent_rank,
    tax_ranks
  )

  out <- fptt |>
    dplyr::group_by(dplyr::pick(all_of(parent_rank))) |>
    dplyr::filter(sum(is.na(!!rank_sym)) * sum(!is.na(!!rank_sym)) > min_ops) |>
    targets::tar_group()
  ensure_table(out)
}

#' Calculate a pre-closed-reference table for taxa with few sequences
#'
#' This calculates a pre-closed-reference table for taxa with few sequences,
#' and groups them into approximately equally sized execution groups.
#'
#' @inheritParams full_preclosed_taxon_table
#' @param max_ops (`integer` scalar) the maximum number of operations to put
#' together in a single execution group
#'
#' @return (`data.frame`) a table of sequences for which closed-reference
#' clustering is required. This includes both sequences which are not assigned
#' to the target rank which will act as queries, and identified sequences
#' which belong to the same parent rank as one or more queries. The results
#' includes a `tar_group` column for dynamic branching over execution groups.
#' If there are no queries, a single row is returned to branch over.
#' @export
small_preclosed_taxon_table <- function(
  known_taxon_table,
  asv_taxsort,
  rank,
  parent_rank,
  tax_ranks,
  max_ops = 1e6
) {
  checkmate::assert_number(max_ops, lower = 1)

  rank_sym <- rlang::sym(rank)

  fptt <- full_preclosed_taxon_table(
    known_taxon_table,
    asv_taxsort,
    rank,
    parent_rank,
    tax_ranks
  )

  out <- fptt |>
    dplyr::group_by(dplyr::pick(all_of(parent_rank))) |>
    dplyr::mutate(ops = sum(is.na(!!rank_sym)) * sum(!is.na(!!rank_sym))) |>
    dplyr::filter(ops <= max_ops) |>
    dplyr::select(-all_of("tar_group")) |>
    # Ungroup so distribute_tasks() sees all taxa, not one group at a time.
    dplyr::ungroup()

  ops <- dplyr::select(out, all_of(c(parent_rank, "ops"))) |>
    dplyr::distinct() |>
    dplyr::mutate(tar_group = distribute_tasks(ops, max_ops)) |>
    dplyr::select(-ops)

  dplyr::left_join(out, ops, by = parent_rank) |>
    dplyr::select(-ops) |>
    ensure_table()
}

#' Run closed-reference clustering
#'
#' This function calls [optimotu::closed_ref_cluster()] for one or more
#' parent taxa.
#'
#' @param preclosed_taxon_table (`data.frame`) a table of sequences for which
#' closed-reference clustering is required. Must have columns `seq_idx`
#' (`integer`) giving the indices of sequences in `seq_file`, as well as
#' columns with names matching `rank` and `parent_rank`.
#' @param rank (`character`) the taxonomic rank to calculate clusters for
#' @param parent_rank (`character`) the taxonomic rank to use as a parent
#' @param seq_file (`character`) the path to the sequence file
#' @param seq_file_index (`character` or `fastqindexr_index`) file path or
#'  [`fastqindexr_index`][fastqindexr::create_index()] object for the sequence
#'  file index
#' @param thresholds (`numeric`) a named vector of thresholds to use for
#' closed-reference clustering. These are interpreted as distances on the
#' interval \[0, 100\], where 0 is a perfect match and 100 is a perfect
#' mismatch.
#' @param dist_config (`list`) a list of configuration options for the
#' distance calculation, as returned by [optimotu::dist_config()].
#' @param parallel_config (`list`) a list of configuration options for the
#' parallelization, as returned by [optimotu::parallel_config()].
#'
#' @return (`data.frame`) a two-column table of mapping between query
#' sequences (those where the `rank` column in `preclosed_taxon_table` is `NA`)
#' and reference sequences (those where the `rank` column is not `NA`).
#' @export
do_closed_ref_cluster <- function(
  preclosed_taxon_table,
  rank,
  parent_rank,
  seq_file,
  seq_file_index,
  thresholds,
  dist_config,
  parallel_config
) {
  checkmate::assert_data_frame(preclosed_taxon_table)
  checkmate::assert_string(rank)
  checkmate::assert_string(parent_rank)
  checkmate::assert_names(
    names(preclosed_taxon_table),
    must.include = c("seq_idx", rank, parent_rank)
  )
  checkmate::assert_integerish(preclosed_taxon_table$seq_idx)
  checkmate::assert_character(preclosed_taxon_table[[rank]])
  checkmate::assert_character(preclosed_taxon_table[[parent_rank]])
  checkmate::assert_file_exists(seq_file, "r")
  checkmate::assert(
    checkmate::check_file_exists(seq_file_index, "r"),
    checkmate::check_class(seq_file_index, "fastqindexr_index"),
    combine = "or"
  )
  checkmate::assert_named(thresholds)
  checkmate::assert_numeric(thresholds)
  checkmate::assert_class(dist_config, "optimotu_dist_config")
  checkmate::assert_class(parallel_config, "optimotu_parallel_config")

  if (dplyr::n_distinct(preclosed_taxon_table[[parent_rank]]) > 1L) {
    out <- preclosed_taxon_table |>
      dplyr::group_by(dplyr::pick(all_of(parent_rank))) |>
      dplyr::group_modify(
        function(x, y) {
          do_closed_ref_cluster(
            preclosed_taxon_table = x,
            rank = rank,
            parent_rank = parent_rank,
            seq_file = seq_file,
            seq_file_index = seq_file_index,
            thresholds = thresholds,
            dist_config = dist_config,
            parallel_config = parallel_config
          )
        },
        .keep = TRUE
      ) |>
      dplyr::ungroup() |>
      dplyr::select(-all_of(parent_rank))
    return(out)
  }

  unknowns <- is.na(preclosed_taxon_table[[rank]])
  taxon <- preclosed_taxon_table[[parent_rank]][1]

  if (all(unknowns) || !any(unknowns)) {
    return(
      tibble::tibble(
        seq_id = character(),
        ref_id = character()
      )
    )
  }
  query_seq <- fastqindexr::extract_sequences(
    index = seq_file_index,
    seq_idx = preclosed_taxon_table$seq_idx[unknowns],
    file = seq_file,
    return = "seq"
  )
  ref_seq <- fastqindexr::extract_sequences(
    index = seq_file_index,
    seq_idx = preclosed_taxon_table$seq_idx[!unknowns],
    file = seq_file,
    return = "seq"
  )
  closed_ref_out <- optimotu::closed_ref_cluster(
    query = query_seq,
    ref = ref_seq,
    threshold = optimotu::threshold_as_dist(thresholds[taxon]),
    dist_config = dist_config,
    parallel_config = parallel_config
  )
  closed_ref_out
}

#' Calculate a full pre-denovo taxon table
#'
#' This calculates the full table of inputs to be passed to
#' [do_denovo_cluster()], with a single group for each parent taxon.
#' It is also used internally by [large_predenovo_taxon_table()] and
#' [small_predenovo_taxon_table()].
#' @inheritParams large_predenovo_taxon_table
#' @return (`data.frame`) a table of sequences for which denovo clustering
#' should be performed.
#' @keywords internal
full_predenovo_taxon_table <- function(
  closedref_taxon_table,
  asv_taxsort,
  rank,
  parent_rank,
  tax_ranks
) {
  # avoid R CMD check NOTEs
  seq_idx <- seq_id <- seq_idx_in <- NULL
  checkmate::assert_data_frame(closedref_taxon_table)
  checkmate::assert_data_frame(asv_taxsort)
  checkmate::assert_string(rank)
  checkmate::assert_string(parent_rank)
  checkmate::assert_character(tax_ranks)
  checkmate::assert_subset(rank, tax_ranks)
  checkmate::assert_subset(parent_rank, tax_ranks)
  checkmate::assert_names(
    names(closedref_taxon_table),
    must.include = c(rank, parent_rank)
  )
  checkmate::assert_character(closedref_taxon_table[[rank]])
  checkmate::assert_character(closedref_taxon_table[[parent_rank]])
  checkmate::assert(
    checkmate::check_names(
      names(closedref_taxon_table),
      must.include = "seq_id"
    ),
    checkmate::check_names(
      names(closedref_taxon_table),
      must.include = "seq_idx"
    ),
    combine = "or"
  )
  checkmate::assert_integerish(
    closedref_taxon_table[["seq_idx"]],
    null.ok = TRUE
  )
  checkmate::assert_character(closedref_taxon_table[["seq_id"]], null.ok = TRUE)

  rank_sym <- rlang::sym(rank)
  super_ranks <- superranks(rank, tax_ranks)

  out <- closedref_taxon_table |>
    dplyr::filter(is.na(!!rank_sym))

  out <- join_table_with_taxsort(out, asv_taxsort) |>
    dplyr::arrange(dplyr::pick(all_of(c(super_ranks, "seq_idx")))) |>
    dplyr::select(-any_of("seq_idx_in")) |>
    dplyr::group_by(dplyr::pick(all_of(parent_rank))) |>
    dplyr::filter(dplyr::n() > 1) |>
    targets::tar_group()
  ensure_table(out)
}

#' Calculate a pre-denovo taxon table for taxa with many sequences
#'
#' This calculates a pre-denovo taxon table for taxa with many sequences,
#' which may be efficiently parallelized within a single parent taxon.
#'
#' @param closedref_taxon_table (`data.frame`) a table of sequences for which
#' closed-reference clustering is required. Must have columns `seq_idx`
#' (`integer`) giving the indices of sequences in `seq_file`, as well as
#' columns with names matching `rank` and `parent_rank`.
#' @inheritParams full_preclosed_taxon_table
#' @param min_ops (`integer` scalar) the minimum number of operations required
#' for a single parent taxon to be considered "large"
#'
#' @return (`data.frame`) a table of sequences for which denovo clustering is
#' required. This includes both sequences which are not assigned to the target
#' rank which will act as queries, and identified sequences which belong to the
#' same parent rank as one or more queries. The results includes a `tar_group`
#' column for dynamic branching over each parent taxon. If there are no
#' queries, a single row is returned to branch over.
#' @export
large_predenovo_taxon_table <- function(
  closedref_taxon_table,
  asv_taxsort,
  rank,
  parent_rank,
  tax_ranks,
  min_ops = 1e6
) {
  checkmate::assert_number(min_ops, lower = 1)

  fptt <- full_predenovo_taxon_table(
    closedref_taxon_table,
    asv_taxsort,
    rank,
    parent_rank,
    tax_ranks
  )

  out <- fptt |>
    dplyr::group_by(dplyr::pick(all_of(parent_rank))) |>
    dplyr::filter(dplyr::n() * (dplyr::n() - 1) / 2 > min_ops) |>
    targets::tar_group()

  ensure_table(out)
}

#' Calculate a pre-denovo taxon table for taxa with few sequences
#'
#' This calculates a pre-denovo taxon table for taxa with few sequences,
#' and groups them into approximately equally sized execution groups.
#'
#' @param closedref_taxon_table (`data.frame`) a table of sequences for which
#' closed-reference clustering is required. Must have columns `seq_idx`
#' (`integer`) giving the indices of sequences in `seq_file`, as well as
#' columns with names matching `rank` and `parent_rank`.
#' @inheritParams full_preclosed_taxon_table
#' @param max_ops (`integer` scalar) the maximum number of operations to put
#' together in a single execution group
#'
#' @return (`data.frame`) a table of sequences for which denovo clustering is
#' required. This includes both sequences which are not assigned to the target
#' rank which will act as queries, and identified sequences which belong to the
#' same parent rank as one or more queries. The results includes a `tar_group`
#' column for dynamic branching over execution groups.
#' @export
small_predenovo_taxon_table <- function(
  closedref_taxon_table,
  asv_taxsort,
  rank,
  parent_rank,
  tax_ranks,
  max_ops = 1e6
) {
  checkmate::assert_number(max_ops, lower = 1)

  fptt <- full_predenovo_taxon_table(
    closedref_taxon_table,
    asv_taxsort,
    rank,
    parent_rank,
    tax_ranks
  )

  out <- fptt |>
    dplyr::group_by(dplyr::pick(all_of(parent_rank))) |>
    dplyr::mutate(ops = dplyr::n() * (dplyr::n() - 1) / 2) |>
    dplyr::filter(ops <= max_ops) |>
    dplyr::select(-all_of("tar_group")) |>
    # Ungroup so distribute_tasks() sees all taxa, not one group at a time.
    dplyr::ungroup()

  ops <- dplyr::select(out, all_of(parent_rank), ops) |>
    dplyr::distinct() |>
    dplyr::mutate(tar_group = distribute_tasks(ops, max_ops)) |>
    dplyr::select(-ops)

  dplyr::left_join(out, ops, by = parent_rank) |>
    dplyr::select(-ops) |>
    ensure_table()
}

#' Run denovo clustering
#'
#' This function calls [optimotu::seq_cluster()] for one or more
#' parent taxa.
#'
#' @param predenovo_taxon_table (`data.frame`) a table of sequences for which
#' denovo clustering is required. Must have columns `seq_idx` (`integer`)
#' giving the indices of sequences in `seq_file`, as well a column with
#' names matching `parent_rank`.
#' @param seq_file (`character`) the path to the sequence file
#' @param seq_file_index (`character` or `fastqindexr_index`) file path or
#'  [`fastqindexr_index`][fastqindexr::create_index()] object for the sequence
#'  file index
#' @param rank (`character`) the (first) taxonomic rank to calculate
#' @param parent_rank (`character`) the taxonomic rank to use as a parent
#' @param tax_ranks (`character`) the taxonomic ranks in use, from most- to
#' least-inclusive (e.g., from "kingdom" to "species")
#' @param denovo_thresholds (`list`) a list of thresholds to use for
#' denovo clustering. These are interpreted as distances on the
#' interval \[0, 100\], where 0 is a perfect match and 100 is a perfect
#' mismatch.
#' @param dist_config (`list`) a list of configuration options for the
#' distance calculation, as returned by [optimotu::dist_config()].
#' @param clust_config (`list`) a list of configuration options for the
#' clustering, as returned by [optimotu::clust_config()].
#' @param parallel_config (`list`) a list of configuration options for the
#' parallelization, as returned by [optimotu::parallel_config()].
#'
#' @return (`data.frame`) a table of taxonomic assignments with columns `seq_id`
#' and one column for each rank in `tax_ranks`; ranks from `parent_rank` and up
#' will be `character` vectors, while ranks from `rank` down will be `integer`
#' vectors (representing unnamed cluster indices)
#' @export
do_denovo_cluster <- function(
  predenovo_taxon_table,
  seq_file,
  seq_file_index,
  rank,
  parent_rank,
  tax_ranks,
  denovo_thresholds,
  dist_config,
  clust_config = optimotu::clust_tree(),
  parallel_config = optimotu::parallel_concurrent(
    optimotu.pipeline::local_cpus()
  )
) {
  super_ranks <- superranks(rank, tax_ranks)
  sub_ranks <- subranks(rank, tax_ranks)
  checkmate::assert_file_exists(seq_file, "r")
  checkmate::assert(
    checkmate::check_file_exists(seq_file_index, "r"),
    checkmate::check_class(seq_file_index, "fastqindexr_index"),
    combine = "or"
  )

  if (nrow(predenovo_taxon_table) <= 1L) {
    return(
      c(
        c("seq_id", super_ranks) |>
          (\(x) `names<-`(x, x))() |>
          purrr::map(~ character(0)),
        c(rank, sub_ranks) |>
          (\(x) `names<-`(x, x))() |>
          purrr::map(~ integer(0))
      ) |>
        tibble::as_tibble()
    )
  }

  if (dplyr::n_distinct(predenovo_taxon_table[[parent_rank]]) > 1L) {
    out <- predenovo_taxon_table |>
      dplyr::group_by(dplyr::pick(all_of(super_ranks))) |>
      dplyr::group_modify(
        function(x, y) {
          do_denovo_cluster(
            predenovo_taxon_table = x,
            seq_file = seq_file,
            seq_file_index = seq_file_index,
            rank = rank,
            parent_rank = parent_rank,
            tax_ranks = tax_ranks,
            denovo_thresholds = denovo_thresholds,
            dist_config = dist_config,
            clust_config = clust_config,
            parallel_config = parallel_config
          ) |>
            dplyr::select(-all_of(super_ranks))
        },
        .keep = TRUE
      ) |>
      dplyr::ungroup()
    return(out)
  }

  denovo_seq <- fastqindexr::extract_sequences(
    index = seq_file_index,
    seq_idx = predenovo_taxon_table$seq_idx,
    file = seq_file,
    return = "seq"
  )
  threshold_config <- optimotu::threshold_set(
    tryCatch(
      denovo_thresholds[[unique(predenovo_taxon_table[[parent_rank]])]],
      error = function(e) denovo_thresholds[["_NA_"]]
    ) |>
      optimotu::threshold_as_dist()
  )
  # seq_cluster() with USEARCH rewrites the call to
  # seq_cluster_usearch.DNAStringSet and evals it in the caller, where
  # that unexported method name is not visible. The exported generic
  # with a DNAStringSet dispatches through S3 instead.
  if (identical(dist_config$method, "usearch")) {
    denovo_sc <- optimotu::seq_cluster_usearch(
      seq = Biostrings::DNAStringSet(denovo_seq),
      threshold_config = threshold_config,
      clust_config = clust_config,
      parallel_config = parallel_config,
      usearch = dist_config$usearch,
      usearch_ncpu = dist_config$usearch_ncpu
    )
  } else {
    denovo_sc <- optimotu::seq_cluster(
      seq = denovo_seq,
      dist_config = dist_config,
      threshold_config = threshold_config,
      clust_config = clust_config,
      parallel_config = parallel_config
    )
  }
  denovo_sc |>
    t() |>
    tibble::as_tibble() |>
    dplyr::bind_cols(
      predenovo_taxon_table |>
        dplyr::select(-any_of(c(rank, "tar_group", "seq_idx"))),
      . = _
    )
}

#' Validate `asv_taxsort` column layout
#'
#' @param asv_taxsort (`data.frame`) ASV taxonomic sort table.
#' @return (`character`) join key: `"seq_idx_in"` or `"seq_id"`.
#' @keywords internal
validate_asv_taxsort <- function(asv_taxsort) {
  checkmate::assert_data_frame(asv_taxsort)
  checkmate::assert_names(names(asv_taxsort), must.include = "seq_idx")
  has_idx_in <- "seq_idx_in" %in% names(asv_taxsort)
  has_id <- "seq_id" %in% names(asv_taxsort)
  if (has_idx_in && has_id) {
    stop(
      "`asv_taxsort` must not contain both `seq_idx_in` and `seq_id`."
    )
  }
  if (!has_idx_in && !has_id) {
    stop(
      "`asv_taxsort` must contain exactly one of `seq_idx_in` or",
      "`seq_id` (in addition to `seq_idx`)."
    )
  }
  checkmate::assert_integerish(asv_taxsort[["seq_idx"]])
  if (has_idx_in) {
    checkmate::assert_integerish(asv_taxsort[["seq_idx_in"]])
    return("seq_idx_in")
  }
  checkmate::assert_character(asv_taxsort[["seq_id"]])
  "seq_id"
}

#' Prepare a taxon table for joining with `asv_taxsort`
#'
#' @param taxon_table (`data.frame`) known or closed-reference taxon table.
#' @param join_by (`character`) `"seq_idx_in"` or `"seq_id"`.
#' @return (`data.frame`) table ready for joining with `asv_taxsort`.
#' @keywords internal
prep_table_for_taxsort <- function(taxon_table, join_by) {
  seq_idx <- seq_id <- NULL
  if (join_by == "seq_idx_in") {
    if ("seq_idx" %in% names(taxon_table)) {
      return(dplyr::rename(taxon_table, seq_idx_in = seq_idx))
    }
    if ("seq_id" %in% names(taxon_table)) {
      return(
        dplyr::mutate(
          taxon_table,
          seq_idx_in = readr::parse_number(seq_id)
        )
      )
    }
    stop(
      "taxon table must contain `seq_idx` or `seq_id` when joining on",
      "`seq_idx_in`."
    )
  }
  if (
    "seq_idx" %in% names(taxon_table) && !("seq_id" %in% names(taxon_table))
  ) {
    stop(
      "When `asv_taxsort` uses `seq_id`, the taxon table must include",
      "`seq_id`; a taxon table with only `seq_idx` is not allowed."
    )
  }
  if (!"seq_id" %in% names(taxon_table)) {
    stop(
      "taxon table must contain `seq_id` when joining on `seq_id`."
    )
  }
  taxon_table
}

#' Join a taxon table with `asv_taxsort`
#'
#' @param taxon_table (`data.frame`) known or closed-reference taxon table.
#' @param asv_taxsort (`data.frame`) ASV taxonomic sort table.
#' @return (`data.frame`) joined table.
#' @keywords internal
join_table_with_taxsort <- function(taxon_table, asv_taxsort) {
  join_by <- validate_asv_taxsort(asv_taxsort)
  taxon_table <- prep_table_for_taxsort(taxon_table, join_by)
  dplyr::left_join(taxon_table, asv_taxsort, by = join_by)
}
