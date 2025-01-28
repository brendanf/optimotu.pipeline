

#### targets metaprogramming ####

#' Get all the target names defined in a plan
#' @param plan (`targets::tar_target()` object or (possibly nested) list of such
#' objects) the plan to extract target names from
#' @return a `character` vector of target names
#' @export
get_target_names <- function(plan) {
  if (methods::is(plan, "tar_target")) {
    plan$settings$name
  } else {
    unname(unlist(lapply(plan, get_target_names)))
  }
}

#' Get variants of a target name which has been run through "tar_map"
#' @param plan (a named nested `list` of `targets::tar_target()` objects, as
#' produced by `targets::tar_map()`) the targets plan to extract target names
#' from
#' @param target_name (`character` string) the pre-mapping name of the target
#' to extract variants of
#' @return a `list` of `symbol` giving the names of the target variants
tar_map_symbols <- function(plan, target_name = NULL) {
  if (!is.null(target_name)) plan <- plan[[target_name]]
  rlang::syms(tarchetypes::tar_select_names(plan, everything()))
}

#' Generate quosures which combines static branching targets
#' @param plan (a named, nested `list` of `targets::tar_target()` objects, as
#' produced by `targets::tar_map()`) the targets plan to operate on
#' @param target_name (`character` string) the pre-mapping name of the target
#' to combine
#' @return a `quosure` which combines the targets
#' @export
tar_map_bind_rows <- function(plan, target_name = NULL) {
  rlang::quo(
    dplyr::bind_rows(
      !!!tar_map_symbols(plan, target_name)
    )
  )
}

#' Generate a quosure which combines static branching targets with `vctrs::vec_c()`
#' @rdname tar_map_bind_rows
#' @export
tar_map_c <- function(plan, target_name = NULL) {
  rlang::quo(
    vctrs::vec_c(
      !!!tar_map_symbols(plan, target_name)
    )
  )
}

#' Generate a quosure which combines static branching targets with `list()`
#' @rdname tar_map_bind_rows
#' @export
tar_map_list <- function(plan, target_name = NULL) {
  rlang::quo(
    list(
      !!!tar_map_symbols(plan, target_name)
    )
  )
}
