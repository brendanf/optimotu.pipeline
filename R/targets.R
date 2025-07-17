

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

#' Substitute values into a target object or plan
#' @param target ([target][targets::tar_target()] object) the target(s) to substitute
#' @param values (`list`) a named list of values to substitute into the target
#' @return ([target][targets::tar_target()] object) the target with substituted values
#' @export
tar_substitute <- function(target, values) {
  UseMethod("tar_substitute", target)
}

#' @exportS3Method
tar_substitute.tar_target <- function(target, values) {
  checkmate::assert(
    checkmate::check_list(values, names = "unique"),
    checkmate::check_data_frame(values)
  )
  command <- do.call(substitute, list(expr = target$command$expr[[1]], env = values))
  pattern <- if (is.null(target$settings$pattern)) {
    NULL
  } else {
    tarchetypes::tar_sub_raw(
      expr = target$settings$pattern[[1]],
      values = values
    )[[1]]
  }
  targets::tar_target_raw(
    name = target$settings$name,
    command = command,
    pattern = pattern,
    packages = target$command$packages,
    library = target$command$library,
    format = target$settings$format,
    repository = target$settings$repository,
    iteration = target$settings$iteration,
    error = target$settings$error,
    memory = target$settings$memory,
    garbage_collection = target$settings$garbage_collection,
    deployment = target$settings$deployment,
    priority = target$settings$priority,
    resources = target$settings$resources,
    storage = target$settings$storage,
    retrieval = target$settings$retrieval,
    cue = targets::tar_cue(
      mode = target$cue$mode,
      command = target$cue$command,
      depend = target$cue$depend,
      format = target$cue$format,
      repository = target$cue$repository,
      iteration = target$cue$iteration,
      file = target$cue$file
    ),
    description = target$settings$description
  )
}

#' @exportS3Method
tar_substitute.list <- function(target, values) {
  lapply(target, tar_substitute, values = values)
}

#' Merge named lists within two targets (sub)plans
#'
#' This is most meaningful when the plans were both produced by `targets::tar_map()`
#' and the targets in the plans have the same names.
#'
#' @param plan1 (`targets::tar_target()` object or (possibly nested) list of such
#' objects) the first plan to merge
#' @param plan2 (`targets::tar_target()` object or (possibly nested) list of such
#' objects) the second plan to merge
#' @return a nested `list` of `targets::tar_target()` objects
#' @export
tar_merge <- function(plan1, plan2) {
  if (methods::is(plan1, "tar_target") && methods::is(plan2, "tar_target")) {
    return(list(plan1, plan2))
  }
  if (methods::is(plan1, "tar_target")) {
    return(c(list(plan1), plan2))
  }
  if (methods::is(plan2, "tar_target")) {
    return(c(plan1, list(plan2)))
  }
  if (is.null(names(plan1)) || is.null(names(plan2))) {
    return(c(plan1, plan2))
  }
  all_names <- unique(c(names(plan1), names(plan2)))
  all_names <- all_names[all_names != ""]
  names(all_names) <- all_names
  c(
    plan1[names(plan1) == ""],
    plan2[names(plan2) == ""],
    lapply(all_names, function(name) {
      tar_merge(plan1[[name]], plan2[[name]])
    })
  )
}
