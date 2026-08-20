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
  if (!is.null(target_name)) {
    plan <- plan[[target_name]]
  }
  if (length(plan) == 0) {
    return(list())
  }
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
  command <- do.call(
    substitute,
    list(expr = target$command$expr[[1]], env = values)
  )
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

# Crew remote workers set tar_runtime$target but leave tar_runtime$meta
# unset. tar_meta() / tar_read() are forbidden there, so these helpers
# use the current target's subpipeline instead.
lookup_names_in_meta <- function(deps, meta) {
  unlist(
    lapply(
      deps,
      function(x) {
        if (meta$exists_record(x)) {
          record <- meta$get_record(x)
          if (record$type %in% c("stem", "branch")) {
            record$name
          } else if (identical(record$type, "pattern")) {
            record$children
          } else {
            NULL
          }
        } else {
          NULL
        }
      }
    ),
    use.names = FALSE
  )
}

lookup_names_in_pipeline <- function(deps, pipeline) {
  unlist(
    lapply(
      deps,
      function(x) {
        if (!targets:::pipeline_exists_target(pipeline, x)) {
          return(NULL)
        }
        target <- targets:::pipeline_get_target(pipeline, x)
        if (inherits(target, "tar_pattern")) {
          targets:::target_get_children(target)
        } else if (inherits(target, "tar_builder")) {
          target$name
        } else {
          NULL
        }
      }
    ),
    use.names = FALSE
  )
}

read_runtime_target <- function(name) {
  name <- as.character(name)[[1L]]
  runtime <- targets::tar_runtime_object()
  meta <- runtime$meta
  if (!is.null(meta) && isTRUE(meta$exists_record(name))) {
    record <- meta$get_record(name)
    store <- targets:::record_bootstrap_store(record)
    file <- targets:::record_bootstrap_file(record)
    return(targets:::store_read_object(store, file))
  }
  pipeline <- runtime$target$subpipeline
  if (
    !is.null(pipeline) &&
      isTRUE(targets:::pipeline_exists_target(pipeline, name))
  ) {
    dep <- targets:::pipeline_get_target(pipeline, name)
    return(targets:::target_read_value(dep, pipeline)$object)
  }
  if (!is.null(runtime$target) || !is.null(runtime$store)) {
    stop(
      "Cannot read target '",
      name,
      "' while the pipeline is running: tar_runtime$meta is unset and ",
      "the name is not in the current target's subpipeline."
    )
  }
  targets::tar_read_raw(name)
}

#' Extract names of targets from an expression
#'
#' Resolves symbols in `expr` to pipeline target names. Inside a running
#' pipeline this uses in-memory metadata when available. Crew remote workers
#' do not populate `tar_runtime$meta`; in that case names are taken from the
#' current target's subpipeline. `targets::tar_meta()` is only used outside a
#' pipeline.
#'
#' @param expr (`symbol`, `call`, or other defused expression accepted by
#'   `targets::tar_deps_raw()`) an unevaluated expression containing symbols
#'   which refer to targets in the current pipeline
#'
#' @return (`character` vector) names of targets present in `expr`. If the
#'   targets use dynamic branching so that they are stored as multiple children,
#'   the names of the children are returned.
#' @keywords internal
extract_targets <- function(expr, ...) {
  deps <- targets::tar_deps_raw(expr)
  runtime <- targets::tar_runtime_object()
  if (!is.null(runtime$meta)) {
    return(lookup_names_in_meta(deps, runtime$meta))
  }
  pipeline <- runtime$target$subpipeline
  if (!is.null(pipeline)) {
    return(lookup_names_in_pipeline(deps, pipeline))
  }
  if (!is.null(runtime$target) || !is.null(runtime$store)) {
    stop(
      "Cannot resolve target names while the pipeline is running: ",
      "tar_runtime$meta is unset and the current target has no subpipeline."
    )
  }
  meta <- targets::tar_meta(
    any_of(deps),
    fields = c("name", "type", "children"),
    targets_only = TRUE
  )
  has_children <- !is.na(meta$children)
  c(
    meta[!has_children, "name"],
    unlist(meta[has_children, "children"])
  )
}
