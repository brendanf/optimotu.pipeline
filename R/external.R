#' Find an executable
#'
#' This function tries to find an executable in the system. It first checks
#' any path configured under `executables:` in `pipeline_options.yaml`, then
#' the environment variables, then the system path, and finally the current
#' working directory / `OPTIMOTU_BIN_DIR`.
#'
#' Because resolution depends on R options set by
#' [parse_pipeline_options()], call this on the main process at plan
#' definition time and inject the result into target commands with `!!`
#' (for example `vsearch = !!find_vsearch()`). Do not rely on default
#' arguments that call `find_*()` to run on `crew` workers.
#'
#' @param executable (`character` string) the name of the executable to find.
#' @return (`character` string) the full path to the executable.
#' @export

find_executable <- function(executable) {
  checkmate::assert_character(executable)
  configured <- configured_executables()
  if (executable %in% names(configured)) {
    configured_path <- unname(configured[[executable]])
    if (checkmate::test_file_exists(configured_path, access = "x")) {
      return(normalizePath(configured_path, mustWork = TRUE))
    }
    which_configured <- unname(Sys.which(configured_path))
    if (nzchar(which_configured) && file.exists(which_configured)) {
      return(which_configured)
    }
  }
  out <- Sys.getenv(executable)
  if (nchar(out) == 0 || !file.exists(out)) {
    out <- Sys.getenv(toupper(executable))
  }
  if (nchar(out) == 0 || !file.exists(out)) {
    out <- unname(Sys.which(executable))
  }
  if (nchar(out) == 0 || !file.exists(out)) {
    bindir <- Sys.getenv("OPTIMOTU_BIN_DIR", unset = "bin")
    out <- list.files(
      path = bindir,
      pattern = executable,
      recursive = TRUE,
      full.names = TRUE
    )
  }
  checkmate::assert_file_exists(out, access = "x", .var.name = executable)
  unname(out)
}
