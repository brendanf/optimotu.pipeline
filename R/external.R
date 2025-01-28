#' Find an executable
#'
#' This function tries to find an executable in the system. It first checks the
#' environment variables, then the system path, and finally the current working
#' directory.
#'
#' @param executable (`character` string) the name of the executable to find.
#' @return (`character` string) the full path to the executable.
#' @export

find_executable <- function(executable) {
  checkmate::assert_character(executable)
  out <- Sys.getenv(executable)
  if (nchar(out) == 0 || !file.exists(out)) {
    out <- Sys.getenv(toupper(executable))
  }
  if (nchar(out) == 0 || !file.exists(out)) {
    out <- Sys.which(executable)
  }
  if (nchar(out) == 0 || !file.exists(out)) {
    out <- list.files(path = "bin", pattern = executable, recursive = TRUE, full.names = TRUE)
  }
  checkmate::assert_file_exists(out, access = "x", .var.name = executable)
  out
}
