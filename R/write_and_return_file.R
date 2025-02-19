#### write a file and returning its name ####

#' Create the parent directory of a file if it does not exist
#' @param file (`character` string) file path to create the parent directory for
#' @return the input file path (invisibly)
#' @export
ensure_directory <- function(file) {
  d <- dirname(file)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  invisible(file)
}

#' Write an object to a file and return the file path
#' @param x (any object) object to write
#' @param file (`character` string) file path to write to
#' @param ... additional arguments to pass to the writing function
write_and_return_file <- function(x, file, ...) {
  UseMethod("write_and_return_file")
}

#' @rdname write_and_return_file
#' @param width (integer) the maximum width of each line in the output file
#' @exportS3Method
write_and_return_file.XStringSet <- function(x, file, width = 20001L, ...) {
  ensure_directory(file)
  Biostrings::writeXStringSet(x, file, width = width, ...)
  file
}

#' @rdname write_and_return_file
#' @exportS3Method
#' @param type (`character` string) the type of file to write to; currently
#'   supported are `"rds"` (default) and `"tsv"`
write_and_return_file.data.frame <- function(x, file, type = c("rds", "tsv"), ...) {
  ensure_directory(file)
  type = match.arg(type)
  switch(
    type,
    rds = saveRDS(x, file, ...),
    tsv = readr::write_tsv(x, file, ...),
    stop("Unknown file type: ", type)
  )
  file
}

#' @rdname write_and_return_file
#' @exportS3Method
write_and_return_file.character <- function(x, file, ...) {
  ensure_directory(file)
  writeLines(x, file, ...)
  file
}

#' @rdname write_and_return_file
#' @exportS3Method
write_and_return_file.ggplot <- function(x, file, ...) {
  ensure_directory(file)
  ggplot2::ggsave(file, plot = x, ...)
  file
}

#' @rdname write_and_return_file
#' @exportS3Method
write_and_return_file.default <- function(x, file, ...) {
  ensure_directory(file)
  saveRDS(x, file, ...)
  file
}
