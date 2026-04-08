#### write a file and returning its name ####

#' Create the parent directory of a file if it does not exist
#' @param file (`character` string) file path to create the parent directory for
#' @return the input file path (invisibly)
#' @export
ensure_directory <- function(file) {
  d <- dirname(file)
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
  invisible(file)
}

#' Write an object to a file and return the file path
#' @param x (any object) object to write
#' @param file (`character` string) file path to write to
#' @param ... additional arguments to pass to the writing function
#' @return the file path
#' @export
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
#' @exportS3Method write_and_return_file data.frame
#' @param type (`character` string) the type of file to write to; currently
#'   supported are `"rds"` (default) and `"tsv"`
write_and_return_file.data.frame <- function(
  x,
  file,
  type = c("rds", "tsv"),
  ...
) {
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
write_and_return_file.default <- function(
  x,
  file,
  type = tolower(tools::file_ext(file)),
  ...
) {
  checkmate::assert_choice(type, c("rds", "qs", "qs2", "qd", "qdata"))
  if (type == "rds") {
    ensure_directory(file)
    saveRDS(x, file, ...)
  } else if (type == "qs") {
    if (!requireNamespace("qs", quietly = TRUE)) {
      stop(
        "qs package is required to save a file with format 'qs'.",
        " Please install it using `install.packages('qs')`."
      )
    }
    ensure_directory(file)
    qs::qsave(x, file, ...)
  } else if (type == "qs2") {
    if (!requireNamespace("qs2", quietly = TRUE)) {
      stop(
        "qs2 package is required to save a file with format 'qs2'.",
        " Please install it using `install.packages('qs2')`."
      )
    }
    ensure_directory(file)
    qs2::qs_save(x, file, ...)
  } else if (type == "qdata" || type == "qd") {
    if (!requireNamespace("qs2", quietly = TRUE)) {
      stop(
        "qs2 package is required to save a file with format '",
        type,
        "'. Please install it using `install.packages('qs2')`."
      )
    }
    ensure_directory(file)
    qs2::qd_save(x, file, ...)
  }
  file
}
