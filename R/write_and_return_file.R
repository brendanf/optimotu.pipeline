#### write a file and returning its name ####

#' Tabular formats handled by [write_and_return_file.data.frame()]
#' @keywords internal
#' @noRd
.per_file_table_formats <- function() {
  c(
    "rds",
    "tsv",
    "csv",
    "xlsx",
    "fst",
    "feather",
    "parquet",
    "qs2",
    "qdata"
  )
}

#' Formats that require a data.frame coercion
#' @keywords internal
#' @noRd
.data_frame_dispatch_formats <- function() {
  c("tsv", "csv", "xlsx", "fst", "feather", "parquet")
}

#' Formats that can be written from a matrix without coercion
#' @keywords internal
#' @noRd
.matrix_direct_formats <- function() {
  setdiff(.per_file_table_formats(), .data_frame_dispatch_formats())
}

#' File extension for a tabular output format name
#' @keywords internal
#' @noRd
.format_file_extension <- function(type) {
  switch(
    tolower(type),
    qdata = "qdata",
    qd = "qdata",
    tsv = "tsv",
    csv = "csv",
    xlsx = "xlsx",
    fst = "fst",
    feather = "feather",
    parquet = "parquet",
    qs2 = "qs2",
    rds = "rds",
    tolower(type)
  )
}

#' Infer and normalize the output format for [write_and_return_file()]
#' @param file destination path
#' @param type format name, or `NULL` to use the `file` extension
#' @keywords internal
#' @noRd
.normalize_write_type <- function(file, type = NULL) {
  if (is.null(type)) {
    type <- tolower(tools::file_ext(file))
  } else {
    type <- tolower(type)
  }
  if (identical(type, "qd")) {
    type <- "qdata"
  }
  type
}

#' Require a suggested package for a file format
#' @keywords internal
#' @noRd
.require_format_package <- function(pkg, type) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      sprintf(
        "%s package is required to save a file with format '%s'.",
        pkg,
        type
      ),
      " Please install it using `install.packages('",
      pkg,
      "')`."
    )
  }
}

#' Resolve a list to a named object list for RData export
#' @param x list passed to [write_and_return_file.list()]
#' @param envir environment for symbol lookup
#' @return named list, or `NULL` if resolution failed
#' @keywords internal
#' @noRd
.resolve_list_for_rdata <- function(x, envir) {
  checkmate::assert_list(x)
  if (length(x) == 0L) {
    return(x)
  }
  if (!is.null(names(x)) && any(nzchar(names(x)))) {
    return(x)
  }
  nm <- character(length(x))
  for (i in seq_along(x)) {
    el <- x[[i]]
    if (is.character(el) && length(el) == 1L) {
      nm[i] <- el
    } else if (is.symbol(el) || inherits(el, "name")) {
      nm[i] <- as.character(el)
    } else if (
      inherits(el, "call") &&
        identical(el[[1L]], quote(quote)) &&
        length(el) == 2L
    ) {
      inner <- el[[2L]]
      if (is.symbol(inner) || inherits(inner, "name")) {
        nm[i] <- as.character(inner)
      } else {
        return(NULL)
      }
    } else {
      return(NULL)
    }
  }
  out <- mget(
    nm,
    envir = envir,
    ifnotfound = list(function(name) {
      stop(
        "Object '",
        name,
        "' not found when building RData file.",
        call. = FALSE
      )
    })
  )
  stats::setNames(out, nm)
}

#' Create the parent directory of one or more files if they do not exist
#' @param file (`character` vector) file paths to create the parent
#'   directories for
#' @return the input file path(s) (invisibly)
#' @export
ensure_directory <- function(file) {
  checkmate::assert_character(file)
  for (d in unique(dirname(file))) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE)
    }
  }
  invisible(file)
}

#' Write tabular pipeline outputs in one or more file formats
#'
#' Writes `x` to `file.<ext>` for each format in `formats` (by default
#' [output_table_formats()]). The path `file` should not include a file
#' extension.
#'
#' @param x (`data.frame`, `matrix`, or `tibble`) tabular data to write
#' @param file (`character`) output path without extension
#' @param formats (`character`) format names; see [output_formats()]
#' @param ... passed to [write_and_return_file()]
#' @return `character` vector of paths written
#' @export
write_tabular_outputs <- function(
  x,
  file,
  formats = output_table_formats(),
  ...
) {
  checkmate::assert_character(file, len = 1)
  file <- sub("\\.[^.]+$", "", file)
  checkmate::assert_character(formats, min.len = 1)
  vapply(
    formats,
    function(fmt) {
      ext <- .format_file_extension(fmt)
      write_and_return_file(
        x,
        paste0(file, ".", ext),
        type = fmt,
        ...
      )
    },
    FUN.VALUE = character(1)
  )
}

#' Write an object to a file and return the file path
#' @param x (any object) object to write
#' @param file (`character` string) file path to write to
#' @param ... additional arguments passed to methods; several methods accept
#'   a `type` argument (see method documentation)
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
#' @param type (`character`) tabular format: `"rds"`, `"tsv"`, `"csv"`,
#'   `"xlsx"`, `"fst"`, `"feather"`, `"parquet"`, `"qs2"`, or `"qdata"`
#'   (`"qd"` is an alias). If `NULL`, inferred from `file` extension (not
#'   case-sensitive).
#' @exportS3Method write_and_return_file data.frame
write_and_return_file.data.frame <- function(
  x,
  file,
  type = NULL,
  ...
) {
  type <- .normalize_write_type(file, type)
  checkmate::assert_choice(type, .per_file_table_formats())
  ensure_directory(file)
  switch(
    type,
    rds = saveRDS(x, file, ...),
    tsv = readr::write_tsv(x, file, ...),
    csv = readr::write_csv(x, file, ...),
    xlsx = {
      .require_format_package("writexl", type)
      writexl::write_xlsx(x, file, ...)
    },
    fst = {
      .require_format_package("fst", type)
      fst::write_fst(x, file, ...)
    },
    feather = {
      .require_format_package("arrow", type)
      arrow::write_feather(x, file, ...)
    },
    parquet = {
      .require_format_package("arrow", type)
      arrow::write_parquet(x, file, ...)
    },
    qs2 = {
      .require_format_package("qs2", type)
      qs2::qs_save(x, file, ...)
    },
    qdata = {
      .require_format_package("qs2", type)
      qs2::qd_save(x, file, ...)
    }
  )
  file
}

#' @rdname write_and_return_file
#' @exportS3Method write_and_return_file matrix
write_and_return_file.matrix <- function(x, file, type = NULL, ...) {
  type <- .normalize_write_type(file, type)
  if (type %in% .matrix_direct_formats()) {
    return(write_and_return_file.default(x, file, type = type, ...))
  }
  checkmate::assert_choice(type, .per_file_table_formats())
  write_and_return_file.data.frame(
    as.data.frame(x, stringsAsFactors = FALSE),
    file,
    type = type,
    ...
  )
}

#' @rdname write_and_return_file
#' @param envir (`environment`) environment for resolving symbol names when
#'   writing RData from an unnamed list of names or symbols
#' @exportS3Method write_and_return_file list
write_and_return_file.list <- function(
  x,
  file,
  type = NULL,
  envir = parent.frame(),
  ...
) {
  type <- .normalize_write_type(file, type)
  if (type %in% c("rdata", "rda")) {
    resolved <- .resolve_list_for_rdata(x, envir = envir)
    if (!is.null(resolved)) {
      ensure_directory(file)
      save(
        list = names(resolved),
        file = file,
        envir = list2env(resolved, parent = emptyenv())
      )
      return(file)
    }
  }
  write_and_return_file.default(x, file, type = type, ...)
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
  type = NULL,
  ...
) {
  type <- .normalize_write_type(file, type)
  checkmate::assert_choice(
    type,
    c("rds", "qs", "qs2", "qd", "qdata", .per_file_table_formats())
  )
  if (type %in% .data_frame_dispatch_formats()) {
    return(
      write_and_return_file.data.frame(
        as.data.frame(x, stringsAsFactors = FALSE),
        file,
        type = type,
        ...
      )
    )
  }
  if (type == "rds") {
    ensure_directory(file)
    saveRDS(x, file, ...)
  } else if (type == "qs") {
    .require_format_package("qs", type)
    ensure_directory(file)
    qs::qsave(x, file, ...)
  } else if (type == "qs2") {
    .require_format_package("qs2", type)
    ensure_directory(file)
    qs2::qs_save(x, file, ...)
  } else if (type == "qdata") {
    .require_format_package("qs2", type)
    ensure_directory(file)
    qs2::qd_save(x, file, ...)
  }
  file
}
