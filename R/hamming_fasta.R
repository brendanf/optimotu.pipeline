#' Packed Hamming width of one sequence string
#'
#' Counts characters other than lowercase `[a-z]` (A2M insertions), matching
#' `PackedSequenceSet` in optimotu.
#'
#' @param seq (`character(1)`) ungapped or aligned sequence (newlines already
#'   removed)
#' @return (`integer(1)`) packed width
#' @keywords internal
#' @noRd
hamming_packed_width <- function(seq) {
  nchar(gsub("[a-z]", "", seq, perl = TRUE), type = "chars")
}

#' Open a FASTA file for text streaming (gzip or plain)
#' @param path (`character(1)`) file path
#' @return a connection
#' @keywords internal
#' @noRd
fasta_text_connection <- function(path) {
  if (grepl("[.]gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
}

#' Whether a FASTA is Hamming-ready (shared packed width)
#'
#' Streams `path` and computes packed width per record (characters other
#' than lowercase `[a-z]`, i.e. A2M insertions). Returns `TRUE` if every
#' sequence shares the same positive packed width. An empty file, or a
#' file with no FASTA records, is not ready. Short-circuits on the first
#' width mismatch.
#'
#' @param path (`character(1)`) FASTA path, possibly gzipped
#' @return (`logical(1)`)
#' @export
fasta_is_hamming_ready <- function(path) {
  checkmate::assert_string(path, min.chars = 1L)
  checkmate::assert_file_exists(path, access = "r")
  con <- fasta_text_connection(path)
  on.exit(close(con), add = TRUE)

  packed_ref <- NA_integer_
  n_seq <- 0L
  buf <- character()
  in_seq <- FALSE

  flush_seq <- function() {
    if (!isTRUE(in_seq)) {
      return(TRUE)
    }
    width <- hamming_packed_width(paste(buf, collapse = ""))
    n_seq <<- n_seq + 1L
    if (is.na(packed_ref)) {
      packed_ref <<- width
      return(TRUE)
    }
    identical(width, packed_ref)
  }

  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) {
      break
    }
    if (startsWith(line, ">")) {
      if (!flush_seq()) {
        return(FALSE)
      }
      buf <- character()
      in_seq <- TRUE
    } else if (isTRUE(in_seq)) {
      buf <- c(buf, line)
    }
  }
  if (!flush_seq()) {
    return(FALSE)
  }
  n_seq > 0L && !is.na(packed_ref) && packed_ref > 0L
}
