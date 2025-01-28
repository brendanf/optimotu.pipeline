#' Try to find the usearch executable
#' @return (`character` string) the full path to the usearch executable.
#' @export
find_usearch <- function() {
  find_executable("usearch")
}

#' Build a usearch database (UDB) file using USEARCH
#' @param infile (`character` string) input file name
#' @param outfile (`character` string) output file name
#' @param type (`character` string) type of database to build
#' @param usearch (`character` string) path to the usearch executable
#' @return (`character` string) the name of the output file
#' @export
build_udb <- function(infile, outfile, type = c("usearch", "sintax", "ublast"),
                      usearch = find_usearch()) {
  type <- match.arg(type)
  command <- paste0("-makeudb_", type)
  args <- c(
    command, infile,
    "-output", outfile
  )
  result <- system2(usearch, args)
  stopifnot(result == 0)
  outfile
}

#' Build a usearch database (UDB) file using USEARCH, removing blacklisted
#' sequences
#' @inheritParams build_udb
#' @param blacklist (`character` vector) names  of sequences to remove from the
#' database
#' @return (`character` string) the name of the output file
#' @export
build_filtered_udb <- function(
    infile,
    outfile,
    type = c("usearch", "sintax", "ublast"),
    blacklist,
    usearch = find_usearch()
) {
  # make sure we have valid arguments
  type <- match.arg(type)
  command <- paste0("-makeudb_", type)
  stopifnot(system2(usearch, "--version") == 0)

  # make a temp file and a temp fifo
  blf <- withr::local_tempfile(fileext = ".txt")
  tf <- withr::local_tempfile(fileext = ".fasta")
  writeLines(blacklist, blf)
  stopifnot(system2("mkfifo", tf) == 0)

  # first usearch call removes the blacklisted sequences
  system2(
    usearch,
    args = c(
      "--fastx_getseqs", infile,
      "--labels", blf,
      "--label_substr_match",
      "--notmatched", tf
    ),
    wait = FALSE
  )
  # second usearch call creates the udb file
  result = system2(
    usearch,
    args = c(
      command, tf,
      "--output", outfile
    )
  )
  stopifnot(result == 0)
  outfile
}
