#' Splits a (possibly gzipped) fastx file into subfiles
#'
#' Can handle even large files with minimal memory.
#' @param infile (`character(1)`) path to input file
#' @param n (`integer(1)`) number of subfiles to split into
#' @param outroot (`character(1)`) path to directory to write subfiles to
#' @param compress (`logical(1)`) whether to compress output files
#' @return (`character(n)`) paths to the output files
#' @export
fastx_split <- function(infile, n, outroot = tempfile(), compress = FALSE) {
  checkmate::assert_string(infile)
  checkmate::assert_file(infile, access = "r")
  checkmate::assert_int(n, lower = 1, upper = 64)
  checkmate::assert_path_for_output(outroot)
  checkmate::assert_flag(compress)

  is_fastq <- grepl(fastq_regex, infile)
  is_gz <- endsWith(infile, ".gz")

  if (n == 1 && is_gz == compress) return(infile)

  suffix <- if (is_fastq) ".fastq" else ".fasta"
  if (compress) suffix <- paste0(suffix, ".gz")

  outfiles <- replicate(n, tempfile(fileext = suffix, tmpdir = outroot))
  if (is_fastq) {
    fastq_split(infile, outfiles, compress = compress)
  } else {
    fasta_split(infile, outfiles, compress = compress)
  }
}

#' Rejoins some split fastx files
#'
#' Of course this is most useful if something happened to them in between.
#' If no lines are missing or reordered in the infiles, then the outfile
#' will end up in the same order as the original file which was split.
#' Otherwise this will almost certainly not happen.
#' @param infiles (`character(n)`) paths to the input files
#' @param outfile (`character(1)`) path to the output file
#' @return (`character(1)`) path to the output file
#' @export
fastx_combine <- function(infiles, outfile) {
  checkmate::assert_file(infiles, "r")
  checkmate::assert_path_for_output(outfile, overwrite = TRUE)
  is_fastq <- grepl(fastq_regex, infiles)
  stopifnot(all(is_fastq) | all(!is_fastq))
  is_fastq <- all(is_fastq)

  is_gz <-endsWith(infiles, ".gz")
  stopifnot(all(is_gz) | all(!is_gz))
  is_gz <- all(is_gz)

  compress <- endsWith(outfile, ".gz")

  if (is_fastq) {
    fastq_combine(infiles, outfile, compress = compress)
  } else {
    fasta_combine(infiles, outfile, compress = compress)
  }
}
