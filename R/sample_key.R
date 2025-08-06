
#' Convert a file name to a sample key
#' @param filename (`character`) the name of a file, possibly including path
#' @return (`character`) the sample key
#' @export
file_to_sample_key <- function(filename) {
  sub("_(fwd|rev)_R[12](_(filt|trim))?\\.fastq\\.gz", "", basename(filename))
}
