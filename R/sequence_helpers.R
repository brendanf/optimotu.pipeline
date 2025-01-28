#### generic sequence helpers ####
# helper functions to work with sequences sets that may be XStringSet,
# (named) character, fasta/fastq file, or data.frame

fasta_regex <- "\\.(fas?|fasta)(\\.gz)?$"
fastq_regex <- "\\.f(ast)?q(\\.gz)?$"

#' Guess the column name in a data frame which refers to the sequence ID
#' @param d (`data.frame`) data frame to search
#' @return (`character`) name of the column
#' @export
find_name_col <- function(d) {
  stopifnot(is.data.frame(d))
  if ("seq_id" %in% names(d)) return("seq_id")
  if ("name" %in% names(d)) return("name")
  if ("ASV" %in% names(d)) return("ASV")
  if ("OTU" %in% names(d)) return("OTU")
  if (ncol(d) == 2) {
    if ("seq" %in% names(d)) return(names(d)[names(d) != "seq"])
    if ("sequence" %in% names(d)) return(names(d)[names(d) != "sequence"])
    return(names(d)[1])
  }
  stop("unable to determine sequence name column:", names(d))
}

#' Guess the column name in a data frame which refers to the sequence itself
#' @param d (`data.frame`) data frame to search
#' @return (`character`) name of the column
#' @export
find_seq_col <- function(d) {
  stopifnot(is.data.frame(d))
  if ("seq" %in% names(d)) return("seq")
  if ("sequence" %in% names(d)) return("sequence")
  if (ncol(d) == 2) {
    if ("seq_id" %in% names(d)) return(names(d)[names(d) != "seq_id"])
    if ("name" %in% names(d)) return(names(d)[names(d) != "name"])
    if ("ASV" %in% names(d)) return(names(d)[names(d) != "ASV"])
    if ("OTU" %in% names(d)) return(names(d)[names(d) != "OTU"])
    return(names(d)[2])
  }
  stop("unable to determine sequence column:", names(d))
}

#' Write sequences as a fasta file
#'
#' @param seq (`data.frame`, `character`, or `XStringSet`) sequences to write
#' @param fname (`character` string) file path to write to
#' @param ... additional arguments to pass to the writing function
#' @export
write_sequence <- function(seq, fname, ...) {
  UseMethod("write_sequence", seq)
}

#' @rdname write_sequence
#' @exportS3Method
#' @param seq_col (`character` string) name of the column in `seq` which
#' contains the sequences
#' @param name_col (`character` string) name of the column in `seq` which
#' contains the sequence names
write_sequence.data.frame <- function(seq, fname, seq_col = find_seq_col(seq),
                                      name_col = find_name_col(seq), ...) {
  dplyr::select(seq, !!name_col, !!seq_col) |>
    tibble::deframe() |>
    Biostrings::DNAStringSet() |>
    write_and_return_file(fname, ...)
}

#' @rdname write_sequence
#' @exportS3Method
write_sequence.character <- function(seq, fname, ...) {
  Biostrings::DNAStringSet(seq) |>
    write_and_return_file(fname, ...)
}

#' @rdname write_sequence
#' @exportS3Method
write_sequence.XStringSet <- function(seq, fname, ...) {
  write_and_return_file(seq, fname, ...)
}

#' Select elements from a sequence set
#' @param seq (`data.frame`, `character` string representing a file name, or
#' any object with named elements) sequences to select from
#' @param which (`integer`, `logical`, or `character`) elements to select
#' @param negate (`logical`) if TRUE, select all elements except those specified
#' @param ... additional arguments to pass to the selection function
#' @return (same type of object as `seq`) selected elements
#' @export
select_sequence <- function(seq, which, negate = FALSE, ...) {
  UseMethod("select_sequence", seq)
}

#' @rdname select_sequence
#' @exportS3Method select_sequence data.frame
#' @param name_col (`character` string) name of the column in `seq` which
#' contains the sequence names
#' @importFrom rlang .data
select_sequence.data.frame <- function(seq, which, negate = FALSE, name_col = find_name_col(seq), ...) {
  if (isTRUE(negate)) {
    if (is.integer(which)  ||
        (is.numeric(which) && all(which == round(which)))) return(seq[-which,])
    if (is.logical(which)) return(seq[!which,])
    if (is.character(which)) return(dplyr::filter(seq, !.data[[name_col]] %in% which))
    stop("'which' should be integer, logical, or character")
  } else if (isFALSE(negate)) {
    if (is.integer(which)  ||
        (is.numeric(which) && all(which == round(which)))) return(seq[which,])
    if (is.logical(which)) return(seq[which,])
    if (is.character(which)) return(dplyr::filter(seq, .data[[name_col]] %in% which))
    stop("'which' should be integer, logical, or character")
  }
  stop("'negate' must be TRUE or FALSE")
}

#' @rdname select_sequence
#' @exportS3Method
select_sequence.character <- function(seq, which, negate = FALSE, ...) {
  if (length(seq) == 1 && file.exists(seq)) {
    seq <- as.character(Biostrings::readBStringSet(seq))
  }
  select_sequence.default(seq, which = which, negate = negate, ...)
}

#' @rdname select_sequence
#' @exportS3Method
select_sequence.default <- function(seq, which, negate = FALSE, ...) {
  if (isTRUE(negate)) {
    if (is.integer(which) || (is.numeric(which) && all(which == round(which)))) return(seq[-which])
    if (is.logical(which)) return(seq[!which])
    if (is.character(which)) return(seq[setdiff(names(seq), which)])
    stop("'which' should be integer, logical, or character")
  } else if (isFALSE(negate)) {
    return(seq[which])
  }
  stop("'negate' must be true or false")
}

#' Get the number of sequences in a sequence set
#' @param seq (`data.frame`, `character` string representing a file name,
#' [`XStringSet`][Biostrings::XStringSet-class], or any object with multiple elements)
#' sequences to count
#' @param ... additional arguments to pass to methods
#' @return (`integer`) number of sequences
#' @export
sequence_size <- function(seq, ...) {
  UseMethod("sequence_size", seq)
}

#' @rdname sequence_size
#' @exportS3Method
sequence_size.XStringSet <- function(seq, ...) {
  length(seq)
}

#' @rdname sequence_size
#' @exportS3Method
sequence_size.character <- function(seq, ...) {
  if (length(seq) > 0 && all(file.exists(seq))) {
    if (all(grepl(fastq_regex, seq))) {
      return(
        lapply(seq, Biostrings::fastq.seqlengths) |>
          vapply(length, 1L)
      )
    } else if (all(grepl(fasta_regex, seq))) {
      return(
        lapply(seq, Biostrings::fasta.seqlengths) |>
          vapply(length, 1L)
      )
    }
  }
  sequence_size.default(seq, ...)
}

#' @rdname sequence_size
#' @exportS3Method
sequence_size.default <- function(seq, ...) {
  vctrs::vec_size(seq)
}

seqhash <- function(...) {
  digest::getVDigest("spookyhash")(...)
}

#' Hash sequences
#' @param seq (`data.frame`, `character` string representing a file name,
#' named `character` vector, or `Biostrings::XStringSet`) sequences to hash
#' @param use_names (`logical`) if `TRUE`, the result hashes will be named using
#' the names of the sequences.
#' @param ... additional arguments to pass to the hashing function
#' @return (`integer` vector) hash codes
#' @export
hash_sequences <- function(seq, use_names = TRUE, ...) {
  UseMethod("hash_sequences", seq)
}

#' @rdname hash_sequences
#' @exportS3Method
hash_sequences.character <- function(seq, use_names = TRUE, ...) {
  if (length(seq) == 1 && file.exists(seq)) {
    # This requires loading everything in memory; would be better to do a
    # batch thing
    if (grepl(fastq_regex, seq)) {
      hash_sequences(Biostrings::readQualityScaledDNAStringSet(seq), use_names, ...)
    } else if (all(grepl(fasta_regex, seq))) {
      hash_sequences(Biostrings::readBStringSet(seq), use_names, ...)
    } else {
      stop("Cannot determine file type for ", seq)
    }
  } else {
    out <- seqhash(seq)
    if (isTRUE(use_names)) {
      names(out) <- names(seq)
    }
    out
  }
}

#' @rdname hash_sequences
#' @exportS3Method
hash_sequences.XStringSet <- function(seq, use_names = TRUE, ...) {
  out <- seqhash(as.character(seq))
  if (isTRUE(use_names)) {
    names(out) <- names(seq)
  }
  out
}

#' @rdname hash_sequences
#' @exportS3Method
#' @param seq_col (`character` string) name of the column in `seq` which
#' contains the sequences
#' @param name_col (`character` string) name of the column in `seq` which
#' contains the sequence names
hash_sequences.data.frame <- function(
    seq,
    use_names = TRUE,
    seq_col = find_seq_col(seq),
    name_col = if (isTRUE(use_names)) find_name_col(seq) else NULL,
    ...
) {
  out <- seqhash(seq[[seq_col]])
  if (isTRUE(use_names)) {
    names(out) <- seq[[name_col]]
  }
  out
}

#### sequence naming ####

# Internal function; force a string to be ASCII
ascii_clean <- function(s) {
  gsub(
    s,
    pattern = "[^A-Za-z0-9[:punct:]]",
    replacement = "x",
    useBytes = TRUE,
    perl = TRUE
  )
}

#' Generate names for a set of sequences
#'
#' Names produced are like "ASV0001", "ASV0002", ... The prefix ("ASV" in
#' the example) is customizable, and the number will be zero-padded to fit
#' all values.
#'
#' @param n (`integer`) number of sequences to name
#' @param prefix (`character`) prefix to use for the names
#' @return (`character` vector) names
#' @export

make_seq_names <- function(n, prefix) {
  sprintf(
    sprintf("%s%%0%dd", prefix, max(floor(log10(n)) + 1L, 0L)),
    seq_len(n)
  )
}

#' Name a set of sequences
#'
#' Sequences are named using `make_seq_names()`.
#'
#' If the input is a file name, the file will be modified.
#'
#' @param seq (`data.frame`, `character`, or `Biostrings::XStringSet`) sequences
#' to name
#' @param ... additional arguments to pass to methods
#' @inheritParams make_seq_names
#' @return (same type of object as `seq`) named sequences
#' @export
name_seqs <- function(seq, prefix, ...) {
  UseMethod("name_seqs", seq)
}

#' @rdname name_seqs
#' @exportS3Method
name_seqs.XStringSet <- function(seq, prefix, ...) {
  names(seq) <- make_seq_names(length(seq), prefix)
  seq
}

#' @rdname name_seqs
#' @exportS3Method
name_seqs.character <- function(seq, prefix, ...) {
  if (length(seq) == 1 && file.exists(seq)) {
    width <- floor(log10(sequence_size(seq))) + 1
    tf <- withr::local_tempfile()
    file.copy(seq, tf)
    if (grepl(fasta_regex, seq)) {
      # don't trust line counts in fasta
      command <- sprintf(
        "awk '/^>/{printf(\">%s%%0%ii\\n\", ++i); next}; {print}'",
        prefix,
        width
      )
    } else if  (grepl(fastq_regex, seq)) {
      command <- sprintf(
        "awk 'NR%4==1{printf(\"@%s%%0%ii\\n\", ++i); next}; {print}'",
        prefix,
        width
      )
    } else {
      stop("Cannot determine file type for ", seq)
    }
    if (endsWith(seq, ".gz")) {
      command <- paste("zcat", tf, "|", command, "| gzip -c - >", seq)
    } else {
      command <- paste(command, "<", tf, ">", seq)
    }
    result <- system(command)
    stopifnot(result == 0L)
  } else {
    names(seq) <- make_seq_names(length(seq), prefix)
  }
  seq
}

#' @rdname name_seqs
#' @exportS3Method
#' @param id_col (`character` string) name of the column in `seq` which
#' contains the sequence IDs
name_seqs.data.frame <- function(seq, prefix, id_col = prefix, ...) {
  seq[[id_col]] <- make_seq_names(nrow(seq), prefix)
  seq
}

#' @rdname name_seqs
#' @exportS3Method
name_seqs.matrix <- function(seq, prefix, ...) {
  colnames(seq) <- make_seq_names(ncol(seq), prefix)
  seq
}
