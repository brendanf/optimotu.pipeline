#' Finalize a sample table after it has been read or inferred
#'
#' This includes adding the `orient` column if it is not already present,
#' generating file names for the trimmed and filtered reads, and validation.
#'
#' @param sample_table (`tibble`) the sample table to finalize.
#' @return A finalized sample table
#' @keywords internal
finalize_sample_table <- function(
  sample_table
) {
  checkmate::assert_tibble(sample_table)
  checkmate::check_names(
    names(sample_table),
    must.include = c("sample", "seqrun", "fastq_R1", "fastq_R2",
                     "neg_control", "pos_control")
  )

  switch(
    read_orientation(),
    fwd = sample_table$orient <- "fwd",
    rev = sample_table$orient <- "rev",
    mixed =
      sample_table <- tidyr::crossing(sample_table, orient = c("fwd", "rev")),
    custom = if (isFALSE(do_custom_sample_table())) {
      stop("option 'orient: custom' requires a custom sample table is given.")
    } else if (!"orient" %in% names(sample_table)) {
      stop("option 'orient: custom' required a column named 'orient' in the",
           " custom sample table, with values consisting of 'fwd', 'rev', and",
           " 'mixed'")
    },
    stop("unknown value for option 'orient'; should be 'fwd', 'rev', 'mixed',",
         " or 'custom'")
  )

  if (do_rarefy()) {
    sample_table <- tidyr::crossing(sample_table, rarefy_meta(dots = FALSE)) |>
      dplyr::mutate(
        source_R1 = fastq_R1,
        source_R2 = fastq_R2,
        sample_key = ifelse(
          rarefy_text == "full",
          paste(seqrun, sample, sep = "_"),
          paste(seqrun, sample, rarefy_text, sep = "_")
        ),
        fastq_R1 = ifelse(
          rarefy_text == "full",
          fastq_R1,
          sprintf("%s/%s_R1.fastq.gz", rarefy_path(), sample_key)
        ),
        fastq_R2 = ifelse(
          rarefy_text == "full",
          fastq_R2,
          sprintf("%s/%s_R2.fastq.gz", rarefy_path(), sample_key)
        )
      )
  } else {
    sample_table <- dplyr::mutate(
      sample_table,
      sample_key = paste(seqrun, sample, sep = "_")
    )
  }

  sample_table <- sample_table |>
    # generate filenames for trimmed and filtered reads
    dplyr::mutate(
      trim_R1 = file.path(
        trim_path(),
        paste(sample_key, orient, "R1_trim.fastq.gz", sep = "_")
      ),
      trim_R2 = file.path(
        trim_path(),
        paste(sample_key, orient, "R2_trim.fastq.gz", sep = "_")
      ),
      filt_R1 = file.path(
        filt_path(),
        paste(sample_key, orient, "R1_filt.fastq.gz", sep = "_")
      ),
      filt_R2 = file.path(
        filt_path(),
        paste(sample_key, orient, "R2_filt.fastq.gz", sep = "_")
      ),
      sample_key = file_to_sample_key(filt_R1) # to be sure
    )

  # spike_strength is used along with the nonspike/spike ratio to convert from
  # read number to "weight"
  if (!("spike_weight") %in% names(sample_table))
    sample_table$spike_weight <- 1

  assertthat::assert_that(
    !any(is.na(sample_table$seqrun)),
    !any(is.na(sample_table$sample)),
    is.numeric(sample_table$spike_weight),
    !any(duplicated(sample_table[c("fastq_R1", "orient")])),
    !any(duplicated(sample_table[c("fastq_R2", "orient")])),
    !any(duplicated(sample_table$trim_R1)),
    !any(duplicated(sample_table$trim_R2)),
    !any(duplicated(sample_table$filt_R1)),
    !any(duplicated(sample_table$filt_R2))
  )

  options(
    optimotu.pipeline.sample_table_hash = targets:::hash_object(sample_table)
  )
  sample_table
}

#' Read sample table from file
#'
#' This function reads a custom sample table from CSV, TSV, or Excel files.
#'
#' Required columns in the sample table are `sample`, `seqrun`, `fastq_R1`, and
#' `fastq_R2`. These are all `character` vectors, giving the sample name,
#' sequencing run identifier, and paths to the two read files, respectively.
#' Optional columns are `neg_control` and `pos_control`, both of which should be
#' `logical` values (or something which may be interpreted as such) indicating
#' whether the sample is a negative or positive control, respectively; as well
#' as `orient`, which may by `fwd`, `rev`, or `mixed` to indicate the
#' orientation of the reads in each sample. Additionally trimming and filtering
#' parameters may be specified in additional columns, with names matching the
#' parameter names in `pipeline_options.yaml`.
#'
#' @param sample_table_file (`character`) the path to the sample table file.
#'
#' @return A `tibble` containing the sample table.
#' @keywords internal
read_sample_table <- function(sample_table_file = custom_sample_table()) {
  if (endsWith(sample_table_file, ".csv")) {
    sample_table <-
      suppressWarnings(
        readr::read_csv(
          sample_table_file,
          col_types = readr::cols(
            seqrun = readr::col_character(),
            sample = readr::col_character(),
            neg_control = readr::col_character(),
            pos_control = readr::col_character(),
            fastq_R1 = readr::col_character(),
            fastq_R2 = readr::col_character(),
            orient = readr::col_character(),
            .default = readr::col_character()
          )
        )
      )
  } else if (endsWith(sample_table_file, ".tsv")) {
    sample_table <-
      suppressWarnings(
        readr::read_tsv(
          sample_table_file,
          col_types = readr::cols(
            seqrun = readr::col_character(),
            sample = readr::col_character(),
            neg_control = readr::col_character(),
            pos_control = readr::col_character(),
            fastq_R1 = readr::col_character(),
            fastq_R2 = readr::col_character(),
            orient = readr::col_character(),
            .default = readr::col_character()
          )
        )
      )
  } else if (
    endsWith(sample_table_file, ".xls") ||
    endsWith(sample_table_file, ".xlsx")
  ) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("The 'readxl' package is required to read Excel files. Please",
           " install it.")
    }
    sample_table <-
      suppressWarnings(
        readxl::read_excel(
          sample_table_file,
          col_names = TRUE,
          col_types ="text"
        )
      )
  } else {
    stop("Unsupported sample table format. Please provide a .csv, .tsv, .xls,",
         " or .xlsx file.")
  }
  sample_table$fastq_R1 <- file.path(raw_path(), sample_table$fastq_R1)
  sample_table$fastq_R2 <- file.path(raw_path(), sample_table$fastq_R2)

  checkmate::check_data_frame(
    sample_table,
    col.names = "named"
  )
  checkmate::check_names(
    names(sample_table),
    must.include = c("sample", "seqrun", "fastq_R1", "fastq_R2")
  )
  checkmate::assert_character(sample_table$sample, any.missing = FALSE)
  checkmate::assert_character(sample_table$seqrun, any.missing = FALSE)
  checkmate::assert_file_exists(sample_table$fastq_R1, access = "r")
  checkmate::assert_file_exists(sample_table$fastq_R2, access = "r")
  if ("pos_control" %in% names(sample_table)) {
    checkmate::assert_subset(sample_table$pos_control, c(true_vals, false_vals))
    sample_table$pos_control <- sample_table$pos_control %in% true_vals
  } else {
    sample_table$pos_control <- FALSE
  }
  if ("neg_control" %in% names(sample_table)) {
    checkmate::assert_subset(sample_table$neg_control, c(true_vals, false_vals))
    sample_table$neg_control <- sample_table$neg_control %in% true_vals
  } else {
    sample_table$neg_control <- FALSE
  }
  if ("orient" %in% names(sample_table)) {
    checkmate::assert_subset(sample_table$orient, c("fwd", "rev", "mixed"))
    if (read_orientation() != "custom") {
      warning(
        "custom sample table '", sample_table_file,
        "' includes an 'orient' column, but option 'orient' is '",
        read_orientation(), "'. The 'orient' column will be ignored."
      )
      sample_table$orient <- NULL
    } else {
      if (any(sample_table$orient == "mixed")) {
        sample_table <-
          dplyr::left_join(
            sample_table,
            tibble::tibble(
              orient = c("fwd", "rev", "mixed", "mixed"),
              new_orient = c("fwd", "rev", "fwd", "rev")
            ),
            by = "orient",
            multiple = "all"
          ) |>
          dplyr::mutate(orient = new_orient, .keep = "unused")
      }
    }
  }
  sample_table <- dplyr::mutate(
    sample_table,
    dplyr::across(
      any_of(c("truncQ_R1", "truncQ_R2", "cut_R1", "cut_R2")),
      \(x) lapply(strsplit(as.character(x), ","), as.numeric)
    )
  )
  sample_table
}


#' Infer sample table from raw sequence file names
#'
#' This function infers a sample table from the directory structure and file
#' names of paired FASTQ files.
#'
#' This function is designed for the case where:
#'
#' 1) Each sequencing run is in a separate subdirectory of `raw_path`. The
#'   actual read files can be nested in further subdirectories below this.
#' 2) The read files are named according to the Illumina convention, with
#'   `R1` and `R2` suffixes for the forward and reverse reads, respectively.
#' 3) The first part of the file name is the sample name.
#' 4) Between the sample name and the `R1`/`R2` suffix, there may be sample and
#'   lane indices of the form "S{nnn}_L{nnn}".  These are ignored.
#' 5) Sample name, sample/lane indices, and `R1`/`R2` suffixes may be separated
#'   by underscores or dots.
#' 6) After the `R1`/`R2` suffix, there may be an additional `_nnn` suffix,
#'   which is ignored.
#' 7) The file extension should match the regular expression in
#'   `file_extension`. The default value matches `.fastq`, `.fq`, `.fastq.gz`,
#'   and `.fq.gz`, which should cover the majority of cases.
#' 8) Sample names including the text "NEGPCR", "NEGEXT", "NEGATIVE",
#'   "NEGCONTROL", "NEG_CONTROL", or "BLANK" (matched case-insensitively) are
#'   assumed to be negative controls, and those including "MOCK", "AMPTK",
#'   "POSITIVE", "POSCONTROL", or "POS_CONTROL" are assumed to be positive
#'   controls.
#'
#' @param raw_path (`character`) the path to the directory containing the raw
#' sequence files.
#' @param file_extension (`character`) a regular expression that matches the
#' file extension of the sequence files.
#'
#' @return A `tibble` with columns `fastq_R1`, `fastq_R2`, `seqrun`, `sample`,
#' `neg_control`, and `pos_control`. The `fastq_R1` and `fastq_R2` columns
#' contain the paths to the forward and reverse read files, respectively.
#' The `seqrun` column contains the name of the sequencing run (the directory
#' containing the read files), and the `sample` column contains the sample
#' name extracted from the file names. The `neg_control` and `pos_control`
#' columns are logical vectors indicating whether the sample is a negative or
#' positive control, respectively.
#' @keywords internal
infer_sample_table <- function(
    raw_path = raw_path(),
    file_extension = read_file_extension()
) {
  tibble::tibble(
    fastq_R1 = sort(list.files(
      raw_path,
      paste0(".*R1(_\\d{3})?[.]", file_extension),
      recursive = TRUE
    )),
    fastq_R2 = sort(list.files(
      raw_path,
      paste0(".*R2(_\\d{3})?[.]", file_extension),
      recursive = TRUE
    ))
  ) |>
    # parse filenames
    tidyr::extract(
      fastq_R1,
      into = c("seqrun", "sample"),
      regex = paste0(
        "([^/]+)/(?:.*/)?(.+?)[._](?:S\\d+_L00\\d[._])?R1(?:_001)?[.]",
        file_extension
      ),
      remove = FALSE
    ) |>
    dplyr::mutate(
      fastq_R1 = file.path(raw_path, fastq_R1),
      fastq_R2 = file.path(raw_path, fastq_R2),
      neg_control = grepl(
        "NEGPCR|NEGEXT|NEGATIVE|NEG_?CONTROL|BLANK",
        sample,
        ignore.case = TRUE
      ),
      pos_control = grepl(
        "MOCK|AMPTK|POSITIVE|POS_?CONTROL",
        sample,
        ignore.case = TRUE
      )
    )
}

.optimotu_cache <- new.env()

#' Return the sample table
#'
#' This function reads the custom sample table from disk or infers it from the
#' raw sequence files, depending on the configuration in
#' `pipeline_options.yaml`. It caches the result to avoid reading the sample
#' table multiple times in the same R session. It is best to avoid calling this
#' function on remote workers.
#'
#' @param ... arguments to be used for `targets` dependency tracking. Ignored
#' by this function.
#' @return A `tibble` with the sample table, as described in
#' `read_sample_table()`
#' @export
sample_table <- function(...) {
  if (is.null(.optimotu_cache$sample_table)) {
    if (do_custom_sample_table()) {
      .optimotu_cache$sample_table <- read_sample_table(
        custom_sample_table()
      )
    } else {
      .optimotu_cache$sample_table <- infer_sample_table(
        raw_path = raw_path(),
        file_extension = read_file_extension()
      )
    }
    .optimotu_cache$sample_table <-
      finalize_sample_table(.optimotu_cache$sample_table)
  }
  .optimotu_cache$sample_table
}

#' @rdname sample_table
#' @export
sample_table_hash <- function() {
  targets:::hash_object(sample_table())
}

#' @rdname paths
#' @export
path <- function() {
  "."
}

#' @rdname paths
#' @export
meta_path <- function() {
  fp <- file.path(path(), "meta")
  if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
  fp
}

#' @rdname paths
#' @export
seq_path <- function() {
  fp <- file.path(path(), "sequences")
  if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
  fp
}

#' @rdname paths
#' @export
raw_path <- function() {
  fp <- file.path(seq_path(), "01_raw")
  if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
  fp
}

#' @rdname paths
#' @export
rarefy_path <- function() {
  fp <- file.path(seq_path(), "02_rarefied")
  if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
  fp
}

#' @rdname paths
#' @export
trim_path <- function() {
  fp <- file.path(seq_path(), "03_trimmed")
  if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
  fp
}

#' @rdname paths
#' @export
filt_path <- function() {
  fp <- file.path(seq_path(), "04_filtered")
  if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
  fp
}

#' @rdname paths
#' @export
asv_path <- function() {
  fp <- file.path(seq_path(), "05_denoised")
  if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
  fp
}

#' @rdname paths
#' @export
aligned_path <- function() {
  fp <- file.path(seq_path(), "06_aligned")
  if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
  fp
}

#' @rdname paths
#' @export
protax_path <- function() {
  fp <- file.path(seq_path(), "07_protax")
  if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
  fp
}

#' @rdname paths
#' @export
output_path <- function() {
    fp <- file.path(path(), "output")
  if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
  fp
}

#' @rdname paths
#' @export
log_path <- function() {
  fp <- file.path(path(), "logs")
  if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
  fp
}
