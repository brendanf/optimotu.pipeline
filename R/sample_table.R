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
  } else if (endsWith(sample_table_file, ".xls") || endsWith(sample_table_file, ".xlsx")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("The 'readxl' package is required to read Excel files. Please install it.")
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
    stop("Unsupported sample table format. Please provide a .csv, .tsv, .xls, or .xlsx file.")
  }

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
  checkmate::assert_file_exists(file.path(raw_path, sample_table$fastq_R1), access = "r")
  checkmate::assert_file_exists(file.path(raw_path, sample_table$fastq_R2), access = "r")
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
#' 7) The file extension should match the regular expression in `file_extension`.
#'   The default value matches `.fastq`, `.fq`, `.fastq.gz`, and `.fq.gz`, which
#'   should cover the majority of cases.
#' 8) Sample named including the test "NEGPCR", "NEGEXT", "NEGATIVE",
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
#' @export
infer_sample_table <- function(
    raw_path = "sequences/01_raw",
    file_extension = read_file_extension()
) {
  tibble::tibble(
    fastq_R1 = sort(list.files(raw_path, paste0(".*R1(_\\d{3})?[.]", file_extension), recursive = TRUE)),
    fastq_R2 = sort(list.files(raw_path, paste0(".*R2(_\\d{3})?[.]", file_extension), recursive = TRUE))
  ) |>
    # parse filenames
    tidyr::extract(
      fastq_R1,
      into = c("seqrun", "sample"),
      regex = paste0("([^/]+)/(?:.*/)?(.+?)[._](?:S\\d+_L00\\d[._])?R1(?:_001)?[.]", file_extension),
      remove = FALSE
    ) |>
    dplyr::mutate(
      neg_control = grepl("NEGPCR|NEGEXT|NEGATIVE|NEG_?CONTROL|BLANK", sample, ignore.case = TRUE),
      pos_control = grepl("MOCK|AMPTK|POSITIVE|POS_?CONTROL", sample, ignore.case = TRUE)
    )
}

