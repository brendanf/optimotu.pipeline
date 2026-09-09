#' Parse the "sfile" written by infernal's cmalign
#' @param file (`character`) the path to the file
#' @return a `data.frame` with columns `idx`, `seq_id`, `match_len`, `cm_from`,
#' `cm_to`, `trunc`, `bit_sc`, `avg_pp`, `time_band_calc`, `time_alignment`,
#' `time_total`, and `mem_mb`
#' @export
read_sfile <- function(file) {
  # avoid R CMD check NOTE: no visible binding for global variable
  text <- is_widths <- part <- NULL

  tibble::tibble(
    text = readLines(file),
    is_widths = grepl("^#[- ]+$", text),
    part = cumsum(is_widths)
  ) |>
    dplyr::filter(is_widths | !startsWith(text, "#")) |>
    dplyr::group_split(part, .keep = FALSE) |>
    purrr::discard(\(x) nrow(x) == 1) |>
    purrr::map_dfr(
      \(x) {
        paste(x$text, collapse = "\n") |>
          readr::read_fwf(
            col_positions = stringr::str_locate_all(x$text[1], "-+")[[1]] |>
              tibble::as_tibble() |>
              tibble::add_column(
                col_names = c(
                  "idx",
                  "seq_id",
                  "match_len",
                  "cm_from",
                  "cm_to",
                  "trunc",
                  "bit_sc",
                  "avg_pp",
                  "time_band_calc",
                  "time_alignment",
                  "time_total",
                  "mem_mb"
                )
              ) |>
              do.call(readr::fwf_positions, args = _),
            skip = 1,
            col_types = "iciiicdddddd"
          )
      }
    )
}

#' Extract the consensus columns from a multiple sequence alignment
#'
#' Columns marked with `.` in the Infernal reference (`RF`) annotation are
#' treated as inserts and removed via [Biostrings::colmask()].
#'
#' @param aln either:
#'   * an `inferrnal` `StockholmMultipleAlignment` (or subclass such as
#'     `StockholmDNAMultipleAlignment`), with a `GC` annotation named `"RF"`; or
#'   * a `list` with elements `alignment` (a
#'     [`MultipleAlignment`][Biostrings::MultipleAlignment-class]) and `GC`
#'     containing `RF` as a [`BString`][Biostrings::XString-class] (legacy
#'     `inferrnal` return shape)
#' @return a `DNAStringSet` with the consensus columns. RNA alignments are
#'   converted to DNA (`U` → `T`).
#' @export
consensus_columns <- function(aln) {
  if (methods::is(aln, "StockholmMultipleAlignment")) {
    checkmate::assert_names(names(aln@GC), must.include = "RF")
    rf <- aln@GC[["RF"]]
    # colmask<- on Stockholm* rebuilds S4 slots incorrectly; drop annotations.
    alignment <- methods::as(aln, "MultipleAlignment")
  } else {
    checkmate::assert_list(aln)
    checkmate::assert_names(names(aln), must.include = c("alignment", "GC"))
    checkmate::assert_names(names(aln$GC), must.include = "RF")
    checkmate::assert_class(aln$alignment, "MultipleAlignment")
    rf <- aln$GC$RF
    alignment <- aln$alignment
  }
  checkmate::assert_class(rf, "BString")
  dots <- gregexpr("[.]+", as.character(rf))[[1]]
  Biostrings::colmask(alignment) <-
    IRanges::IRanges(start = dots, width = attr(dots, "match.length"))
  xss <- methods::as(alignment, "XStringSet")
  if (methods::is(xss, "RNAStringSet")) {
    Biostrings::DNAStringSet(xss)
  } else {
    methods::as(xss, "DNAStringSet")
  }
}
