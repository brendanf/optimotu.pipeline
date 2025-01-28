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
            col_positions =  stringr::str_locate_all(x$text[1], "-+")[[1]] |>
              tibble::as_tibble() |>
              tibble::add_column(
                col_names = c("idx", "seq_id", "match_len", "cm_from", "cm_to",
                              "trunc", "bit_sc", "avg_pp", "time_band_calc",
                              "time_alignment", "time_total", "mem_mb")
              ) |>
              do.call(readr::fwf_positions, args = _),
            skip = 1,
            col_types = "iciiicdddddd"
          )
      }
    )
}

# TODO: implement for newer version of inferrnal StockholmMultipleAlignment
#       (or move to inferrnal)

#' Extract the consensus columns from a multiple sequence alignment
#' @param aln (`list`) a list with elements `alignment` (a `MultipleAlignment`
#' object) and `GC` (a `BString` object with the consensus sequence)
#' @return a `DNAStringSet` with the consensus columns
#' @export
consensus_columns <- function(aln) {
  checkmate::assert_names(names(aln), must.include = c("alignment", "GC"))
  checkmate::assert_names(names(aln$GC), must.include = "RF")
  checkmate::assert_class(aln$alignment, "MultipleAlignment")
  checkmate::assert_class(aln$GC$RF, "BString")
  dots <- gregexpr("[.]+", aln$GC$RF)[[1]]
  Biostrings::colmask(aln$alignment) <-
    IRanges::IRanges(start = dots, width = attr(dots, "match.length"))
  methods::as(aln$alignment, "DNAStringSet")
}
