#' Parse taxonomy from tabular (TSV) lines
#' @param lines (`character`) lines of a TSV file (with or without header)
#' @param ranks (`character` vector) the taxonomic ranks to parse
#' from most inclusive (e.g., kingdom) to least inclusive (e.g., species)
#' @param rank_map (two-column `data.frame`) a mapping of taxonomic ranks
#' as used in the project (first column) to their names in the taxonomy file
#' (second column).  If not provided, the ranks are assumed to be the same.
#' @noRd
parse_taxonomy_tabular <- function(lines, ranks = tax_ranks(), rank_map = NULL) {
  if (length(lines) == 0) {
    return(tibble::tibble(seq_id = character(), !!!stats::setNames(
      rep(list(character()), length(ranks)), ranks
    )))
  }
  n_ranks <- length(ranks)
  text <- paste(lines, collapse = "\n")

  # Try reading with header
  tab <- readr::read_tsv(
    I(text),
    show_col_types = FALSE,
    na = character()
  )
  tab <- dplyr::rename(tab, seq_id = 1)

  # If rank_map provided, rename file rank columns to project rank names
  if (!is.null(rank_map)) {
    checkmate::assert_data_frame(rank_map, ncols = 2, min.rows = 1)
    file_ranks <- rank_map[[2]]
    project_ranks <- rank_map[[1]]
    for (i in seq_len(nrow(rank_map))) {
      if (file_ranks[[i]] %in% names(tab)) {
        names(tab)[names(tab) == file_ranks[[i]]] <- project_ranks[[i]]
      }
    }
  }

  has_compliant_header <-
    all(ranks %in% names(tab)) || "taxonomy" %in% names(tab)

  if (has_compliant_header) {
    return(parse_taxonomy_tabular_with_header(tab, ranks))
  }

  # No compliant header: re-read as headerless and detect format
  tab0 <- readr::read_tsv(
    I(text),
    col_names = FALSE,
    show_col_types = FALSE,
    na = character()
  )
  n_cols <- ncol(tab0)

  if (n_cols == 2L) {
    return(parse_taxonomy_tabular_headerless_taxonomy(tab0, ranks))
  }
  if (n_cols == 1L + n_ranks) {
    return(parse_taxonomy_tabular_headerless_ranks(tab0, ranks))
  }

  stop(
    "Could not detect tabular format. With no header, expected either ",
    "2 columns (seq_id, taxonomy) or ", 1L + n_ranks,
    " columns (seq_id + ", n_ranks, " rank columns in order); got ", n_cols,
    " columns."
  )
}

parse_taxonomy_tabular_with_header <- function(tab, ranks) {
  # avoid R CMD check NOTE about global variables due to NSE
  seq_id <- taxonomy <- NULL
  rank_cols_present <- ranks %in% names(tab)
  if (all(rank_cols_present)) {
    out <- dplyr::select(tab, seq_id, dplyr::all_of(ranks))
    return(dplyr::mutate(out, dplyr::across(
      dplyr::any_of(ranks),
      \(x) dplyr::na_if(trimws(x), "")
    )))
  }
  if (!"taxonomy" %in% names(tab)) {
    stop(
      "Tabular file must have either a column for each requested rank (",
      paste(ranks, collapse = ", "),
      ") or a column named 'taxonomy'."
    )
  }
  parsed <- parse_taxonomy_header(tab[["taxonomy"]], ranks)
  if (!is.null(parsed)) {
    out <- dplyr::mutate(parsed, seq_id = tab[["seq_id"]])
    return(dplyr::select(out, seq_id, dplyr::any_of(ranks)))
  }
  tidyr::separate_wider_delim(
    tab,
    taxonomy,
    delim = stringr::regex("[,;|]"),
    too_few = "align_start",
    too_many = "drop",
    names = ranks
  ) |>
  dplyr::mutate(
    dplyr::across(
      dplyr::any_of(ranks),
      \(x) dplyr::na_if(trimws(x), "")
    )
  )
}

parse_taxonomy_tabular_headerless_ranks <- function(tab0, ranks) {
  # avoid R CMD check NOTE about global variables due to NSE
  seq_id <- NULL
  names(tab0) <- c("seq_id", ranks)
  out <- dplyr::mutate(tab0, dplyr::across(
    dplyr::any_of(ranks),
    \(x) dplyr::na_if(trimws(as.character(x)), "")
  ))
  dplyr::select(out, seq_id, dplyr::all_of(ranks))
}

parse_taxonomy_tabular_headerless_taxonomy <- function(tab0, ranks) {
  seq_id <- as.character(tab0[[1]])
  taxonomy <- as.character(tab0[[2]])
  tax_parts <- strsplit(taxonomy, "[,;|]")
  n_ranks <- length(ranks)
  out <- matrix(NA_character_, nrow = length(seq_id), ncol = n_ranks)
  for (i in seq_along(tax_parts)) {
    vals <- trimws(tax_parts[[i]])
    n <- min(length(vals), n_ranks)
    if (n > 0) {
      out[i, seq_len(n)] <- vals[seq_len(n)]
    }
  }
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  names(out) <- ranks
  out <- dplyr::mutate(out, dplyr::across(
    dplyr::everything(),
    \(x) dplyr::na_if(x, "")
  ))
  tibble::as_tibble(dplyr::bind_cols(seq_id = seq_id, out))
}

#' Parse taxonomy from a fasta, fastq, or tabular file
#'
#' @description Accepts various formats of taxonomy for fasta/fastq headers, including:
#' - `"sintax"`: `(ID);(other_fields;)*tax=d:(kingdom),p:(phylum),c:(class),o:(order),f:(family),g:(genus),s:(species)`
#' - `"BOLD"`: `(ID)|(marker_name)|(country)|(kingdom),(phylum),(class),(order),(family),(subfamily),(tribe),(genus),(species),(subspecies)`
#' - `"UNITE"`: `(ID)|k__(kingdom);p__(phylum);c__(class);o__(order);f__(family);g__(genus);s__(species)|(SH_ID)`
#'
#' For tabular files, the first column is assumed to be the sequence ID.
#' There should also be a named column for each taxonomic rank, or
#' alternatively, a single column named "taxonomy" that is either a
#' comma-, semicolon-, or pipe-delimited string of taxonomic ranks, or matches
#' one of the above formats for fasta headers.
#' @param file (`character`) path to the file to parse
#' @param ranks (`character` vector) the taxonomic ranks to parse, in order
#' from most inclusive (e.g., kingdom) to least inclusive (e.g., species)
#' @param rank_map (two-column `data.frame`) a mapping of taxonomic ranks
#' as used in the project (first column) to their names in the taxonomy file
#' (second column).  If not provided, the ranks are assumed to be the same.
#' @return a `data.frame` with column `seq_id` and named columns for each
#' taxonomic rank in `ranks`
#' @export
parse_reference_taxonomy <- function(file, ranks = tax_ranks(), rank_map = NULL) {
  checkmate::assert_file(file, "r")
  if (grepl("\\.fa(s(ta)?)?", file)) {
    result <- parse_taxonomy_header(
      names(Biostrings::fasta.seqlengths(file)), ranks, rank_map
    )
    if (is.null(result)) {
      stop("Unrecognized sequence header format. Supported: sintax, BOLD, UNITE.")
    }
    result
  } else if (grepl("\\.f(a(st)?)?q", file)) {
    result <- parse_taxonomy_header(fastq_names(file), ranks, rank_map)
    if (is.null(result)) {
      stop("Unrecognized sequence header format. Supported: sintax, BOLD, UNITE.")
    }
    result
  } else {
    parse_taxonomy_tabular(readLines(file), ranks, rank_map)
  }
}

is_sintax_header <- function(header) {
  sintax_regex <-
    "^([^;]+;)+tax=([dkpcofgst]:[^,;]*,)*[dkpcofgst]:[^,;]*(;.*)?$"
  if (length(header) <= 100 || all(grepl(sintax_regex, header[1:100]))) {
    all(grepl(sintax_regex, header))
  } else {
    FALSE
  }
}

sintax_ranks <- c(d = "domain", k = "kingdom", p = "phylum",
                  c = "class", o = "order", f = "family",
                  g = "genus", s = "species", t = "strain")

parse_sintax_header <- function(header) {
  # avoid R CMD check NOTE about global variables due to NSE
  seq_id <- taxonomy <- NULL
  tibble::tibble(header = header) |>
    tidyr::separate_wider_regex(
      header,
      patterns = c(
        seq_id = "[^;]+",
        ";(?:.*;)?tax=",
        taxonomy = "(?:[dkpcofgst]:[^,;]*,)*[dkpcofgst]:[^,;]*",
        "(?:;.*)?"
      )
    ) |>
    tidyr::separate_rows(
      taxonomy,
      sep = ","
    ) |>
    dplyr::mutate(
      rank = sintax_ranks[substr(taxonomy, 1, 1)],
      taxonomy = sub("^[dkpcofgst]:", "", taxonomy)
    ) |>
    tidyr::pivot_wider(
      id_cols = seq_id,
      names_from = rank,
      values_from = taxonomy
    ) |>
    dplyr::select(seq_id, dplyr::any_of(unname(sintax_ranks)))
}

is_bold_header <- function(header) {
  bold_regex <- "^([^|,;]+\\|){2,3}([^,;|]*,)*[^,;|]*$"
  if (length(header) <= 100 || all(grepl(bold_regex, header[1:100]))) {
    all(grepl(bold_regex, header))
  } else {
    FALSE
  }
}

is_unite_header <- function(header) {
  unite_regex <- "^[^|]+\\|[kpcofgs]__[^;,]+(;[kpcofgs]__[^;,]+)*\\|[^|]+$"
  if (length(header) <= 100 || all(grepl(unite_regex, header[1:100]))) {
    all(grepl(unite_regex, header))
  } else {
    FALSE
  }
}

bold_ranks <- c(
  "kingdom", "phylum", "class", "order", "family",
  "subfamily", "tribe", "genus", "species", "subspecies"
)

parse_bold_header <- function(header) {
  # avoid R CMD check NOTE about global variables due to NSE
  seq_id <- taxonomy <- NULL
  parts <- strsplit(header, "|", fixed = TRUE)
  seq_id <- vapply(parts, `[`, 1L, FUN.VALUE = character(1))
  taxonomy_str <- vapply(parts, utils::tail, 1L, FUN.VALUE = character(1))
  tibble::tibble(seq_id = seq_id, taxonomy = taxonomy_str) |>
    tidyr::separate_wider_delim(
      taxonomy,
      delim = ",",
      names = bold_ranks,
      too_many = "drop",
      too_few = "align_start"
    ) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::any_of(bold_ranks),
        \(x) dplyr::na_if(trimws(x), "")
      )
    ) |>
    dplyr::select(seq_id, dplyr::any_of(bold_ranks))
}

unite_ranks <- c(
  k = "kingdom", p = "phylum", c = "class", o = "order", f = "family",
  g = "genus", s = "species"
)

parse_unite_header <- function(header) {
  # avoid R CMD check NOTE about global variables due to NSE
  taxonomy <- NULL
  parts <- strsplit(header, "|", fixed = TRUE)
  seq_id <- vapply(parts, `[`, 1L, FUN.VALUE = character(1))
  taxonomy_str <- vapply(parts, `[`, 2L, FUN.VALUE = character(1))
  tibble::tibble(seq_id = seq_id, taxonomy = taxonomy_str) |>
    tidyr::separate_rows(
      taxonomy,
      sep = ";"
    ) |>
    dplyr::mutate(
      rank = unite_ranks[substr(taxonomy, 1, 1)],
      taxonomy = dplyr::na_if(sub("^[kpcofgs]__", "", taxonomy), "")
    ) |>
    tidyr::pivot_wider(
      id_cols = seq_id,
      names_from = rank,
      values_from = taxonomy
    ) |>
    dplyr::select(seq_id, dplyr::any_of(unname(unite_ranks)))
}

#' Parse taxonomy from sequence headers (sintax, BOLD, or UNITE format)
#' @param header (`character`) vector of sequence header strings
#' @inheritParams parse_taxonomy
#' @return A tibble with `seq_id` and rank columns, or `NULL` if format unrecognized
#' @keywords internal
#' @noRd
parse_taxonomy_header <- function(header, ranks = tax_ranks(), rank_map = NULL) {
  # avoid R CMD check NOTE about global variables due to NSE
  seq_id <- NULL
  out <- if (is_sintax_header(header)) {
    parse_sintax_header(header)
  } else if (is_bold_header(header)) {
    parse_bold_header(header)
  } else if (is_unite_header(header)) {
    parse_unite_header(header)
  } else {
    return(NULL)
  }
  if (!is.null(rank_map)) {
    checkmate::assert_data_frame(rank_map, ncols = 2, min.rows = 1)
    checkmate::assert_character(rank_map[[1]], any.missing = FALSE)
    checkmate::assert_character(rank_map[[2]], any.missing = FALSE)
    checkmate::assert_subset(rank_map[[1]], ranks)
    ranks_from <- rank_map[[2]][rank_map[[2]] %in% names(out)]
    ranks_to <- rank_map[[1]][rank_map[[2]] %in% names(out)]
    out <- out[c("seq_id", ranks_from)]
    names(out) <- c("seq_id", ranks_to)
  }
  dplyr::select(out, seq_id, dplyr::any_of(ranks))
}
