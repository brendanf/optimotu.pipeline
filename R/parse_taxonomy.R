#' Parse taxonomy from tabular (TSV) file
#' @param file (`character`) path to a TSV file (with or without header)
#' @param ranks (`character` vector) the taxonomic ranks to parse
#' from most inclusive (e.g., kingdom) to least inclusive (e.g., species)
#' @param rank_map (two-column `data.frame`) a mapping of taxonomic ranks
#' as used in the project (first column) to their names in the taxonomy file
#' (second column).  If not provided, the ranks are assumed to be the same.
#' @noRd
parse_taxonomy_tabular <- function(
  file,
  ranks = tax_ranks(),
  rank_map = NULL
) {
  checkmate::assert(
    checkmate::check_file_exists(file, "r"),
    checkmate::check_class(file, "AsIs")
  )
  n_ranks <- length(ranks)

  file_ranks <- ranks
  if (!is.null(rank_map)) {
    checkmate::assert_data_frame(rank_map, min.cols = 2)
    checkmate::assert_subset(ranks, rank_map[[1]])
    file_ranks <- rank_map[[2]][match(ranks, rank_map[[1]])]
  }

  # Read first line to determine if file has header
  first_line <- (if (methods::is(file, "AsIs")) {
    checkmate::assert_character(file)
    if (length(file) == 0) {
      return(
        tibble::tibble(
          seq_id = character(),
          !!!stats::setNames(
            replicate(n_ranks, character(), simplify = FALSE),
            ranks
          )
        )
      )
    }
    first_line <- strsplit(file, "\n")[[1]][1]
  } else {
    first_line <- readLines(file, n = 1)
  }) |>
    strsplit("\t") |>
    unlist()
  n_cols <- length(first_line)
  has_compliant_header <-
    all(file_ranks %in% utils::tail(first_line, -1)) ||
    "taxonomy" %in% utils::tail(first_line, -1)

  if (has_compliant_header) {
    out <- parse_taxonomy_tabular_with_header(file, file_ranks)
  } else if (n_cols == 2L) {
    out <- parse_taxonomy_tabular_headerless_taxonomy(file, file_ranks)
  } else if (n_cols == 1L + n_ranks) {
    out <- parse_taxonomy_tabular_headerless_ranks(file, file_ranks)
  } else if (n_cols == 3L) {
    out <- parse_taxonomy_tabular_headerless_taxonomy(
      file,
      file_ranks,
      middle_column = TRUE
    )
  } else {
    stop(
      "Could not detect tabular format. With no header, expected either ",
      "2 columns (seq_id, taxonomy), 3 columns (seq_id, [ignored], taxonomy), ",
      "or 1 + ",
      n_ranks,
      " columns (seq_id + ",
      n_ranks,
      " rank columns in ",
      "order); got ",
      n_cols,
      " columns."
    )
  }
  if (!is.null(rank_map)) {
    names(out)[-1] <- rank_map[[1]][match(names(out)[-1], rank_map[[2]])]
  }
  dplyr::mutate(
    out,
    dplyr::across(
      dplyr::any_of(ranks),
      \(x) dplyr::na_if(trimws(x), "")
    )
  )
}

parse_taxonomy_tabular_with_header <- function(file, ranks) {
  # avoid R CMD check NOTE about global variables due to NSE
  seq_id <- taxonomy <- NULL
  tab <- readr::read_tsv(file, col_types = c(.default = "c")) |>
    dplyr::rename(seq_id = 1)

  if (all(ranks %in% names(tab)[-1])) {
    dplyr::select(tab, seq_id, dplyr::all_of(ranks)) |>
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(ranks),
          \(x) dplyr::na_if(trimws(x), "")
        )
      )
  } else if ("taxonomy" %in% names(tab)) {
    parsed <- parse_taxonomy_header(tab[["taxonomy"]], ranks)
    if (is.null(parsed)) {
      tidyr::separate_wider_delim(
        tab,
        cols = taxonomy,
        delim = stringr::regex("[,;|]"),
        names = ranks,
        too_few = "align_start",
        too_many = "drop"
      )
    } else {
      dplyr::bind_cols(
        dplyr::select(tab, seq_id),
        dplyr::select(parsed, any_of(ranks))
      )
    }
  } else {
    stop(
      "Tabular file must have either a column for each requested rank (",
      paste(ranks, collapse = ", "),
      ") or a column named 'taxonomy'."
    )
  }
}

parse_taxonomy_tabular_headerless_ranks <- function(file, ranks) {
  readr::read_tsv(
    file,
    col_types = c(.default = "c"),
    col_names = c("seq_id", ranks)
  )
}

parse_taxonomy_tabular_headerless_taxonomy <- function(
  file,
  ranks,
  middle_column = FALSE
) {
  checkmate::assert_flag(middle_column)

  tab <- readr::read_tsv(
    file,
    col_types = if (middle_column) "c-c" else "cc",
    col_names = c("seq_id", "taxonomy")
  )
  taxonomy <- parse_taxonomy_header(tab[["taxonomy"]], ranks)
  if (is.null(taxonomy)) {
    tidyr::separate_wider_delim(
      tab,
      cols = taxonomy,
      delim = stringr::regex("[,;|]"),
      names = ranks,
      too_few = "align_start",
      too_many = "drop"
    )
  } else {
    dplyr::bind_cols(
      dplyr::select(tab, "seq_id"),
      dplyr::select(taxonomy, any_of(ranks))
    )
  }
}

#' Parse taxonomy from a fasta, fastq, or tabular file
#'
#' @details Accepts various formats of taxonomy for fasta/fastq headers, including:
#' - `"sintax"`: `(ID);(other_fields;)*tax=d:(kingdom),p:(phylum),c:(class),o:(order),f:(family),g:(genus),s:(species)`
#' - `"BOLD"`: `(ID)|(marker_name)|(country)|(kingdom),(phylum),(class),(order),(family),(subfamily),(tribe),(genus),(species),(subspecies)`
#' - `"UNITE"`: `(ID)|k__(kingdom);p__(phylum);c__(class);o__(order);f__(family);g__(genus);s__(species)|(SH_ID)`
#' - `"BayesANT"`: `(ID) Root;(rank1);{...};(rankN)`
#'
#' For tabular files, the first column is assumed to be the sequence ID.
#' There should also be a named column for each taxonomic rank, or
#' alternatively, a single column named "taxonomy" that is either a
#' comma-, semicolon-, or pipe-delimited string of taxonomic ranks, or matches
#' one of the above formats for fasta headers. An additional column between the
#' sequence ID and taxonomic columns is allowed but ignored (for instance the
#' rank in a protax seqid2tax file).
#' @param file (`character`) path to the file to parse
#' @param ranks (`character` vector) the taxonomic ranks to parse, in order
#' from most inclusive (e.g., kingdom) to least inclusive (e.g., species)
#' @param rank_map (two-column `data.frame`) a mapping of taxonomic ranks
#' as used in the project (first column) to their names in the taxonomy file
#' (second column).  If not provided, the ranks are assumed to be the same.
#' @return a `data.frame` with column `seq_id` and named columns for each
#' taxonomic rank in `ranks`
#' @export
parse_reference_taxonomy <- function(
  file,
  ranks = tax_ranks(),
  rank_map = NULL
) {
  checkmate::assert_file(file, "r")
  if (grepl("\\.fa(s(ta)?)?(.gz)?$", file)) {
    result <- parse_taxonomy_header(
      names(Biostrings::fasta.seqlengths(file)),
      ranks,
      rank_map
    )
    if (is.null(result)) {
      stop(
        "Unrecognized sequence header format. Supported: SINTAX, BOLD, UNITE,",
        " and BayesANT."
      )
    }
    result
  } else if (grepl("\\.f(a(st)?)?q", file)) {
    result <- parse_taxonomy_header(fastq_names(file), ranks, rank_map)
    if (is.null(result)) {
      stop(
        "Unrecognized sequence header format. Supported: SINTAX, BOLD, UNITE,",
        " and BayesANT."
      )
    }
    result
  } else {
    parse_taxonomy_tabular(file, ranks, rank_map)
  }
}

#' Extract `seq_id` values from annotated FASTA headers
#'
#' Detects supported taxonomy header formats (SINTAX, BOLD, UNITE, BayesANT,
#' and Protax-style `id<tab>taxonomy`) and returns the sequence identifier
#' from each header, in input order. Use with
#' `optimotu::optimize_thresholds(seq_names = ...)` when FASTA headers
#' carry taxonomy and are not bare `seq_id` values.
#'
#' @param headers (`character`) FASTA header strings (without leading `>`)
#' @param allow_plain (`logical(1)`) if `TRUE`, unrecognized headers (for
#'   example bare IDs after `hmmalign`) are returned unchanged instead of
#'   erroring
#' @return (`character`) `seq_id` for each header
#' @export
fasta_header_seq_ids <- function(headers, allow_plain = FALSE) {
  checkmate::assert_character(headers, any.missing = FALSE)
  checkmate::assert_flag(allow_plain)
  if (length(headers) == 0L) {
    return(character())
  }
  # Same leading-field rules as parse_*_header(); avoid full taxonomy parse.
  if (is_sintax_header(headers)) {
    sub(";.*$", "", headers)
  } else if (is_bold_header(headers) || is_unite_header(headers)) {
    vapply(strsplit(headers, "|", fixed = TRUE), `[`, character(1), 1L)
  } else if (is_bayesant_header(headers)) {
    sub("[[:space:]].*$", "", headers)
  } else if (any(grepl("\t", headers, fixed = TRUE))) {
    # Aligned Protax refs.aln: short ID, tab, taxonomy string
    sub("\t.*$", "", headers)
  } else if (isTRUE(allow_plain)) {
    headers
  } else {
    stop(
      "Unrecognized sequence header format. Supported: SINTAX, BOLD, UNITE,",
      " BayesANT, and tab-separated Protax (`id<tab>taxonomy`)."
    )
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

sintax_ranks <- c(
  d = "domain",
  k = "kingdom",
  p = "phylum",
  c = "class",
  o = "order",
  f = "family",
  g = "genus",
  s = "species",
  t = "strain"
)

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

is_bayesant_header <- function(header) {
  bayesant_regex <- "^[^[:space:];]+[[:space:]]+Root(?:;[^;]*)+$"
  if (length(header) <= 100 || all(grepl(bayesant_regex, header[1:100]))) {
    all(grepl(bayesant_regex, header))
  } else {
    FALSE
  }
}

bold_ranks <- c(
  "kingdom",
  "phylum",
  "class",
  "order",
  "family",
  "subfamily",
  "tribe",
  "genus",
  "species",
  "subspecies"
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
  k = "kingdom",
  p = "phylum",
  c = "class",
  o = "order",
  f = "family",
  g = "genus",
  s = "species"
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

parse_bayesant_header <- function(header, ranks = tax_ranks()) {
  # avoid R CMD check NOTE about global variables due to NSE
  seq_id <- taxonomy <- pos <- NULL
  out <- tibble::tibble(header = header) |>
    tidyr::separate_wider_regex(
      header,
      patterns = c(
        seq_id = "[^[:space:];]+",
        "[[:space:]]+",
        taxonomy = ".*"
      )
    ) |>
    tidyr::separate_rows(
      taxonomy,
      sep = ";"
    ) |>
    dplyr::group_by(seq_id) |>
    dplyr::mutate(pos = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::filter(pos > 1) |>
    dplyr::mutate(rank = ranks[pos - 1]) |>
    dplyr::filter(!is.na(rank)) |>
    dplyr::mutate(taxonomy = dplyr::na_if(trimws(taxonomy), "")) |>
    tidyr::pivot_wider(
      id_cols = seq_id,
      names_from = rank,
      values_from = taxonomy
    )
  dplyr::select(out, seq_id, dplyr::any_of(ranks))
}

#' Parse taxonomy from sequence headers
#' @param header (`character`) vector of sequence header strings
#' @inheritParams parse_taxonomy
#' @return A tibble with `seq_id` and rank columns, or `NULL` if format
#' is unrecognized
#' @keywords internal
#' @noRd
parse_taxonomy_header <- function(
  header,
  ranks = tax_ranks(),
  rank_map = NULL
) {
  # avoid R CMD check NOTE about global variables due to NSE
  seq_id <- NULL
  out <- if (is_sintax_header(header)) {
    parse_sintax_header(header)
  } else if (is_bold_header(header)) {
    parse_bold_header(header)
  } else if (is_unite_header(header)) {
    parse_unite_header(header)
  } else if (is_bayesant_header(header)) {
    parse_bayesant_header(header, ranks)
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

#' Taxonomy table for one or more outgroup reference libraries
#'
#' Reads headers (or an optional sidecar) with
#' [parse_reference_taxonomy()] and prefixes `seq_id` as
#' `{library}:{seq_id}`. Sequence bodies are not read. Missing ranks are
#' `NA`. Duplicate ids within one library are an error. Unrecognized
#' headers name the library in the error.
#'
#' @param libraries (`list`) as returned by [outgroup_libraries()]
#' @param ranks (`character`) taxonomic ranks to keep
#' @return a `tibble` with `ref_id` and a column for each rank in `ranks`
#' @export
outgroup_library_taxonomy <- function(libraries, ranks = tax_ranks()) {
  checkmate::assert_list(libraries, min.len = 1L)
  checkmate::assert_character(ranks, min.len = 1L, any.missing = FALSE)
  rows <- lapply(libraries, function(lib) {
    checkmate::assert_list(lib)
    checkmate::assert_string(lib$name, min.chars = 1L)
    src <- lib$taxonomy
    if (is.null(src) || !nzchar(src)) {
      src <- lib$sequences
    }
    checkmate::assert_string(src, min.chars = 1L)
    parsed <- tryCatch(
      parse_reference_taxonomy(src, ranks = ranks),
      error = function(e) {
        stop(
          "Outgroup library '",
          lib$name,
          "': ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
    if (anyDuplicated(parsed$seq_id)) {
      stop(
        "Duplicate sequence ids in outgroup library '",
        lib$name,
        "'.",
        call. = FALSE
      )
    }
    parsed$ref_id <- paste(lib$name, parsed$seq_id, sep = ":")
    parsed$seq_id <- NULL
    missing_ranks <- setdiff(ranks, names(parsed))
    for (rank in missing_ranks) {
      parsed[[rank]] <- NA_character_
    }
    dplyr::select(parsed, "ref_id", dplyr::all_of(ranks))
  })
  dplyr::bind_rows(rows)
}

#' Keep the closest hit across outgroup libraries
#'
#' @param hits (`data.frame` or list of them) with `id_col`, `ref_id`,
#'   and `dist`
#' @param taxonomy (`data.frame`) from [outgroup_library_taxonomy()]
#' @param id_col (`character(1)`) query id column (`seq_idx` or `seq_id`)
#' @return `hits` reduced to one row per query, joined to `taxonomy`
#' @export
reduce_library_best_hits <- function(hits, taxonomy, id_col = "seq_idx") {
  # NSE columns used below
  dist <- NULL
  checkmate::assert_string(id_col, min.chars = 1L)
  if (is.list(hits) && !is.data.frame(hits)) {
    hits <- dplyr::bind_rows(hits)
  }
  checkmate::assert_data_frame(hits)
  checkmate::assert_names(
    names(hits),
    must.include = c(id_col, "ref_id", "dist")
  )
  checkmate::assert_data_frame(taxonomy)
  checkmate::assert_names(names(taxonomy), must.include = "ref_id")
  if (nrow(hits) == 0L) {
    return(dplyr::left_join(hits, taxonomy, by = "ref_id"))
  }
  dplyr::slice_min(hits, dist, by = dplyr::all_of(id_col), with_ties = FALSE) |>
    dplyr::left_join(taxonomy, by = "ref_id")
}

#' Parse streamed multi-library vsearch hit labels
#'
#' Headers are `{library}:{original header}`. Mixed BOLD and SINTAX
#' labels are split on the first `:` and parsed per library.
#'
#' @param headers (`character`) vsearch target labels
#' @return (`character`) `ref_id` values matching
#'   [outgroup_library_taxonomy()]
#' @export
outgroup_hit_ref_ids <- function(headers) {
  checkmate::assert_character(headers, any.missing = FALSE)
  if (length(headers) == 0L) {
    return(character())
  }
  if (any(!grepl(":", headers, fixed = TRUE))) {
    stop(
      "Expected library-prefixed hit labels of the form 'library:header'.",
      call. = FALSE
    )
  }
  prefix <- sub(":.*$", "", headers)
  out <- character(length(headers))
  split_i <- split(seq_along(headers), prefix)
  for (lib in names(split_i)) {
    idx <- split_i[[lib]]
    parsed <- tryCatch(
      fasta_header_seq_ids(headers[idx]),
      error = function(e) {
        stop(
          "Outgroup library '",
          lib,
          "': ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
    out[idx] <- parsed
  }
  out
}
