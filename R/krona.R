#' @noRd
krona_unknown_tip_ranks <- function(ranks = tax_ranks()) {
  utils::tail(ranks[-1L], 3L)
}

#' Krona node data columns for the last few unknown ranks
#'
#' Replaces the old hardcoded `sp` / `gen` / `fam` attributes. Keys are the
#' rank names themselves.
#'
#' @param ranks (`character`) taxonomic ranks, most inclusive first
#' @return named `list` suitable for [krona_xml_nodes()] `node_data_format`
#' @export
krona_node_data_format <- function(ranks = tax_ranks()) {
  checkmate::assert_character(ranks, min.len = 2L, any.missing = FALSE)
  tip3 <- krona_unknown_tip_ranks(ranks)
  out <- list(
    f = c("focc", "fread", "fotu"),
    nocc = rep("nocc", 3L),
    nread = rep("nread", 3L),
    notu = rep("notu", 3L)
  )
  for (r in tip3) {
    out[[r]] <- paste0(r, "_unknown_", c("focc", "fread", "fotu"))
  }
  out
}

#' Krona HTML attribute and color header
#'
#' XML lines for the `<attributes>` block plus the default color mapping
#' (tip-most unknown rank).
#'
#' @param ranks (`character`) taxonomic ranks, most inclusive first
#' @return (`character`) XML lines
#' @export
krona_html_attributes <- function(ranks = tax_ranks()) {
  checkmate::assert_character(ranks, min.len = 2L, any.missing = FALSE)
  tip3 <- krona_unknown_tip_ranks(ranks)
  c(
    '<attributes magnitude="f">',
    '<attribute display="Weighted fraction">f</attribute>',
    '<attribute display="Total occurences">nocc</attribute>',
    '<attribute display="Total reads">nread</attribute>',
    '<attribute display="Total OTUs">notu</attribute>',
    sprintf(
      paste0(
        '<attribute display="Weighted fraction belonging to unknown %s">',
        '%s</attribute>'
      ),
      tip3,
      tip3
    ),
    '</attributes>',
    sprintf(
      paste0(
        '<color attribute="%s" valueStart="0" valueEnd="1" ',
        'hueStart="120" hueEnd="0" default="true"></color>'
      ),
      tip3[[length(tip3)]]
    )
  )
}

krona_empty_table <- function(ranks) {
  unknown_ranks <- ranks[-1L]
  unknown_cols <- stats::setNames(
    rep(list(numeric()), length(unknown_ranks) * 3L),
    paste0(
      rep(unknown_ranks, each = 3L),
      "_unknown_",
      c("fread", "fotu", "focc")
    )
  )
  tibble::tibble(
    rank = rank2factor(character(), ranks),
    taxon = character(),
    parent_taxonomy = character(),
    !!!unknown_cols,
    nread = integer(),
    nocc = integer(),
    notu = integer(),
    child_unknown_fread = numeric(),
    child_unknown_fotu = numeric(),
    child_unknown_focc = numeric(),
    fread = numeric(),
    focc = numeric(),
    fotu = numeric()
  )
}

#' Generate data to write a KronaTools XML file
#'
#' Builds parent-taxonomy strings and per-rank unknown fractions for every
#' rank except the root. Taxa whose names start with `pseudo` count as
#' unknown. [remove_mycobank_number()] is applied to `genus` and `species`
#' when those ranks are present.
#'
#' @param otu_taxonomy (`data.frame`) OTU taxonomy with rank columns named in
#'   `ranks`, plus `nread` and `nsample`
#' @param ranks (`character`) taxonomic ranks, most inclusive first. Defaults
#'   to [tax_ranks()].
#' @return `data.frame` with columns `rank`, `taxon`, `parent_taxonomy`,
#'   `{rank}_unknown_{fread,fotu,focc}` for each non-root rank, `nread`,
#'   `nocc`, `notu`, `child_unknown_{fread,fotu,focc}`, `fread`, `focc`,
#'   `fotu`
#' @export
generate_krona_data <- function(otu_taxonomy, ranks = tax_ranks()) {
  # avoid R CMD check NOTE: no visible binding for global variable
  taxon <- parent <- nread <- nsample <- nocc <- notu <- NULL

  checkmate::assert_character(ranks, min.len = 2L, any.missing = FALSE)
  checkmate::assert_data_frame(otu_taxonomy)
  checkmate::assert_names(
    names(otu_taxonomy),
    must.include = c(ranks, "nread", "nsample")
  )

  if (nrow(otu_taxonomy) == 0L) {
    return(krona_empty_table(ranks))
  }

  unknown_ranks <- ranks[-1L]
  tax <- otu_taxonomy
  for (r in intersect(c("genus", "species"), ranks)) {
    tax[[r]] <- remove_mycobank_number(tax[[r]])
  }

  parent_path <- tax[[ranks[[1L]]]]
  for (i in seq_along(ranks)[-1L]) {
    tax[[paste0(ranks[[i]], "_parent")]] <- parent_path
    parent_path <- paste(parent_path, tax[[ranks[[i]]]], sep = ",")
  }
  for (r in unknown_ranks) {
    tax[[paste0(r, "_unknown")]] <- startsWith(tax[[r]], "pseudo")
  }
  for (r in ranks) {
    names(tax)[names(tax) == r] <- paste0(r, "_taxon")
  }

  unknown_vars <- paste0(unknown_ranks, "_unknown")
  out <- tax |>
    tidyr::pivot_longer(
      cols = c(
        dplyr::all_of(paste0(ranks, "_taxon")),
        dplyr::all_of(paste0(unknown_ranks, "_parent"))
      ),
      names_to = c("rank", ".value"),
      names_sep = "_",
      names_transform = list(rank = function(x) rank2factor(x, ranks))
    ) |>
    dplyr::mutate(taxon = chartr("_", " ", taxon)) |>
    dplyr::group_by(rank, taxon, parent) |>
    dplyr::summarize(
      dplyr::across(
        dplyr::all_of(unknown_vars),
        list(
          fread = ~ sum(nread * .) / sum(nread),
          fotu = ~ sum(.) / dplyr::n(),
          focc = ~ sum(nsample * .) / sum(nsample)
        ),
        .names = "{.col}_{.fn}"
      ),
      nread = sum(nread),
      nocc = sum(nsample),
      notu = dplyr::n(),
      .groups = "drop"
    )

  for (metric in c("fread", "fotu", "focc")) {
    child_col <- paste0("child_unknown_", metric)
    out[[child_col]] <- NA_real_
    for (i in seq_along(ranks)) {
      src_rank <- if (i < length(ranks)) {
        ranks[[i + 1L]]
      } else {
        ranks[[i]]
      }
      src_col <- paste0(src_rank, "_unknown_", metric)
      idx <- as.character(out$rank) == ranks[[i]]
      out[[child_col]][idx] <- out[[src_col]][idx]
    }
  }

  out |>
    dplyr::group_by(rank) |>
    dplyr::mutate(
      fread = nread / sum(nread),
      focc = nocc / sum(nocc),
      fotu = notu / sum(notu)
    ) |>
    dplyr::ungroup() |>
    dplyr::rename(parent_taxonomy = parent)
}

#' Internal function: convert a list of data to the XML format to be sent to KronaTools
#'
#' Converts `list(type1 = c(val1a, val1b), type2 = c(val2a, val2b))` to
#' `<type1><val>{val1a}</val><val>{val1b}</val></type1><type2><val>{val2a}</val><val>{val2b}</val></type2>`
#'
#' @param data_format (`list`) a named list whose elements are values to be
#' formatted as data for a KronaTools XML file
#' @return (`character`) the XML-formatted data
xml_format <- function(data_format) {
  lapply(data_format, vapply, sprintf, "", fmt = "<val>{%s}</val>") |>
    vapply(paste, "", collapse = ",") |>
    purrr::imap_chr(sprintf, fmt = "<%2$s>%1$s</%2$s>") |>
    paste(collapse = "\n")
}


#' Write a KronaTools file
#' @param data (`data.frame`) as generated by `generate_krona_data()`
#' @param .rank (`character` or `factor`) the rank to generate nodes for
#' @param maxrank (`factor`) the maximum rank to generate nodes for
#' @param outfile (`character` or `connection`) the file to write the XML to
#' @param pre (`character`) text to write before the XML nodes
#' @param post (`character`) text to write after the XML nodes
#' @param taxonomy (`character`) the taxonomic prefix to filter the data by
#' @param node_data_format (`list`) a named list whose elements are values to be
#' formatted as data for a KronaTools XML file
#' @param node_xml_format (`character`) the XML format to use for each node
#' @param ... additional arguments; passed during recursion but currently ignored
#' @return (`character`) the name of the file written
#' @export
krona_xml_nodes <- function(
  data,
  .rank,
  maxrank = rank2factor(tip_rank()),
  outfile,
  pre = NULL,
  post = NULL,
  taxonomy = paste(known_taxa(), collapse = ","),
  node_data_format = NULL,
  node_xml_format = xml_format(node_data_format),
  ...
) {
  # avoid R CMD check NOTE: no visible binding for global variable
  parent_taxonomy <- taxon <- NULL

  if (is.character(.rank)) {
    .rank <- rank2factor(.rank)
  }
  con <- outfile
  if (!methods::is(con, "connection")) {
    con <- file(con, open = "w")
    on.exit(close(con))
  }
  my_data <- data
  if (!is.null(taxonomy)) {
    my_data <- dplyr::filter(data, startsWith(parent_taxonomy, taxonomy))
  }
  xml <- dplyr::filter(my_data, rank == .rank) |>
    dplyr::transmute(
      taxon = taxon,
      taxonomy = ifelse(
        is.na(parent_taxonomy),
        taxon,
        paste(parent_taxonomy, taxon, sep = ",")
      ),
      pre = glue::glue(
        '<node name="{taxon}">',
        node_xml_format,
        .sep = "\n"
      ),
      post = "</node>"
    )
  if (!is.null(pre)) {
    writeLines(pre, con)
  }
  if (.rank == maxrank) {
    writeLines(paste(xml$pre, xml$post, sep = "\n"), con)
  } else {
    purrr::pwalk(
      xml,
      krona_xml_nodes,
      data = my_data,
      .rank = subranks(.rank)[1],
      maxrank = maxrank,
      outfile = con,
      ...,
      node_xml_format = node_xml_format
    )
  }
  if (!is.null(post)) {
    writeLines(post, con)
  }
  outfile
}
