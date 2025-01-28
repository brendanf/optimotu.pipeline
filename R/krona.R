# TODO: let this use different ranks
#' Generate data to write a KronaTools XML file
#' @param otu_taxonomy (`data.frame`) with (at least) columns `kingdom`, `phylum`,
#' `class`, `order`, `family`, `genus`, `species`, `nread` and `nsample`
#' @return `data.frame` with columns `rank`, `taxon`, `parent_taxonomy`,
#' `phylum_unknown_fread`, `phylum_unknown_fotu`, `phylum_unknown_focc`,
#' `class_unknown_fread`, `class_unknown_fotu`, `class_unknown_focc`,
#' `order_unknown_fread`, `order_unknown_fotu`, `order_unknown_focc`,
#' `family_unknown_fread`, `family_unknown_fotu`, `family_unknown_focc`,
#' `genus_unknown_fread`, `genus_unknown_fotu`, `genus_unknown_focc`,
#' `species_unknown_fread`, `species_unknown_fotu`, `species_unknown_focc`,
#' `nread`, `nocc`, `notu`, `child_unknown_fread`, `child_unknown_fotu`,
#' `child_unknown_focc`, `fread`, `focc`, `fotu`
#' @export
generate_krona_data <- function(otu_taxonomy) {
  # avoid R CMD check NOTE: no visible binding for global variable
  genus <- species <- kingdom <- phylum_parent <- phylum <- class_parent <-
    order_parent <- family_parent <- family <- genus_parent <- kingdom_taxon <-
    species_parent <- taxon <- parent <- phylum_unknown <- species_unknown <-
    nread <- nsample <- nocc <- notu <- NULL

  if (nrow(otu_taxonomy) == 0) {
    tibble::tibble(
      rank = rank2factor(character()),
      taxon = character(),
      parent_taxonomy = character(),
      phylum_unknown_fread = numeric(),
      phylum_unknown_fotu = numeric(),
      phylum_unknown_focc = numeric(),
      class_unknown_fread = numeric(),
      class_unknown_fotu = numeric(),
      class_unknown_focc = numeric(),
      order_unknown_fread = numeric(),
      order_unknown_fotu = numeric(),
      order_unknown_focc = numeric(),
      family_unknown_fread = numeric(),
      family_unknown_fotu = numeric(),
      family_unknown_focc = numeric(),
      genus_unknown_fread = numeric(),
      genus_unknown_fotu = numeric(),
      genus_unknown_focc = numeric(),
      species_unknown_fread = numeric(),
      species_unknown_fotu = numeric(),
      species_unknown_focc = numeric(),
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
  } else {
    otu_taxonomy |>
      dplyr::mutate(
        genus = remove_mycobank_number(genus),
        species = remove_mycobank_number(species),
        phylum_parent = kingdom,
        class_parent = paste(phylum_parent, phylum, sep = ","),
        order_parent = paste(class_parent, class, sep = ","),
        family_parent = paste(order_parent, order, sep = ","),
        genus_parent = paste(family_parent, family, sep = ","),
        species_parent = paste(genus_parent, genus, sep = ","),
        dplyr::across(
          .cols = phylum:species,
          .fns = \(x) startsWith(x, "pseudo"),
          .names = "{.col}_unknown"
        )
      ) |>
      dplyr::rename_with(.fn = paste0, .cols = kingdom:species, "_taxon") |>
      tidyr::pivot_longer(
        kingdom_taxon:species_parent,
        names_to = c("rank", ".value"),
        names_sep = "_",
        names_transform = list(rank = rank2factor)
      ) |>
      dplyr::mutate(taxon = chartr("_", " ", taxon)) |>
      dplyr::group_by(rank, taxon, parent) |>
      dplyr::summarize(
        dplyr::across(
          phylum_unknown:species_unknown,
          list(
            fread = ~sum(nread * .)/sum(nread),
            fotu = ~sum(.)/dplyr::n(),
            focc = ~sum(nsample*.)/(sum(nsample))
          ),
          .names = "{.col}_{.fn}"
        ),
        nread = sum(nread),
        nocc = sum(nsample),
        notu = dplyr::n(),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        child_unknown_fread = dplyr::case_when(
          rank == "kingdom" ~ phylum_unknown_fread,
          rank == "phylum" ~ class_unknown_fread,
          rank == "class" ~ order_unknown_fread,
          rank == "order" ~ family_unknown_fread,
          rank == "family" ~ genus_unknown_fread,
          TRUE ~ species_unknown_fread
        ),
        child_unknown_focc = dplyr::case_when(
          rank == "kingdom" ~ phylum_unknown_focc,
          rank == "phylum" ~ class_unknown_focc,
          rank == "class" ~ order_unknown_focc,
          rank == "order" ~ family_unknown_focc,
          rank == "family" ~ genus_unknown_focc,
          TRUE ~ species_unknown_focc
        ),
        child_unknown_fotu = dplyr::case_when(
          rank == "kingdom" ~ phylum_unknown_fotu,
          rank == "phylum" ~ class_unknown_fotu,
          rank == "class" ~ order_unknown_fotu,
          rank == "order" ~ family_unknown_fotu,
          rank == "family" ~ genus_unknown_fotu,
          TRUE ~ species_unknown_fotu
        )
      ) |>
      dplyr::group_by(rank) |>
      dplyr::mutate(
        fread = nread/sum(nread),
        focc = nocc/sum(nocc),
        fotu = notu/sum(notu)
      ) |>
      dplyr::rename(parent_taxonomy = parent)
  }
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
    purrr::imap_chr(sprintf, fmt="<%2$s>%1$s</%2$s>") |>
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

  if (is.character(.rank)) .rank <- rank2factor(.rank)
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
      taxonomy = ifelse(is.na(parent_taxonomy), taxon, paste(parent_taxonomy, taxon, sep = ",")),
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
  if (!is.null(post)) writeLines(post, con)
  outfile
}
