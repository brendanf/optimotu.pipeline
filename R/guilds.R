#' Default path to the Carlos/lifestyle guild table
#'
#' @return (`character`) relative path used by the pipeline
#' @export
lifestyle_guild_path <- function() {
  "data/lifestyle/Fung_LifeStyle_Data.RDS"
}

#' Whether the optional lifestyle guild RDS is present
#'
#' @param path (`character`) path to the RDS file
#' @return (`logical`) `TRUE` if `path` exists
#' @export
has_lifestyle_guild_db <- function(path = lifestyle_guild_path()) {
  file.exists(path)
}

read_guild_table <- function(path) {
  checkmate::assert_string(path)
  ext <- tolower(tools::file_ext(path))
  raw <- switch(
    ext,
    rds = readRDS(path),
    csv = readr::read_csv(path, show_col_types = FALSE),
    tsv = ,
    txt = readr::read_tsv(path, show_col_types = FALSE),
    fst = {
      if (!requireNamespace("fst", quietly = TRUE)) {
        stop("Reading an .fst guild database requires the fst package.")
      }
      fst::read_fst(path)
    },
    stop(
      "Unsupported guild database file extension '.",
      ext,
      "'. Use rds, tsv, csv, txt, or fst.",
      call. = FALSE
    )
  )
  checkmate::assert_data_frame(raw)
  if (!"guild" %in% names(raw)) {
    stop(
      "Guild database '",
      path,
      "' must have a 'guild' column.",
      call. = FALSE
    )
  }
  if (!"searchkey" %in% names(raw)) {
    if (!"taxon" %in% names(raw)) {
      stop(
        "Guild database '",
        path,
        "' must have a 'searchkey' column, or 'taxon' so one can be built.",
        call. = FALSE
      )
    }
    raw$searchkey <- paste0("@", sub("[_ ]", "@", raw$taxon), "@")
  }
  tibble::as_tibble(raw)
}

#' Load a guild database specified in `pipeline_options.yaml`
#'
#' @param source (`character`) one of `download`, `lifestyle`, or `file`
#' @param path (`character`) file path; ignored when `source` is `download`
#' @return `data.frame` suitable for `FUNGuildR::funguild_assign()`
#' @export
load_guild_database <- function(source, path = NULL) {
  source <- match.arg(source, c("download", "lifestyle", "file"))
  switch(
    source,
    download = {
      if (!requireNamespace("FUNGuildR", quietly = TRUE)) {
        stop("Package FUNGuildR is required for the 'funguild' database.")
      }
      FUNGuildR::get_funguild_db()
    },
    lifestyle = lifestyle_guild_db(path),
    file = read_guild_table(path)
  )
}

#' Prepare OTU taxonomy for FUNGuild-style assignment
#'
#' Strips trailing Mycobank-style `_[0-9]+` suffixes from `genus` and
#' `species` when those ranks exist, then unites rank columns into
#' `Taxonomy`.
#'
#' @param otu_taxonomy (`data.frame`) OTU taxonomy with rank columns
#' @param ranks (`character`) taxonomic ranks, most inclusive first
#' @return `data.frame` with a `Taxonomy` column and rank columns removed
#' @export
prepare_guild_taxonomy <- function(otu_taxonomy, ranks = tax_ranks()) {
  checkmate::assert_data_frame(otu_taxonomy)
  checkmate::assert_character(ranks, min.len = 1L, any.missing = FALSE)
  checkmate::assert_names(names(otu_taxonomy), must.include = ranks)

  out <- otu_taxonomy
  strip <- intersect(c("genus", "species"), ranks)
  if (length(strip) > 0L) {
    out <- dplyr::mutate(
      out,
      dplyr::across(
        dplyr::all_of(strip),
        \(x) sub("([A-Z].+)_[0-9]+", "\\1", x)
      )
    )
  }
  tidyr::unite(
    out,
    "Taxonomy",
    dplyr::all_of(ranks),
    sep = ","
  )
}

#' Convert a Carlos/lifestyle RDS into a FUNGuild-style table
#'
#' Uses only columns from `path` (typically `taxon`, `guild`,
#' `citationSource`). Does not join Protax `taxonomy_new`. Genus-level
#' rows are synthesized from species annotations when a genus-level
#' `taxon` is not already present.
#'
#' @param path (`character`) path to `Fung_LifeStyle_Data.RDS` or similar
#' @return `data.frame` with FUNGuild columns `taxon`, `taxonomicLevel`,
#'   `trophicMode`, `guild`, `citationSource`, `searchkey`
#' @export
lifestyle_guild_db <- function(path) {
  # avoid R CMD check NOTE: no visible binding for global variable
  taxon <- guild <- genus <- citationSource <- NULL

  checkmate::assert_string(path)
  raw <- readRDS(path)
  checkmate::assert_data_frame(raw)
  checkmate::assert_names(names(raw), must.include = c("taxon", "guild"))
  if (!"citationSource" %in% names(raw)) {
    raw$citationSource <- NA_character_
  }

  raw <- dplyr::mutate(
    raw,
    genus = sub(" .*", "", taxon),
    guild = sub("Lichenized_Saprotroph", "Lichenized Saprotroph", guild),
    guild = sub(
      "Lichen_Parasite_Saprotroph",
      "Lichen_Parasite Saprotroph",
      guild
    )
  )

  species_rows <- dplyr::transmute(
    raw,
    taxon,
    taxonomicLevel = ifelse(grepl(" ", taxon, fixed = TRUE), 20L, 13L),
    trophicMode = NA_character_,
    guild = chartr(" ", ",", guild),
    citationSource,
    searchkey = paste0("@", sub("[_ ]", "@", taxon), "@")
  )

  genus_rows <- dplyr::summarize(
    raw,
    guild = paste(
      setdiff(
        unique(unlist(strsplit(guild, "[ ,]"))),
        c("NA", NA_character_)
      ),
      collapse = ","
    ),
    .by = genus
  ) |>
    dplyr::transmute(
      taxon = genus,
      taxonomicLevel = 13L,
      trophicMode = NA_character_,
      guild,
      citationSource = "combined from species-level annotations",
      searchkey = paste0("@", taxon, "@")
    ) |>
    dplyr::anti_join(raw, by = "taxon")

  dplyr::bind_rows(species_rows, genus_rows)
}
