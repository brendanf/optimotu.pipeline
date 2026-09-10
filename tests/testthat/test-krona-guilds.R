fungal_ranks <- c(
  "kingdom",
  "phylum",
  "class",
  "order",
  "family",
  "genus",
  "species"
)

metazoa_ranks <- c(
  "kingdom",
  "phylum",
  "class",
  "order",
  "family",
  "subfamily",
  "tribe",
  "genus",
  "species"
)

fungal_otus <- function() {
  tibble::tibble(
    seq_id = c("otu1", "otu2"),
    kingdom = c("Fungi", "Fungi"),
    phylum = c("Ascomycota", "Ascomycota"),
    class = c("Sordariomycetes", "Sordariomycetes"),
    order = c("Hypocreales", "Hypocreales"),
    family = c("Nectriaceae", "Nectriaceae"),
    genus = c("Fusarium_123", "pseudoNectriaceae"),
    species = c("Fusarium_oxysporum_456", "pseudoNectriaceae"),
    nread = c(90L, 10L),
    nsample = c(3L, 1L)
  )
}

metazoa_otus <- function() {
  tibble::tibble(
    seq_id = "otu1",
    kingdom = "Animalia",
    phylum = "Arthropoda",
    class = "Insecta",
    order = "Diptera",
    family = "Culicidae",
    subfamily = "Culicinae",
    tribe = "Aedini",
    genus = "Aedes",
    species = "Aedes_aegypti",
    nread = 50L,
    nsample = 2L
  )
}

test_that("generate_krona_data is rank-generic for fungal ranks", {
  out <- generate_krona_data(fungal_otus(), ranks = fungal_ranks)
  expect_true(all(
    paste0(
      rep(fungal_ranks[-1], each = 3L),
      "_unknown_",
      c("fread", "fotu", "focc")
    ) %in%
      names(out)
  ))
  expect_false("kingdom_unknown_fread" %in% names(out))

  genus_row <- dplyr::filter(out, rank == "genus", taxon == "Fusarium")
  expect_equal(nrow(genus_row), 1L)
  expect_equal(genus_row$nread, 90L)
  expect_equal(genus_row$species_unknown_fread, 0)

  unknown_genus <- dplyr::filter(
    out,
    rank == "genus",
    taxon == "pseudoNectriaceae"
  )
  expect_equal(unknown_genus$nread, 10L)
  expect_equal(unknown_genus$genus_unknown_fotu, 1)

  family_row <- dplyr::filter(out, rank == "family", taxon == "Nectriaceae")
  expect_equal(family_row$nread, 100L)
  expect_equal(family_row$notu, 2L)
  expect_equal(family_row$genus_unknown_fread, 0.1)
  expect_equal(family_row$child_unknown_fread, 0.1)
})

test_that("generate_krona_data handles metazoa ranks including subfamily/tribe", {
  out <- generate_krona_data(metazoa_otus(), ranks = metazoa_ranks)
  expect_true(all(
    c("subfamily_unknown_fread", "tribe_unknown_fread") %in% names(out)
  ))
  expect_setequal(
    as.character(unique(out$rank)),
    metazoa_ranks
  )
  tribe <- dplyr::filter(out, rank == "tribe")
  expect_equal(tribe$taxon, "Aedini")
  expect_equal(
    tribe$parent_taxonomy,
    "Animalia,Arthropoda,Insecta,Diptera,Culicidae,Culicinae"
  )
})

test_that("generate_krona_data returns a typed empty table", {
  empty <- fungal_otus()[0, ]
  out <- generate_krona_data(empty, ranks = fungal_ranks)
  expect_equal(nrow(out), 0L)
  expect_true("phylum_unknown_fread" %in% names(out))
  expect_s3_class(out$rank, "factor")
})

test_that("krona_node_data_format uses the last three unknown ranks", {
  fungal <- krona_node_data_format(fungal_ranks)
  expect_named(
    fungal,
    c("f", "nocc", "nread", "notu", "family", "genus", "species")
  )
  expect_equal(
    fungal$species,
    c("species_unknown_focc", "species_unknown_fread", "species_unknown_fotu")
  )

  metazoa <- krona_node_data_format(metazoa_ranks)
  expect_named(
    metazoa,
    c("f", "nocc", "nread", "notu", "tribe", "genus", "species")
  )
})

test_that("prepare_guild_taxonomy unites arbitrary ranks and strips suffixes", {
  out <- prepare_guild_taxonomy(fungal_otus(), ranks = fungal_ranks)
  expect_false("genus" %in% names(out))
  expect_true("Taxonomy" %in% names(out))
  expect_equal(
    out$Taxonomy[[1]],
    paste(
      "Fungi",
      "Ascomycota",
      "Sordariomycetes",
      "Hypocreales",
      "Nectriaceae",
      "Fusarium",
      "Fusarium_oxysporum",
      sep = ","
    )
  )

  meta <- prepare_guild_taxonomy(metazoa_otus(), ranks = metazoa_ranks)
  expect_match(meta$Taxonomy, "Culicinae,Aedini,Aedes,")
})

test_that("lifestyle_guild_db builds FUNGuild rows without Protax taxonomy", {
  path <- withr::local_tempfile(fileext = ".rds")
  saveRDS(
    tibble::tibble(
      taxon = c("Foo bar", "Foo"),
      guild = c("Pathogen Saprotroph", "Pathogen"),
      citationSource = c("paper", "paper")
    ),
    path
  )
  out <- lifestyle_guild_db(path)
  expect_named(
    out,
    c(
      "taxon",
      "taxonomicLevel",
      "trophicMode",
      "guild",
      "citationSource",
      "searchkey"
    )
  )
  expect_equal(out$taxonomicLevel[out$taxon == "Foo bar"], 20L)
  expect_equal(out$taxonomicLevel[out$taxon == "Foo"], 13L)
  expect_equal(out$guild[out$taxon == "Foo bar"], "Pathogen,Saprotroph")
  expect_false(any(c("classification", "rank", "prior") %in% names(out)))
})

test_that("has_lifestyle_guild_db reflects file presence", {
  missing <- file.path(tempdir(), "no-such-lifestyle.rds")
  expect_false(has_lifestyle_guild_db(missing))
  path <- withr::local_tempfile(fileext = ".rds")
  saveRDS(tibble::tibble(taxon = "Foo", guild = "Pathogen"), path)
  expect_true(has_lifestyle_guild_db(path))
})
