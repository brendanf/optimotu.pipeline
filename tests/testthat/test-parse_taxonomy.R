# Fixture headers (one per format) for parse_taxonomy tests
sintax_headers <- c(
  "seq1;tax=k:Fungi,p:Ascomycota,c:Dothideomycetes,o:Pleosporales,f:Pleosporaceae,g:Alternaria,s:Alternaria alternata",
  "seq2;other=foo;tax=d:Eukaryota,k:Protista,p:Ochrophyta,c:Bacillariophyta,o:Naviculales,f:Naviculaceae,g:Navicula,s:Navicula sp."
)
bold_headers <- c(
  "BOLD:ABC123|COI-5P|Canada|Animalia,Arthropoda,Insecta,Coleoptera,Carabidae,,,Pterostichus,Pterostichus sp.",
  "BOLD:DEF456|COI-5P|USA|Animalia,Chordata,Mammalia,Carnivora,Felidae,,,Felis,Felis catus"
)
unite_headers <- c(
  "SH123456.01FU|k__Fungi;p__Ascomycota;c__Dothideomycetes;o__Pleosporales;f__Pleosporaceae;g__Alternaria;s__Alternaria alternata|SH123456",
  "SH789012.01FU|k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales;f__Agaricaceae;g__Agaricus;s__Agaricus bisporus|SH789012"
)
bayesant_headers <- c(
  "seqBA1 Root;Fungi;Ascomycota;Dothideomycetes;Pleosporales;Pleosporaceae;Alternaria;Alternaria alternata",
  "seqBA2 Root;Fungi;Basidiomycota;Agaricomycetes;Agaricales;Agaricaceae;Agaricus;Agaricus bisporus"
)

# Rank map for testing
# Encodes BOLD ranks to the set supported by sintax
rank_map <- data.frame(
  BOLD = c(
    "class",
    "order",
    "family",
    "subfamily",
    "tribe",
    "genus",
    "species"
  ),
  sintax = c(
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
  )
)

# Headers with BOLD ranks encoded as sintax ranks
bold_sintax_headers <- c(
  "BOLD:ABC123;tax=k:Insecta,p:Coleoptera,c:Carabidae,o:Pterostichinae,f:Pterostichini,g:Pterostichus,s:Pterostichus sp.",
  "BOLD:DEF456;tax=k:Mammalia,p:Carnivora,c:Felidae,g:Felis,s:Felis catus"
)

test_that("is_sintax_header detects sintax format", {
  expect_true(optimotu.pipeline:::is_sintax_header(sintax_headers))
  expect_false(optimotu.pipeline:::is_sintax_header(bold_headers))
  expect_false(optimotu.pipeline:::is_sintax_header(unite_headers))
})

test_that("is_bold_header detects BOLD format", {
  expect_true(optimotu.pipeline:::is_bold_header(bold_headers))
  expect_false(optimotu.pipeline:::is_bold_header(sintax_headers))
  expect_false(optimotu.pipeline:::is_bold_header(unite_headers))
})

test_that("is_unite_header detects UNITE format", {
  expect_true(optimotu.pipeline:::is_unite_header(unite_headers))
  expect_false(optimotu.pipeline:::is_unite_header(sintax_headers))
  expect_false(optimotu.pipeline:::is_unite_header(bold_headers))
})

test_that("is_bayesant_header detects BayesANT format", {
  expect_true(optimotu.pipeline:::is_bayesant_header(bayesant_headers))
  expect_false(optimotu.pipeline:::is_bayesant_header(sintax_headers))
  expect_false(optimotu.pipeline:::is_bayesant_header(bold_headers))
  expect_false(optimotu.pipeline:::is_bayesant_header(unite_headers))
})

test_that("parse_sintax_header returns seq_id and rank columns", {
  ranks <- tax_ranks()
  out <- optimotu.pipeline:::parse_sintax_header(sintax_headers)
  checkmate::expect_names(
    names(out),
    must.include = c("seq_id", "kingdom", "species")
  )
  expect_equal(out$seq_id, c("seq1", "seq2"))
  expect_equal(out$kingdom[1], "Fungi")
  expect_equal(out$species[1], "Alternaria alternata")
})

test_that("parse_bold_header returns seq_id and rank columns", {
  ranks <- tax_ranks()
  out <- optimotu.pipeline:::parse_bold_header(bold_headers)
  checkmate::expect_names(
    names(out),
    must.include = c("seq_id", "kingdom", "genus")
  )
  expect_equal(out$seq_id[1], "BOLD:ABC123")
  expect_equal(out$kingdom[1], "Animalia")
  expect_equal(out$genus[1], "Pterostichus")
})

test_that("parse_unite_header returns seq_id and rank columns", {
  ranks <- tax_ranks()
  out <- optimotu.pipeline:::parse_unite_header(unite_headers)
  checkmate::expect_names(
    names(out),
    must.include = c("seq_id", "kingdom", "species")
  )
  expect_equal(out$seq_id[1], "SH123456.01FU")
  expect_equal(out$kingdom[1], "Fungi")
  expect_equal(out$species[1], "Alternaria alternata")
})

test_that("parse_bayesant_header returns seq_id and rank columns", {
  ranks <- tax_ranks()
  out <- optimotu.pipeline:::parse_bayesant_header(bayesant_headers, ranks)
  checkmate::expect_names(
    names(out),
    must.include = c("seq_id", "kingdom", "species")
  )
  expect_equal(out$seq_id, c("seqBA1", "seqBA2"))
  expect_equal(out$kingdom, c("Fungi", "Fungi"))
  expect_equal(out$species, c("Alternaria alternata", "Agaricus bisporus"))
})

test_that("parse_taxonomy_header detects format and dispatches", {
  ranks <- tax_ranks()
  out_sintax <- optimotu.pipeline:::parse_taxonomy_header(sintax_headers, ranks)
  expect_false(is.null(out_sintax))
  checkmate::expect_names(names(out_sintax), must.include = c("seq_id", ranks))
  expect_equal(out_sintax$seq_id[1], "seq1")

  out_bold <- optimotu.pipeline:::parse_taxonomy_header(bold_headers, ranks)
  expect_false(is.null(out_bold))
  expect_equal(out_bold$seq_id[1], "BOLD:ABC123")

  out_unite <- optimotu.pipeline:::parse_taxonomy_header(unite_headers, ranks)
  expect_false(is.null(out_unite))
  expect_equal(out_unite$seq_id[1], "SH123456.01FU")

  out_bayesant <- optimotu.pipeline:::parse_taxonomy_header(
    bayesant_headers,
    ranks
  )
  expect_false(is.null(out_bayesant))
  expect_equal(out_bayesant$seq_id[1], "seqBA1")
})

test_that("parse_taxonomy_header returns NULL for unrecognized format", {
  bad_headers <- c("just_an_id", "another_plain_id")
  expect_null(optimotu.pipeline:::parse_taxonomy_header(
    bad_headers,
    tax_ranks()
  ))
})

test_that("parse_taxonomy_tabular with rank columns", {
  lines <- c(
    "seq_id\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies",
    "s1\tFungi\tAscomycota\tDothideomycetes\tPleosporales\tPleosporaceae\tAlternaria\tAlternaria alternata"
  )
  out <- optimotu.pipeline:::parse_taxonomy_tabular(I(lines), tax_ranks())
  expect_equal(nrow(out), 1)
  expect_equal(out$seq_id, "s1")
  expect_equal(out$kingdom, "Fungi")
  expect_equal(out$species, "Alternaria alternata")
})

test_that("parse_taxonomy_tabular with taxonomy column (comma-delimited)", {
  lines <- c(
    "seq_id\ttaxonomy",
    "s1\tFungi,Ascomycota,Dothideomycetes,Pleosporales,Pleosporaceae,Alternaria,Alternaria alternata"
  )
  out <- optimotu.pipeline:::parse_taxonomy_tabular(I(lines), tax_ranks())
  expect_equal(nrow(out), 1)
  expect_equal(out$seq_id, "s1")
  expect_equal(out$kingdom, "Fungi")
  expect_equal(out$species, "Alternaria alternata")
})

test_that("parse_taxonomy_tabular with taxonomy column (sintax-style)", {
  lines <- c(
    "id\ttaxonomy",
    "s1\tseq1;tax=k:Fungi,p:Ascomycota,c:Dothideomycetes,o:Pleosporales,f:Pleosporaceae,g:Alternaria,s:Alternaria alternata"
  )
  out <- optimotu.pipeline:::parse_taxonomy_tabular(I(lines), tax_ranks())
  expect_equal(nrow(out), 1)
  expect_equal(out$seq_id, "s1")
  expect_equal(out$kingdom, "Fungi")
})

test_that("parse_taxonomy_tabular empty lines returns empty tibble", {
  out <- optimotu.pipeline:::parse_taxonomy_tabular(
    I(character(0)),
    tax_ranks()
  )
  expect_equal(nrow(out), 0)
  expect_true(all(c("seq_id", tax_ranks()) %in% names(out)))
})

test_that("parse_taxonomy_tabular errors when format not detected", {
  lines <- c("id\tfoo\tbar\tbaf", "1\ta\tb")
  expect_error(
    optimotu.pipeline:::parse_taxonomy_tabular(I(lines), tax_ranks()),
    "Could not detect tabular format"
  )
})

test_that("parse_taxonomy_tabular headerless: seq_id + rank columns in order", {
  ranks <- tax_ranks()
  lines <- c(
    "s1\tFungi\tAscomycota\tDothideomycetes\tPleosporales\tPleosporaceae\tAlternaria\tAlternaria alternata",
    "s2\tProtista\tOchrophyta\tBacillariophyta\tNaviculales\tNaviculaceae\tNavicula\tNavicula sp."
  )
  out <- optimotu.pipeline:::parse_taxonomy_tabular(I(lines), ranks)
  expect_equal(nrow(out), 2)
  expect_equal(names(out), c("seq_id", ranks))
  expect_equal(out$seq_id, c("s1", "s2"))
  expect_equal(out$kingdom, c("Fungi", "Protista"))
  expect_equal(out$species, c("Alternaria alternata", "Navicula sp."))
})

test_that("parse_taxonomy_tabular headerless: seq_id + taxonomy column", {
  lines <- c(
    "s1\tFungi,Ascomycota,Dothideomycetes,Pleosporales,Pleosporaceae,Alternaria,Alternaria alternata",
    "s2\tProtista;Ochrophyta;Bacillariophyta;Naviculales;Naviculaceae;Navicula;Navicula sp."
  )
  out <- optimotu.pipeline:::parse_taxonomy_tabular(I(lines), tax_ranks())
  expect_equal(nrow(out), 2)
  expect_equal(out$seq_id, c("s1", "s2"))
  expect_equal(out$kingdom, c("Fungi", "Protista"))
  expect_equal(out$species, c("Alternaria alternata", "Navicula sp."))
})

test_that("parse_taxonomy_header with rank map", {
  out <- optimotu.pipeline:::parse_taxonomy_header(
    bold_sintax_headers,
    c("class", "order", "family", "subfamily", "tribe", "genus", "species"),
    rank_map
  )
  checkmate::expect_names(
    names(out),
    must.include = c(
      "seq_id",
      "class",
      "order",
      "family",
      "subfamily",
      "tribe",
      "genus",
      "species"
    )
  )
  expect_equal(out$seq_id[1], "BOLD:ABC123")
  expect_equal(out$class[1], "Insecta")
  expect_equal(out$order[1], "Coleoptera")
  expect_equal(out$family[1], "Carabidae")
  expect_equal(out$subfamily[1], "Pterostichinae")
  expect_equal(out$tribe[1], "Pterostichini")
  expect_equal(out$genus[1], "Pterostichus")
  expect_equal(out$species[1], "Pterostichus sp.")
})

test_that("parse_reference_taxonomy detects format and dispatches", {
  dna <- c("AGCTCTACTATCTATTACTGTGGGCTCGTAG", "GCTAGCGATGCTGAGCGTATCGATGC")
  test_file <- withr::local_tempfile(fileext = ".fa")
  stats::setNames(dna, sintax_headers) |>
    Biostrings::DNAStringSet() |>
    Biostrings::writeXStringSet(test_file)
  out_sintax <- optimotu.pipeline:::parse_reference_taxonomy(
    test_file,
    tax_ranks()
  )
  expect_false(is.null(out_sintax))
  checkmate::expect_names(
    names(out_sintax),
    must.include = c("seq_id", tax_ranks())
  )
  checkmate::expect_names(names(out_sintax), must.include = tax_ranks())
  expect_equal(out_sintax$seq_id[1], "seq1")
  expect_equal(out_sintax$kingdom[1], "Fungi")
  expect_equal(out_sintax$species[1], "Alternaria alternata")

  test_file <- withr::local_tempfile(fileext = ".fa")
  stats::setNames(dna, bold_headers) |>
    Biostrings::DNAStringSet() |>
    Biostrings::writeXStringSet(test_file)
  out_bold <- optimotu.pipeline:::parse_reference_taxonomy(
    test_file,
    tax_ranks()
  )
  expect_false(is.null(out_bold))
  checkmate::expect_names(
    names(out_bold),
    must.include = c("seq_id", tax_ranks())
  )
  expect_equal(out_bold$seq_id[1], "BOLD:ABC123")
  expect_equal(out_bold$kingdom[1], "Animalia")
  expect_equal(out_bold$genus[1], "Pterostichus")

  test_file <- withr::local_tempfile(fileext = ".fa")
  stats::setNames(dna, unite_headers) |>
    Biostrings::DNAStringSet() |>
    Biostrings::writeXStringSet(test_file)
  out_unite <- optimotu.pipeline:::parse_reference_taxonomy(
    test_file,
    tax_ranks()
  )
  expect_false(is.null(out_unite))
  checkmate::expect_names(
    names(out_unite),
    must.include = c("seq_id", tax_ranks())
  )
  expect_equal(out_unite$seq_id[1], "SH123456.01FU")
  expect_equal(out_unite$kingdom[1], "Fungi")
  expect_equal(out_unite$species[1], "Alternaria alternata")

  test_file <- withr::local_tempfile(fileext = ".fa")
  stats::setNames(dna, bayesant_headers) |>
    Biostrings::DNAStringSet() |>
    Biostrings::writeXStringSet(test_file)
  out_bayesant <- optimotu.pipeline:::parse_reference_taxonomy(
    test_file,
    tax_ranks()
  )
  expect_false(is.null(out_bayesant))
  checkmate::expect_names(
    names(out_bayesant),
    must.include = c("seq_id", tax_ranks())
  )
  expect_equal(out_bayesant$seq_id[1], "seqBA1")
  expect_equal(out_bayesant$kingdom[1], "Fungi")
  expect_equal(out_bayesant$species[1], "Alternaria alternata")

  test_file <- withr::local_tempfile(fileext = ".tsv")
  writeLines(
    c(
      "seq_id\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies",
      "s1\tFungi\tAscomycota\tDothideomycetes\tPleosporales\tPleosporaceae\tAlternaria\tAlternaria alternata"
    ),
    test_file
  )
  out_tabular <- optimotu.pipeline:::parse_reference_taxonomy(
    test_file,
    tax_ranks()
  )
  expect_false(is.null(out_tabular))
  checkmate::expect_names(
    names(out_tabular),
    must.include = c("seq_id", tax_ranks())
  )
  expect_equal(out_tabular$seq_id[1], "s1")
  expect_equal(out_tabular$kingdom[1], "Fungi")
  expect_equal(out_tabular$species[1], "Alternaria alternata")
})
