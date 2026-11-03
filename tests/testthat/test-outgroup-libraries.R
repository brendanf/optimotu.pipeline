test_that("outgroup_library_taxonomy prefixes BOLD and gzipped SINTAX ids", {
  ranks <- c("kingdom", "phylum", "genus", "species")
  bold <- withr::local_tempfile(fileext = ".fasta")
  writeLines(
    c(
      paste0(
        ">seqA|BIN1|Animalia,Arthropoda,Insecta,Diptera,FamA,Sub,Tribe,",
        "GenusA,speciesA"
      ),
      "ACGTACGTACGT"
    ),
    bold
  )
  bact <- withr::local_tempfile(fileext = ".fasta.gz")
  con <- gzfile(bact, "wt")
  writeLines(
    c(
      paste0(
        ">seqB;tax=k:Bacteria,p:Bacillota,c:Bacilli,o:Caryophanales,",
        "f:FamB,g:GenusB,s:speciesB;gene=ctaD"
      ),
      "TGCATGCATGCA"
    ),
    con
  )
  close(con)
  libraries <- list(
    list(name = "bold", sequences = bold, taxonomy = NULL),
    list(name = "bacteria", sequences = bact, taxonomy = NULL)
  )
  tax <- outgroup_library_taxonomy(libraries, ranks = ranks)
  expect_equal(tax$ref_id, c("bold:seqA", "bacteria:seqB"))
  expect_equal(tax$kingdom, c("Animalia", "Bacteria"))
  expect_equal(tax$genus, c("GenusA", "GenusB"))
})

test_that("outgroup_library_taxonomy errors on duplicate ids", {
  ranks <- c("kingdom", "genus")
  fasta <- withr::local_tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">seqA|BIN1|Animalia,Arthropoda,Insecta,Diptera,FamA,Sub,Tribe,GenusA,spA",
      "ACGT",
      ">seqA|BIN2|Animalia,Arthropoda,Insecta,Diptera,FamA,Sub,Tribe,GenusA,spB",
      "TGCA"
    ),
    fasta
  )
  expect_error(
    outgroup_library_taxonomy(
      list(list(name = "bold", sequences = fasta, taxonomy = NULL)),
      ranks = ranks
    ),
    "Duplicate sequence ids"
  )
})

test_that("reduce_library_best_hits keeps the closest hit", {
  hits <- dplyr::bind_rows(
    tibble::tibble(seq_idx = 1L, ref_id = "bold:seqA", dist = 0.2),
    tibble::tibble(seq_idx = 1L, ref_id = "bacteria:seqB", dist = 0.05)
  )
  taxonomy <- tibble::tibble(
    ref_id = c("bold:seqA", "bacteria:seqB"),
    kingdom = c("Animalia", "Bacteria")
  )
  out <- reduce_library_best_hits(hits, taxonomy, id_col = "seq_idx")
  expect_equal(nrow(out), 1L)
  expect_equal(out$ref_id, "bacteria:seqB")
  expect_equal(out$kingdom, "Bacteria")
})

test_that("outgroup_hit_ref_ids parses mixed prefixed headers", {
  headers <- c(
    "bold:seqA|BIN1|Animalia,Arthropoda,Insecta,Diptera,FamA,Sub,Tribe,GenusA,spA",
    paste0(
      "bacteria:seqB;tax=k:Bacteria,p:Bacillota,c:Bacilli,o:Caryophanales,",
      "f:FamB,g:GenusB,s:speciesB;gene=ctaD"
    )
  )
  expect_equal(
    outgroup_hit_ref_ids(headers),
    c("bold:seqA", "bacteria:seqB")
  )
})

test_that("streamed uchime_ref and usearch_global see both libraries", {
  skip_if_no_vsearch()
  rng <- function(n, seed) {
    set.seed(seed)
    paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")
  }
  left <- rng(300L, 1L)
  right <- rng(300L, 2L)
  chimera <- paste0(substr(left, 1L, 150L), substr(right, 151L, 300L))
  bold <- withr::local_tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">seqA|BIN|Animalia,Arthropoda,Insecta,Diptera,FamA,Sub,Tribe,GenusA,spA",
      left
    ),
    bold
  )
  bact <- withr::local_tempfile(fileext = ".fasta.gz")
  con <- gzfile(bact, "wt")
  writeLines(
    c(
      ">seqB;tax=k:Bacteria,p:Bacillota,c:Bacilli,o:Caryophanales,f:FamB,g:GenusB,s:speciesB",
      right
    ),
    con
  )
  close(con)
  query <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(">chim", chimera), query)
  libraries <- list(
    list(name = "bold", sequences = bold, taxonomy = NULL),
    list(name = "bacteria", sequences = bact, taxonomy = NULL)
  )
  both <- vsearch_uchime_ref(
    query,
    libraries,
    id_only = TRUE,
    vsearch = find_vsearch()
  )
  expect_equal(both, "chim")
  expect_equal(
    vsearch_uchime_ref(query, bold, id_only = TRUE, vsearch = find_vsearch()),
    character()
  )
  expect_equal(
    vsearch_uchime_ref(query, bact, id_only = TRUE, vsearch = find_vsearch()),
    character()
  )
  hits <- vsearch_usearch_global(
    query,
    libraries,
    threshold = 0.5,
    vsearch = find_vsearch()
  )
  expect_true(nrow(hits) >= 1L)
  expect_true(all(
    outgroup_hit_ref_ids(hits$cluster) %in% c("bold:seqA", "bacteria:seqB")
  ))
})
