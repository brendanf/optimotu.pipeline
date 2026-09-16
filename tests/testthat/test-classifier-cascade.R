test_that("assigned_ranks drops scoped ranks", {
  expect_equal(
    assigned_ranks(c("kingdom"), c("kingdom", "phylum", "class")),
    c("phylum", "class")
  )
  expect_equal(
    assigned_ranks(character(), c("kingdom", "phylum")),
    c("kingdom", "phylum")
  )
})

test_that("ids_needing_reclassify uses known-prob gate and treats unk as failure", {
  tax <- tibble::tibble(
    seq_idx = c(1L, 1L, 2L, 3L, 4L),
    rank = c("family", "family", "family", "family", "genus"),
    taxon = c("FamA", "unk", "unk", NA_character_, "GenA"),
    prob = c(0.95, 0.9, 0.99, NA_real_, 0.99)
  )
  need <- ids_needing_reclassify(
    tax_prob = tax,
    rank = "family",
    min_prob = 0.9,
    all_ids = 1:5,
    id_col = "seq_idx"
  )
  # 1 keeps FamA >= 0.9; 2 is unk; 3 missing taxon; 4 no family row; 5 no rows
  expect_equal(sort(need), c(2L, 3L, 4L, 5L))
})

test_that("ids_needing_reclassify unions best_hit scope screen", {
  tax <- tibble::tibble(
    seq_idx = c(1L, 2L),
    rank = c("family", "family"),
    taxon = c("FamA", "FamB"),
    prob = c(0.99, 0.99)
  )
  screen <- tibble::tibble(
    seq_idx = c(1L, 2L, 3L),
    kingdom = c("Animalia", "Plantae", "Fungi")
  )
  need <- ids_needing_reclassify(
    tax_prob = tax,
    rank = "family",
    min_prob = 0.9,
    all_ids = 1:3,
    id_col = "seq_idx",
    screen = screen,
    scope_ranks = "kingdom",
    scope_taxa = "Animalia"
  )
  # 1 confident + Animalia keep; 2 confident but Plantae screen; 3 no tax + Fungi
  expect_equal(sort(need), c(2L, 3L))
})

test_that("ids_needing_reclassify ignores best_hit exception labels", {
  tax <- tibble::tibble(
    seq_idx = 1L,
    rank = "family",
    taxon = "FamA",
    prob = 0.99
  )
  screen <- tibble::tibble(
    seq_idx = 1L,
    kingdom = "unspecified"
  )
  need <- ids_needing_reclassify(
    tax_prob = tax,
    rank = "family",
    min_prob = 0.9,
    all_ids = 1L,
    id_col = "seq_idx",
    screen = screen,
    scope_ranks = "kingdom",
    scope_taxa = "Animalia"
  )
  expect_equal(need, integer())
})

test_that("combine_tax_prob replaces whole ids and tags classifier", {
  primary <- tibble::tibble(
    seq_idx = c(1L, 1L, 2L),
    rank = c("family", "genus", "family"),
    taxon = c("FamA", "GenA", "FamB"),
    prob = c(0.9, 0.8, 0.7)
  )
  fallback <- tibble::tibble(
    seq_idx = c(2L, 2L),
    rank = c("kingdom", "family"),
    taxon = c("Plantae", "FamX"),
    prob = c(0.95, 0.6)
  )
  out <- combine_tax_prob(
    primary,
    fallback,
    reclassify_ids = 2L,
    id_col = "seq_idx",
    primary_name = "protax",
    fallback_name = "epa"
  )
  expect_equal(sort(unique(out$seq_idx)), c(1L, 2L))
  expect_equal(unique(out$classifier[out$seq_idx == 1L]), "protax")
  expect_equal(unique(out$classifier[out$seq_idx == 2L]), "epa")
  expect_false(any(out$seq_idx == 2L & out$taxon == "FamB"))
  expect_true(any(out$seq_idx == 2L & out$taxon == "Plantae"))
})
