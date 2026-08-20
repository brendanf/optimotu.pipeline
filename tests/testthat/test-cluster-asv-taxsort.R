tax_ranks_fixture <- c(
  "kingdom",
  "phylum",
  "class",
  "order",
  "family",
  "genus",
  "species"
)

known_taxon_idx_fixture <- function() {
  tibble::tibble(
    seq_idx = c(1L, 2L, 3L),
    genus = rep("Homo", 3L),
    species = c("sapiens", NA, NA)
  )
}

known_taxon_id_fixture <- function() {
  tibble::tibble(
    seq_id = c("ASV0001", "ASV0002", "ASV0003"),
    genus = rep("Homo", 3L),
    species = c("sapiens", NA, NA)
  )
}

asv_taxsort_idx_fixture <- function() {
  tibble::tibble(
    seq_idx_in = 1:3,
    seq_idx = 3:1,
    kingdom = rep("Euk", 3L),
    phylum = rep("Chordata", 3L),
    class = rep("Mammalia", 3L),
    order = rep("Primates", 3L),
    family = rep("Hominidae", 3L)
  )
}

asv_taxsort_id_fixture <- function() {
  tibble::tibble(
    seq_id = c("ASV0001", "ASV0002", "ASV0003"),
    seq_idx = 3:1,
    kingdom = rep("Euk", 3L),
    phylum = rep("Chordata", 3L),
    class = rep("Mammalia", 3L),
    order = rep("Primates", 3L),
    family = rep("Hominidae", 3L)
  )
}

closedref_taxon_idx_fixture <- function() {
  known_taxon_idx_fixture() |>
    dplyr::filter(is.na(species))
}

closedref_taxon_id_fixture <- function() {
  known_taxon_id_fixture() |>
    dplyr::filter(is.na(species))
}

test_that("full_preclosed joins on seq_idx_in with taxon table seq_idx", {
  out <- full_preclosed_taxon_table(
    known_taxon_table = known_taxon_idx_fixture(),
    asv_taxsort = asv_taxsort_idx_fixture(),
    rank = "species",
    parent_rank = "genus",
    tax_ranks = tax_ranks_fixture
  )
  expect_true("seq_idx" %in% names(out))
  expect_true("seq_idx_in" %in% names(out))
  expect_gt(nrow(out), 0L)
})

test_that("full_preclosed joins on seq_idx_in with taxon table seq_id only", {
  out <- full_preclosed_taxon_table(
    known_taxon_table = known_taxon_id_fixture(),
    asv_taxsort = asv_taxsort_idx_fixture(),
    rank = "species",
    parent_rank = "genus",
    tax_ranks = tax_ranks_fixture
  )
  expect_true("seq_idx" %in% names(out))
  expect_false("seq_idx_in" %in% names(out))
  expect_gt(nrow(out), 0L)
})

test_that("full_preclosed joins on seq_id", {
  out <- full_preclosed_taxon_table(
    known_taxon_table = known_taxon_id_fixture(),
    asv_taxsort = asv_taxsort_id_fixture(),
    rank = "species",
    parent_rank = "genus",
    tax_ranks = tax_ranks_fixture
  )
  expect_true(all(c("seq_id", "seq_idx") %in% names(out)))
  expect_false("seq_idx_in" %in% names(out))
  expect_gt(nrow(out), 0L)
})

test_that("full_preclosed errors when asv_taxsort uses seq_id but taxon table has seq_idx only", {
  expect_error(
    full_preclosed_taxon_table(
      known_taxon_table = known_taxon_idx_fixture(),
      asv_taxsort = asv_taxsort_id_fixture(),
      rank = "species",
      parent_rank = "genus",
      tax_ranks = tax_ranks_fixture
    ),
    "only `seq_idx` is not allowed"
  )
})

test_that("full_predenovo joins on seq_id", {
  out <- full_predenovo_taxon_table(
    closedref_taxon_table = closedref_taxon_id_fixture(),
    asv_taxsort = asv_taxsort_id_fixture(),
    rank = "species",
    parent_rank = "genus",
    tax_ranks = tax_ranks_fixture
  )
  expect_true("seq_idx" %in% names(out))
  expect_false("seq_idx_in" %in% names(out))
  expect_gt(nrow(out), 0L)
})

test_that("full_predenovo errors when asv_taxsort uses seq_id but taxon table has seq_idx only", {
  expect_error(
    full_predenovo_taxon_table(
      closedref_taxon_table = closedref_taxon_idx_fixture(),
      asv_taxsort = asv_taxsort_id_fixture(),
      rank = "species",
      parent_rank = "genus",
      tax_ranks = tax_ranks_fixture
    ),
    "only `seq_idx` is not allowed"
  )
})

test_that("validate_asv_taxsort rejects ambiguous or incomplete tables", {
  expect_error(
    optimotu.pipeline:::validate_asv_taxsort(
      tibble::tibble(seq_idx = 1L, seq_idx_in = 1L, seq_id = "a")
    ),
    "must not contain both"
  )
  expect_error(
    optimotu.pipeline:::validate_asv_taxsort(tibble::tibble(seq_idx = 1L)),
    "must contain exactly one"
  )
  expect_error(
    optimotu.pipeline:::validate_asv_taxsort(
      tibble::tibble(seq_idx_in = 1L)
    ),
    "must.include"
  )
})

small_table_asv_taxsort <- function(seq_id) {
  tibble::tibble(
    seq_id = seq_id,
    seq_idx = seq_along(seq_id),
    kingdom = "Euk",
    phylum = "Chordata",
    class = "Mammalia",
    order = "Primates",
    family = "Hominidae"
  )
}

expect_one_group_per_parent <- function(out, parent_rank) {
  groups_per_parent <- out |>
    dplyr::distinct(dplyr::pick(dplyr::all_of(c(parent_rank, "tar_group")))) |>
    dplyr::count(dplyr::pick(dplyr::all_of(parent_rank)))
  expect_true(all(groups_per_parent$n == 1L))
}

test_that("small_preclosed_taxon_table splits leftover taxa across tar_groups", {
  known <- tibble::tibble(
    seq_id = sprintf("ASV%04d", 1:8),
    genus = rep(c("Homo", "Pan"), each = 4L),
    species = c(
      "sapiens",
      "sapiens",
      NA,
      NA,
      "troglodytes",
      "troglodytes",
      NA,
      NA
    )
  )
  asv_taxsort <- small_table_asv_taxsort(known$seq_id)
  # Each genus: 2 unknown * 2 known = 4 ops.
  out_split <- small_preclosed_taxon_table(
    known_taxon_table = known,
    asv_taxsort = asv_taxsort,
    rank = "species",
    parent_rank = "genus",
    tax_ranks = tax_ranks_fixture,
    max_ops = 5
  )
  expect_equal(dplyr::n_distinct(out_split$tar_group), 2L)
  expect_one_group_per_parent(out_split, "genus")

  out_packed <- small_preclosed_taxon_table(
    known_taxon_table = known,
    asv_taxsort = asv_taxsort,
    rank = "species",
    parent_rank = "genus",
    tax_ranks = tax_ranks_fixture,
    max_ops = 8
  )
  expect_equal(dplyr::n_distinct(out_packed$tar_group), 1L)
})

test_that("small_predenovo_taxon_table splits leftover taxa across tar_groups", {
  closedref <- tibble::tibble(
    seq_id = sprintf("ASV%04d", 1:6),
    genus = rep(c("Homo", "Pan"), each = 3L),
    species = NA_character_
  )
  asv_taxsort <- small_table_asv_taxsort(closedref$seq_id)
  # Each genus: 3 * (3 - 1) / 2 = 3 ops.
  out_split <- small_predenovo_taxon_table(
    closedref_taxon_table = closedref,
    asv_taxsort = asv_taxsort,
    rank = "species",
    parent_rank = "genus",
    tax_ranks = tax_ranks_fixture,
    max_ops = 4
  )
  expect_equal(dplyr::n_distinct(out_split$tar_group), 2L)
  expect_one_group_per_parent(out_split, "genus")

  out_packed <- small_predenovo_taxon_table(
    closedref_taxon_table = closedref,
    asv_taxsort = asv_taxsort,
    rank = "species",
    parent_rank = "genus",
    tax_ranks = tax_ranks_fixture,
    max_ops = 6
  )
  expect_equal(dplyr::n_distinct(out_packed$tar_group), 1L)
})

test_that("preclosed large/small partition on min_ops, pack on max_ops", {
  known <- tibble::tibble(
    seq_id = sprintf("ASV%04d", 1:10),
    genus = c(rep("Homo", 4L), rep("Pan", 6L)),
    species = c(
      "sapiens",
      "sapiens",
      NA,
      NA,
      "troglodytes",
      "troglodytes",
      "troglodytes",
      NA,
      NA,
      NA
    )
  )
  asv_taxsort <- small_table_asv_taxsort(known$seq_id)
  # Homo: 2 unknown * 2 known = 4; Pan: 3 * 3 = 9.
  out_large <- large_preclosed_taxon_table(
    known_taxon_table = known,
    asv_taxsort = asv_taxsort,
    rank = "species",
    parent_rank = "genus",
    tax_ranks = tax_ranks_fixture,
    min_ops = 5
  )
  out_small <- small_preclosed_taxon_table(
    known_taxon_table = known,
    asv_taxsort = asv_taxsort,
    rank = "species",
    parent_rank = "genus",
    tax_ranks = tax_ranks_fixture,
    min_ops = 5,
    max_ops = 100
  )
  expect_equal(sort(unique(out_large$genus)), "Pan")
  expect_equal(sort(unique(out_small$genus)), "Homo")
  expect_equal(dplyr::n_distinct(out_small$tar_group), 1L)
})

test_that("predenovo large/small partition on min_ops, pack on max_ops", {
  closedref <- tibble::tibble(
    seq_id = sprintf("ASV%04d", 1:8),
    genus = c(rep("Homo", 3L), rep("Pan", 5L)),
    species = NA_character_
  )
  asv_taxsort <- small_table_asv_taxsort(closedref$seq_id)
  # Homo: 3 * 2 / 2 = 3; Pan: 5 * 4 / 2 = 10.
  out_large <- large_predenovo_taxon_table(
    closedref_taxon_table = closedref,
    asv_taxsort = asv_taxsort,
    rank = "species",
    parent_rank = "genus",
    tax_ranks = tax_ranks_fixture,
    min_ops = 5
  )
  out_small <- small_predenovo_taxon_table(
    closedref_taxon_table = closedref,
    asv_taxsort = asv_taxsort,
    rank = "species",
    parent_rank = "genus",
    tax_ranks = tax_ranks_fixture,
    min_ops = 5,
    max_ops = 100
  )
  expect_equal(sort(unique(out_large$genus)), "Pan")
  expect_equal(sort(unique(out_small$genus)), "Homo")
  expect_equal(dplyr::n_distinct(out_small$tar_group), 1L)
})

test_that("large_preclosed_taxon_table packs taxa across tar_groups", {
  known <- tibble::tibble(
    seq_id = sprintf("ASV%04d", 1:12),
    genus = rep(c("Homo", "Pan"), each = 6L),
    species = c(
      "sapiens",
      "sapiens",
      "sapiens",
      NA,
      NA,
      NA,
      "troglodytes",
      "troglodytes",
      "troglodytes",
      NA,
      NA,
      NA
    )
  )
  asv_taxsort <- small_table_asv_taxsort(known$seq_id)
  # Each genus: 3 unknown * 3 known = 9 ops.
  out_split <- large_preclosed_taxon_table(
    known_taxon_table = known,
    asv_taxsort = asv_taxsort,
    rank = "species",
    parent_rank = "genus",
    tax_ranks = tax_ranks_fixture,
    min_ops = 5,
    max_ops = 10
  )
  expect_equal(dplyr::n_distinct(out_split$tar_group), 2L)
  expect_one_group_per_parent(out_split, "genus")

  out_packed <- large_preclosed_taxon_table(
    known_taxon_table = known,
    asv_taxsort = asv_taxsort,
    rank = "species",
    parent_rank = "genus",
    tax_ranks = tax_ranks_fixture,
    min_ops = 5,
    max_ops = 20
  )
  expect_equal(dplyr::n_distinct(out_packed$tar_group), 1L)
})

test_that("large_predenovo_taxon_table packs taxa across tar_groups", {
  closedref <- tibble::tibble(
    seq_id = sprintf("ASV%04d", 1:10),
    genus = rep(c("Homo", "Pan"), each = 5L),
    species = NA_character_
  )
  asv_taxsort <- small_table_asv_taxsort(closedref$seq_id)
  # Each genus: 5 * 4 / 2 = 10 ops.
  out_split <- large_predenovo_taxon_table(
    closedref_taxon_table = closedref,
    asv_taxsort = asv_taxsort,
    rank = "species",
    parent_rank = "genus",
    tax_ranks = tax_ranks_fixture,
    min_ops = 5,
    max_ops = 12
  )
  expect_equal(dplyr::n_distinct(out_split$tar_group), 2L)
  expect_one_group_per_parent(out_split, "genus")

  out_packed <- large_predenovo_taxon_table(
    closedref_taxon_table = closedref,
    asv_taxsort = asv_taxsort,
    rank = "species",
    parent_rank = "genus",
    tax_ranks = tax_ranks_fixture,
    min_ops = 5,
    max_ops = 20
  )
  expect_equal(dplyr::n_distinct(out_packed$tar_group), 1L)
})

test_that("large_preclosed and small_predenovo delegate seq_id join path", {
  out_large <- large_preclosed_taxon_table(
    known_taxon_table = known_taxon_id_fixture(),
    asv_taxsort = asv_taxsort_id_fixture(),
    rank = "species",
    parent_rank = "genus",
    tax_ranks = tax_ranks_fixture,
    min_ops = 1
  )
  expect_true("seq_idx" %in% names(out_large))

  out_small <- small_predenovo_taxon_table(
    closedref_taxon_table = closedref_taxon_id_fixture(),
    asv_taxsort = asv_taxsort_id_fixture(),
    rank = "species",
    parent_rank = "genus",
    tax_ranks = tax_ranks_fixture,
    max_ops = 1e6
  )
  expect_true("seq_idx" %in% names(out_small))
})
