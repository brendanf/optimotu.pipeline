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
