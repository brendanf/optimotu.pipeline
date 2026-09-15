# Three-level scoped LULU (batch / seqrun / global) equivalence tests
#
# Synthetic layout: 2 seqruns x 2 batches. OTUs:
#   1  parent in multiple seqruns
#   2  parent in one seqrun (both batches)
#   3  parent in one batch (seqrun 1, batch 1)
#   4  child of 3 in the same batch
#   5  child of 2 in the same seqrun
#   6  child of 1 across seqruns
#   7  singleton child of 3 in one batch
#   8  OTU in one batch of seqrun 2

make_scoped_fixture <- function() {
  otu <- list(
    # seqrun 1, batch 1
    data.frame(
      sample = c("s1", "s1", "s1", "s1", "s1", "s1", "s1"),
      seq_idx = c(1L, 2L, 3L, 4L, 5L, 6L, 7L),
      nread = c(100L, 80L, 60L, 20L, 25L, 30L, 5L)
    ),
    # seqrun 1, batch 2
    data.frame(
      sample = c("s2", "s2", "s2", "s2"),
      seq_idx = c(1L, 2L, 5L, 6L),
      nread = c(90L, 70L, 22L, 28L)
    ),
    # seqrun 2, batch 1
    data.frame(
      sample = c("s3", "s3", "s3"),
      seq_idx = c(1L, 6L, 8L),
      nread = c(95L, 27L, 10L)
    ),
    # seqrun 2, batch 2
    data.frame(
      sample = c("s4", "s4"),
      seq_idx = c(1L, 6L),
      nread = c(85L, 26L)
    )
  )
  mk <- function(i1, i2, n1, n2, dist = 0.05) {
    data.frame(
      seq_idx1 = as.integer(i1),
      seq_idx2 = as.integer(i2),
      nread1 = as.integer(n1),
      nread2 = as.integer(n2),
      dist = dist
    )
  }
  matches <- list(
    dplyr::bind_rows(
      mk(3, 4, 60, 20),
      mk(3, 7, 60, 5),
      mk(2, 5, 80, 25),
      mk(1, 6, 100, 30),
      mk(2, 4, 80, 20),
      mk(1, 4, 100, 20),
      mk(1, 5, 100, 25),
      mk(1, 2, 100, 80),
      mk(1, 3, 100, 60),
      mk(2, 3, 80, 60),
      mk(3, 5, 60, 25),
      mk(3, 6, 60, 30),
      mk(2, 6, 80, 30),
      mk(2, 7, 80, 5),
      mk(1, 7, 100, 5),
      mk(4, 5, 20, 25),
      mk(4, 6, 20, 30),
      mk(5, 6, 25, 30),
      mk(4, 7, 20, 5),
      mk(5, 7, 25, 5),
      mk(6, 7, 30, 5)
    ),
    dplyr::bind_rows(
      mk(2, 5, 70, 22),
      mk(1, 6, 90, 28),
      mk(1, 5, 90, 22),
      mk(1, 2, 90, 70),
      mk(2, 6, 70, 28),
      mk(5, 6, 22, 28)
    ),
    dplyr::bind_rows(
      mk(1, 6, 95, 27),
      mk(1, 8, 95, 10),
      mk(6, 8, 27, 10)
    ),
    dplyr::bind_rows(
      mk(1, 6, 85, 26)
    )
  )
  list(otu = otu, matches = matches, seqrun_ids = c(1L, 1L, 2L, 2L))
}

run_scoped_split <- function(
  otu,
  matches,
  seqrun_ids,
  max_dist = 0.1,
  min_abundance_ratio = 1,
  min_cooccurrence_ratio = 1,
  use_mean_abundance_ratio = FALSE
) {
  stats <- optimotu.pipeline:::lulu_otu_stats_dfs_impl(otu, seqrun_ids)
  batch_maps <- lapply(seq_along(matches), function(i) {
    optimotu.pipeline:::lulu_map_scoped_dfs_impl(
      stats,
      list(matches[[i]]),
      "batch",
      max_dist,
      min_abundance_ratio,
      min_cooccurrence_ratio,
      use_mean_abundance_ratio
    )
  })
  seqrun_maps <- lapply(unique(seqrun_ids), function(sid) {
    idx <- which(seqrun_ids == sid)
    optimotu.pipeline:::lulu_map_scoped_dfs_impl(
      stats,
      matches[idx],
      "seqrun",
      max_dist,
      min_abundance_ratio,
      min_cooccurrence_ratio,
      use_mean_abundance_ratio
    )
  })
  final_map <- optimotu.pipeline:::lulu_map_scoped_dfs_impl(
    stats,
    matches,
    "global",
    max_dist,
    min_abundance_ratio,
    min_cooccurrence_ratio,
    use_mean_abundance_ratio
  )
  combined <- optimotu.pipeline:::lulu_map_combine_impl(
    stats,
    c(batch_maps, seqrun_maps, list(final_map))
  )
  list(stats = stats, map = combined)
}

run_reference_map <- function(
  otu,
  matches,
  max_dist = 0.1,
  min_abundance_ratio = 1,
  min_cooccurrence_ratio = 1,
  use_mean_abundance_ratio = FALSE
) {
  optimotu.pipeline::lulu_map(
    otu_table = dplyr::bind_rows(otu),
    match_table = dplyr::bind_rows(matches),
    max_dist = max_dist,
    min_abundance_ratio = min_abundance_ratio,
    min_cooccurrence_ratio = min_cooccurrence_ratio,
    use_mean_abundance_ratio = use_mean_abundance_ratio,
    id_is_int = TRUE,
    id_is_sorted = FALSE
  )
}

testthat::test_that("scoped LULU grains classify OTUs by OTU-table partition", {
  fx <- make_scoped_fixture()
  stats <- optimotu.pipeline:::lulu_otu_stats_dfs_impl(fx$otu, fx$seqrun_ids)
  grain <- stats$grain[match(1:8, stats$seq_idx)]
  # 1 multi-seqrun, 2 seqrun, 3 batch, 4 batch, 5 seqrun, 6 multi-seqrun,
  # 7 batch, 8 batch (only in seqrun2 batch1 in this fixture)
  testthat::expect_equal(grain, c(2L, 1L, 0L, 0L, 1L, 2L, 0L, 0L))
  testthat::expect_equal(stats$occurrence[stats$seq_idx == 7L], 1L)
})

testthat::test_that("scoped split+combine matches lulu_map at ratio 1", {
  fx <- make_scoped_fixture()
  split <- run_scoped_split(
    fx$otu,
    fx$matches,
    fx$seqrun_ids,
    min_cooccurrence_ratio = 1
  )
  ref <- run_reference_map(
    fx$otu,
    fx$matches,
    min_cooccurrence_ratio = 1
  )
  split$map <- split$map[order(split$map$seq_idx), ]
  ref <- ref[order(ref$seq_idx), ]
  testthat::expect_equal(split$map$seq_idx, ref$seq_idx)
  testthat::expect_equal(split$map$lulu_idx, ref$lulu_idx)
})

testthat::test_that("scoped split+combine matches lulu_map at ratio 0.5", {
  fx <- make_scoped_fixture()
  # Parent in a finer partition: child 8 also co-occurs with batch parent 3
  # in one sample. Child 8 occurs twice, so co-occurrence ratio with 3 is 0.5.
  fx$otu[[1]] <- dplyr::bind_rows(
    fx$otu[[1]],
    data.frame(sample = "s1b", seq_idx = 8L, nread = 8L),
    data.frame(sample = "s1b", seq_idx = 3L, nread = 50L)
  )
  fx$matches[[1]] <- dplyr::bind_rows(
    fx$matches[[1]],
    data.frame(
      seq_idx1 = 3L,
      seq_idx2 = 8L,
      nread1 = 50L,
      nread2 = 8L,
      dist = 0.05
    )
  )

  split_lo <- run_scoped_split(
    fx$otu,
    fx$matches,
    fx$seqrun_ids,
    min_cooccurrence_ratio = 0.5
  )
  ref_lo <- run_reference_map(
    fx$otu,
    fx$matches,
    min_cooccurrence_ratio = 0.5
  )
  split_lo$map <- split_lo$map[order(split_lo$map$seq_idx), ]
  ref_lo <- ref_lo[order(ref_lo$seq_idx), ]
  testthat::expect_equal(split_lo$map$seq_idx, ref_lo$seq_idx)
  testthat::expect_equal(split_lo$map$lulu_idx, ref_lo$lulu_idx)

  split_1 <- run_scoped_split(
    fx$otu,
    fx$matches,
    fx$seqrun_ids,
    min_cooccurrence_ratio = 1
  )
  ref_1 <- run_reference_map(
    fx$otu,
    fx$matches,
    min_cooccurrence_ratio = 1
  )
  split_1$map <- split_1$map[order(split_1$map$seq_idx), ]
  ref_1 <- ref_1[order(ref_1$seq_idx), ]
  testthat::expect_equal(split_1$map$lulu_idx, ref_1$lulu_idx)
})

testthat::test_that("singleton fast-path matches reference lulu_map", {
  fx <- make_scoped_fixture()
  split <- run_scoped_split(fx$otu, fx$matches, fx$seqrun_ids)
  row7 <- split$map[split$map$seq_idx == 7L, ]
  testthat::expect_equal(nrow(row7), 1L)
  testthat::expect_true(row7$lulu_idx != 7L)
  ref <- run_reference_map(fx$otu, fx$matches)
  testthat::expect_equal(
    split$map$lulu_idx[split$map$seq_idx == 7L],
    ref$lulu_idx[ref$seq_idx == 7L]
  )
})

testthat::test_that("lulu_map_combine path-compresses cross-grain chains", {
  stats <- data.frame(
    seq_idx = 1:4,
    occurrence = c(4L, 3L, 2L, 1L),
    abundance = c(400L, 300L, 200L, 100L),
    grain = c(2L, 1L, 0L, 0L)
  )
  combined <- optimotu.pipeline:::lulu_map_combine_impl(
    stats,
    list(
      data.frame(seq_idx = 4L, lulu_idx = 3L),
      data.frame(seq_idx = 3L, lulu_idx = 2L),
      data.frame(seq_idx = 2L, lulu_idx = 1L)
    )
  )
  testthat::expect_equal(
    combined$lulu_idx[match(1:4, combined$seq_idx)],
    c(1L, 1L, 1L, 1L)
  )
})
