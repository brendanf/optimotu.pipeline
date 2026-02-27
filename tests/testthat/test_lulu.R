# Test data for LULU:
# There are 3 true "parent" sequences, which occur in 0.9, 0.5, and 0.1
# of 1000 samples.
# Conditional on presence, read count is Poisson(1000).
# Each has 3 noisy "child" sequences.
# 1. Occurs in the same sample as parent with probability 0.99, and other
# samples with probability 0.01 * parent probability.  Half read count of parent when present,
# otherwise Poisson(100).
# 2. Occurs in the same sample as parent with probability 0.5, and half read
# count of parent.
# 3. Occurs in the same sample as parent with probability 0.2. Read count
# Normal(0.9 * parent, 0.1).

set.seed(123)
comm_matrix <- matrix(
  0L,
  nrow = 1000,
  ncol = 12,
  dimnames = list(
    1:1000,
    c(
      "p1", "p2", "p3",
      "c1.1", "c2.1", "c3.1",
      "c1.2", "c2.2", "c3.2",
      "c1.3", "c2.3", "c3.3"
    )
  )
)

comm_matrix[, 1] <- rbinom(1000L, 1, 0.9) * rpois(1000, 1000) # 908
comm_matrix[, 2] <- rbinom(1000L, 1, 0.5) * rpois(1000, 2000) # 497
comm_matrix[, 3] <- rbinom(1000L, 1, 0.1) * rpois(1000, 4000) # 93

comm_matrix[, 4] <- ifelse(
  comm_matrix[, 1],
  rbinom(1000L, 1, 0.99) * comm_matrix[, 1] %/% 2,
  rbinom(1000, 1, 0.09) * rpois(1000, 100)
) # 904
comm_matrix[, 5] <- ifelse(
  comm_matrix[, 2],
  rbinom(1000L, 1, 0.99) * comm_matrix[, 2] %/% 2,
  rbinom(1000, 1, 0.01) * rpois(1000, 200)
) # 496
comm_matrix[, 6] <- ifelse(
  comm_matrix[, 3],
  rbinom(1000L, 1, 0.99) * comm_matrix[, 3] %/% 2,
  rbinom(1000, 1, 0.002) * rpois(1000, 400)
) # 91

comm_matrix[, 7] <- rbinom(1000, 1, 0.5) * comm_matrix[, 1] %/% 2 # 434
comm_matrix[, 8] <- rbinom(1000, 1, 0.5) * comm_matrix[, 2] %/% 2 # 253
comm_matrix[, 9] <- rbinom(1000, 1, 0.5) * comm_matrix[, 3] %/% 2 # 47

comm_matrix[, 10] <- as.integer(
  comm_matrix[, 1] * rbinom(1000, 1, 0.2) * rnorm(1000, 0.9, 0.05)
) # 179
comm_matrix[, 11] <- as.integer(
  comm_matrix[, 2] * rbinom(1000, 1, 0.2) * rnorm(1000, 0.9, 0.05)
) # 99
comm_matrix[, 12] <- as.integer(
  comm_matrix[, 3] * rbinom(1000, 1, 0.2) * rnorm(1000, 0.9, 0.05)
) # 18

# Reordered, the indices of the sequences are 1,3,9, 2,4,10, 5,6,11, 7,8,12
comm_matrix <-
  comm_matrix[, order(colSums(sign(comm_matrix)), decreasing = TRUE)]

comm_matrix_long <- tibble::as_tibble(comm_matrix) |>
  tibble::rowid_to_column("sample") |>
  tidyr::pivot_longer(-sample, names_to = "seq_id", values_to = "nread") |>
  dplyr::mutate(seq_id = ordered(seq_id, levels = colnames(comm_matrix))) |>
  dplyr::filter(nread > 0L)

pairs <- dplyr::inner_join(
  comm_matrix_long,
  comm_matrix_long,
  by = "sample",
  suffix = c("1", "2"),
  relationship = "many-to-many"
) |>
  dplyr::filter(seq_id1 < seq_id2) |>
  dplyr::mutate(dist = 0.05)

testthat::test_that("lulu_map_impl works", {
  testthat::expect_equal(
    optimotu.pipeline:::lulu_map_impl(
      pairs$seq_id1,
      pairs$seq_id2,
      pairs$nread1,
      pairs$nread2,
      pairs$dist,
      comm_matrix_long$seq_id,
      comm_matrix_long$nread,
      0.16,
      min_cooccurrence_ratio = 0.95,
      min_abundance_ratio = 0.9,
      use_mean_abundance_ratio = TRUE
    ),
    tibble::tibble(
      seq_idx = 1L:12L,
      lulu_idx = c(1L, 3L, 9L)[as.integer(substring(colnames(comm_matrix), 2L))]
    ) |>
      as.data.frame()
  )
})

# Run example from LULU github

download.file(
  "https://raw.githubusercontent.com/tobiasgf/lulu/master/Example_data/centroids_test.txt",
  destfile = "centroids_test.txt"
)
download.file(
  "https://raw.githubusercontent.com/tobiasgf/lulu/master/Example_data/otutable_test.txt",
  destfile = "otutable_test.txt"
)

# run in terminal
# blastn -query centroids_test.txt\
#        -subject centroids_test.txt\
#        -outfmt "6 qseqid sseqid pident"\
#        -out blast_output.txt\
#        -qcov_hsp_perc 80\
#        -perc_identity 84

lulu_otutab <- read.csv("otutable_test.txt", sep = "\t", header = TRUE, as.is = TRUE, row.names = 1)
lulu_matchlist <- read.table("blast_output.txt", header = FALSE, as.is = TRUE, stringsAsFactors = FALSE)

# run lulu to get reference result
# lulu is very verbose...
sink("/dev/null")
# use minimum_relative_cooccurrence = 1.0 to avoid LULU bug
lulu_result <- lulu::lulu(lulu_otutab, lulu_matchlist, minimum_ratio_type = "min", minimum_relative_cooccurence = 1.0)
sink()

# convert lulu OTU table to long format
long_otutab <- tibble::as_tibble(lulu_otutab, rownames = "seq_id") |>
  tidyr::pivot_longer(-seq_id, names_to = "sample_key", values_to = "nread") |>
  dplyr::filter(nread > 0L)

# generate match table as it would have been calculated by running BLAST on
# each sample separately
long_matchlist <-
  # remove self-hits
  dplyr::filter(lulu_matchlist, V1 != V2) |>
  # swap hits where the second hit is alphabetically earlier, in order to
  # deduplicate the hit list.
  # The optimotu.pipeline implementation does not care which sequence is seq1
  # and which is seq2, but it does need there to be only one hit per pair per sample
  dplyr::mutate(
    V4 = pmin(V1, V2),
    V2 = pmax(V1, V2),
    V1 = V4
  ) |>
  # The example data give slightly different values for hits when query and
  # subject are swapped.  Take the larger similarity/smaller distance.
  dplyr::summarize(V3 = max(V3), .by = c(V1, V2)) |>
  # join with long OTU table to get read counts for each sample for seq1
  dplyr::left_join(long_otutab, y = _, by = c("seq_id" = "V1"), relationship = "many-to-many") |>
  dplyr::rename(seq_id1 = seq_id, seq_id2 = V2) |>
  # join again to get read counts for seq2
  dplyr::inner_join(long_otutab, by = c("sample_key", "seq_id2" = "seq_id")) |>
  dplyr::rename(nread1 = nread.x, nread2 = nread.y) |>
  # convert to fractional distance
  dplyr::mutate(dist = optimotu::threshold_as_dist(V3)) |>
  dplyr::select(seq_id1, seq_id2, nread1, nread2, dist)

# Note: although "long_matchlist" contains the same sequence pairs multiple
# times when they co-occur more than once, it does not include pairs that never
# co-occur.  lulu_matchlist has 320k entries, while long_matchlist has only
# 60k. This illustrates the fact that even though work is being redone,
# searching for matching pairs within each sample combining results takes less
# computation than doing a global all-vs-all search.

long_lulu_map <- optimotu.pipeline::lulu_map(
  otu_table = long_otutab,
  match_table = long_matchlist,
  max_dist = 0.16,
  min_abundance_ratio = 1,
  min_cooccurrence_ratio = 1,
  use_mean_abundance_ratio = FALSE,
  verbose = 0
)

testthat::test_that("lulu_map and lulu give same map result", {
  testthat::expect_equal(
    long_lulu_map$lulu_id,
    lulu_result$otu_map[long_lulu_map$seq_id, "parent_id"]
  )
})

long_lulu_otutab <- optimotu.pipeline::lulu_table(long_lulu_map, long_otutab)

# pivot to a wide table to match LULU result
wide_lulu_otutab <- tidyr::pivot_wider(
  long_lulu_otutab,
  names_from = sample_key,
  values_from = nread,
  values_fill = 0L
) |>
  tibble::column_to_rownames("seq_id")

# LULU keeps empty rows (OTUs). optimotu does not.  Check that this is accurate
testthat::test_that("lulu_table does not include rows which are empty in lulu", {
  testthat::expect_length(intersect(rownames(wide_lulu_otutab), rownames(lulu_result$curated_table)[rowSums(lulu_result$curated_table) == 0]), 0)
})

# Add these rows
wide_lulu_otutab[rownames(lulu_result$curated_table)[rowSums(lulu_result$curated_table) == 0], ] <- 0L

# reorder the table to match the LULU result
wide_lulu_otutab <- wide_lulu_otutab[rownames(lulu_result$curated_table), colnames(lulu_result$curated_table)]

testthat::test_that("lulu_table and lulu give same result", {
  testthat::expect_equal(
    wide_lulu_otutab,
    lulu_result$curated_table
  )
})
