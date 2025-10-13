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
lulu_result <- lulu::lulu(lulu_otutab, lulu_matchlist)
sink()

# convert lulu OTU table to long format with seq_id as a factor
long_otutab <- tibble::as_tibble(lulu_otutab, rownames = "seq_id") |>
  tidyr::pivot_longer(-seq_id, names_to = "sample_key", values_to = "nread") |>
  dplyr::filter(nread > 0L) |>
  dplyr::mutate(n_sample = dplyr::n(), nread_tot = sum(nread), .by = seq_id) |>
  dplyr::arrange(dplyr::desc(n_sample), dplyr::desc(nread_tot)) |>
  dplyr::mutate(seq_id = ordered(seq_id, levels = unique(seq_id)))

# convert sequence IDs in the match list to factors with the same levels as in
# the OTU table
factor_matchlist <- lulu_matchlist |>
  dplyr::mutate(
    V1 = ordered(V1, levels = levels(long_otutab$seq_id)),
    V2 = ordered(V2, levels = levels(long_otutab$seq_id))
  )

# generate match table as it would have been calculated by running BLAST on
# each sample separately
long_matchlist <-
  dplyr::left_join(long_otutab, factor_matchlist, by = c("seq_id" = "V1"), relationship = "many-to-many") |>
  dplyr::rename(seq_id1 = seq_id, seq_id2 = V2) |>
  dplyr::inner_join(long_otutab, by = c("sample_key", "seq_id2" = "seq_id")) |>
  dplyr::filter(seq_id1 > seq_id2) |>
  dplyr::rename(nread1 = nread.x, nread2 = nread.y) |>
  dplyr::mutate(dist = optimotu::threshold_as_dist(V3))

# note: although "long_matchlist" contains the same sequence pairs multiple
# times when they co-occur more than once, it does not include pairs that never
# co-occur.  lulu_matchlist has 320k entries, while long_matchlist has only
# 60k. (long_matchlist also does not include the 2425 self-self matches and does
# not double-count pairs, but that would still leave 150k entries in the lulu
# version)

long_lulu_map <- optimotu.pipeline:::lulu_map_impl(
  match_id1 = long_matchlist$seq_id1,
  match_id2 = long_matchlist$seq_id2,
  match_nread1 = long_matchlist$nread1,
  match_nread2 = long_matchlist$nread2,
  match_dist = long_matchlist$dist,
  seq_idx = long_otutab$seq_id,
  nread = long_otutab$nread,
  max_dist = 0.16
) |>
  dplyr::mutate(
    seq_id = ordered(seq_idx, labels = levels(long_otutab$seq_id)),
    lulu_id = ordered(lulu_idx, levels = seq_idx, labels = levels(seq_id))
  )

long_lulu_otutab <- dplyr::left_join(long_otutab, long_lulu_map, by = "seq_id") |>
  dplyr::summarize(nread = sum(nread), .by = c(seq_id, sample_key)) |>
  tidyr::pivot_wider(names_from = "sample_key", values_from = "nread", values_fill = 0L)


