skip_if_no_vsearch <- function() {
  tc_skip_if_no_exec("vsearch")
}

write_qual_fastq <- function(seqs, file, qual = 40L) {
  q <- lapply(Biostrings::width(seqs), \(n) as.integer(rep(qual, n))) |>
    S4Vectors::List() |>
    Biostrings::PhredQuality()
  Biostrings::writeQualityScaledXStringSet(
    Biostrings::QualityScaledDNAStringSet(seqs, q),
    file
  )
  file
}

# Overlapping pairs: 40 bp unique prefix/suffix + 80 bp overlap.
# Sequences are drawn so they are not palindromic/repetitive, which would
# give vsearch multiple possible merge alignments.
make_overlapping_pairs <- function(n, seed, ids) {
  withr::with_seed(seed, {
    overlap <- paste(
      sample(c("A", "C", "G", "T"), 80L, replace = TRUE),
      collapse = ""
    )
    prefix <- paste(
      sample(c("A", "C", "G", "T"), 40L, replace = TRUE),
      collapse = ""
    )
    suffix <- paste(
      sample(c("A", "C", "G", "T"), 40L, replace = TRUE),
      collapse = ""
    )
  })
  fwd <- Biostrings::DNAStringSet(rep(paste0(prefix, overlap), n))
  rev <- Biostrings::reverseComplement(
    Biostrings::DNAStringSet(rep(paste0(overlap, suffix), n))
  )
  names(fwd) <- ids
  names(rev) <- ids
  list(fwd = fwd, rev = rev, merged = paste0(prefix, overlap, suffix))
}

test_that("merged_filter_options constructs and updates", {
  opts <- merged_filter_options(maxEE = 2, minLen = 100L)
  expect_s3_class(opts, "merged_filter_options")
  expect_equal(opts$maxEE, 2)
  expect_equal(opts$minLen, 100L)
  updated <- stats::update(opts, list(maxEE = 0.5, maxNs = 0L))
  expect_equal(updated$maxEE, 0.5)
  expect_equal(updated$maxNs, 0L)
  expect_equal(updated$minLen, 100L)
})

test_that("vsearch_fastq_merge_pairs handles empty input", {
  expect_equal(
    vsearch_fastq_merge_pairs(
      seq_R1 = character(),
      seq_R2 = character(),
      seq_out = character()
    ),
    character()
  )
})

test_that("vsearch_fastq_merge_pairs skips empty files and clamps extra shards", {
  skip_if_no_vsearch()
  empty_r1 <- withr::local_tempfile(fileext = ".fastq.gz")
  empty_r2 <- withr::local_tempfile(fileext = ".fastq.gz")
  empty_out <- withr::local_tempfile(fileext = ".fastq.gz")
  close(gzfile(empty_r1, "wb"))
  close(gzfile(empty_r2, "wb"))

  expect_equal(
    vsearch_fastq_merge_pairs(
      empty_r1,
      empty_r2,
      empty_out,
      threads = 1,
      shards = 8
    ),
    empty_out
  )
  expect_equal(sequence_size(empty_out), 0L)

  pairs <- make_overlapping_pairs(3L, 1L, sprintf("read%02d", 1:3))
  r1 <- withr::local_tempfile(fileext = ".fastq")
  r2 <- withr::local_tempfile(fileext = ".fastq")
  out <- withr::local_tempfile(fileext = ".fastq")
  write_qual_fastq(pairs$fwd, r1)
  write_qual_fastq(pairs$rev, r2)
  expect_equal(
    vsearch_fastq_merge_pairs(
      r1,
      r2,
      out,
      min_overlap = 16,
      max_mismatch = 2,
      filter_options = merged_filter_options(maxEE = 1),
      threads = 1,
      shards = 8
    ),
    out
  )
  expect_equal(sequence_size(out), 3L)
})

test_that("vsearch_fastq_merge_pairs merges overlapping reads and respects minLen", {
  skip_if_no_vsearch()
  pairs <- make_overlapping_pairs(5L, 1L, sprintf("read%02d", 1:5))
  r1 <- withr::local_tempfile(fileext = ".fastq")
  r2 <- withr::local_tempfile(fileext = ".fastq")
  out <- withr::local_tempfile(fileext = ".fastq")
  write_qual_fastq(pairs$fwd, r1)
  write_qual_fastq(pairs$rev, r2)

  result <- vsearch_fastq_merge_pairs(
    r1,
    r2,
    out,
    min_overlap = 16,
    max_mismatch = 2,
    filter_options = merged_filter_options(maxEE = 1),
    threads = 1,
    shards = 1
  )
  expect_equal(result, out)
  expect_equal(sequence_size(out), 5L)

  out_filt <- withr::local_tempfile(fileext = ".fastq")
  vsearch_fastq_merge_pairs(
    r1,
    r2,
    out_filt,
    min_overlap = 16,
    filter_options = merged_filter_options(maxEE = 1, minLen = 10000L),
    threads = 1,
    shards = 1
  )
  expect_equal(sequence_size(out_filt), 0L)

  out_gz <- withr::local_tempfile(fileext = ".fastq.gz")
  vsearch_fastq_merge_pairs(
    r1,
    r2,
    out_gz,
    min_overlap = 16,
    max_mismatch = 2,
    filter_options = merged_filter_options(maxEE = 1),
    threads = 1,
    shards = 1
  )
  expect_equal(sequence_size(out_gz), 5L)
})

test_that("vsearch_cluster_unoise2 returns empty uc_cluster for empty files", {
  skip_if_no_vsearch()
  empty <- withr::local_tempfile(fileext = ".fasta")
  writeLines(character(), empty)
  out <- vsearch_cluster_unoise2(empty, min_size = 8, threads = 1, shards = 1)
  expect_named(out, empty)
  expect_s3_class(out[[1]], "uc_cluster")
  expect_equal(nrow(out[[1]]$clusters), 0L)
  expect_equal(
    vsearch_cluster_unoise2(character()),
    stats::setNames(list(), character())
  )
})

test_that("UNOISE clustering, seqtable, and read_map agree on synthetic ASVs", {
  skip_if_no_vsearch()
  n_a <- 12L
  n_b <- 10L
  n_c <- 2L
  ids_a <- sprintf("A%02d", seq_len(n_a))
  ids_b <- sprintf("B%02d", seq_len(n_b))
  ids_c <- sprintf("C%02d", seq_len(n_c))
  pairs_a <- make_overlapping_pairs(n_a, 11L, ids_a)
  pairs_b <- make_overlapping_pairs(n_b, 22L, ids_b)
  pairs_c <- make_overlapping_pairs(n_c, 33L, ids_c)

  r1 <- withr::local_tempfile(fileext = ".fastq")
  r2 <- withr::local_tempfile(fileext = ".fastq")
  merged <- withr::local_tempfile(fileext = ".fastq")
  raw <- withr::local_tempfile(fileext = ".fastq")
  trim <- withr::local_tempfile(fileext = ".fastq")
  fwd <- c(pairs_a$fwd, pairs_b$fwd, pairs_c$fwd)
  rev <- c(pairs_a$rev, pairs_b$rev, pairs_c$rev)
  write_qual_fastq(fwd, r1)
  write_qual_fastq(rev, r2)
  write_qual_fastq(fwd, raw)
  write_qual_fastq(fwd, trim)

  vsearch_fastq_merge_pairs(
    r1,
    r2,
    merged,
    min_overlap = 16,
    filter_options = merged_filter_options(maxEE = 1),
    threads = 1,
    shards = 1
  )
  expect_equal(sequence_size(merged), n_a + n_b + n_c)

  uc <- vsearch_cluster_unoise2(
    merged,
    min_size = 8,
    alpha = 2,
    threads = 1,
    shards = 1
  )[[1]]
  expect_s3_class(uc, "uc_cluster")
  expect_true(nrow(uc$clusters) >= 2L)
  centroids <- uc$clusters$seq
  expect_true(pairs_a$merged %in% centroids)
  expect_true(pairs_b$merged %in% centroids)
  expect_false(pairs_c$merged %in% centroids)

  seq_all <- Biostrings::DNAStringSet(unique(c(centroids, pairs_c$merged)))
  names(seq_all) <- seq_along(seq_all)
  dm <- make_denoise_map(stats::setNames(list(uc), "sample1"), seq_all)
  expect_true(all(dm$denoise_idx %in% uc$clusters$clust_idx))
  seqtable <- denoise_map_to_seqtable(dm)
  expect_named(seqtable, c("sample", "seq_idx", "nread"))
  expect_equal(seqtable$sample, rep("sample1", nrow(seqtable)))
  idx_a <- as.integer(BiocGenerics::match(pairs_a$merged, seq_all))
  idx_b <- as.integer(BiocGenerics::match(pairs_b$merged, seq_all))
  nread_a <- seqtable$nread[seqtable$seq_idx == idx_a]
  nread_b <- seqtable$nread[seqtable$seq_idx == idx_b]
  expect_equal(as.integer(nread_a), n_a)
  expect_equal(as.integer(nread_b), n_b)

  smap <- unoise_read_map(
    sample = "sample1",
    fq_raw = raw,
    fq_trim = trim,
    fq_merged = merged,
    uc = uc,
    denoise_map = dm
  )
  expect_equal(nrow(smap), n_a + n_b + n_c)
  expect_equal(smap$sample, rep("sample1", nrow(smap)))
  expect_true(all(as.integer(smap$flags) >= 0x03))
  denoised <- bitwAnd(as.integer(smap$flags), 0x04) != 0L
  expect_equal(sum(denoised), n_a + n_b)
  expect_true(all(is.na(smap$seq_idx[!denoised])))
})

test_that("unoise_read_map vectorizes over samples and matches scalar calls", {
  skip_if_no_vsearch()
  make_sample <- function(n, seed, label) {
    ids <- sprintf("%s%02d", label, seq_len(n))
    pairs <- make_overlapping_pairs(n, seed, ids)
    # Temp files must outlive this helper; create them in the test env.
    r1 <- tempfile(fileext = ".fastq")
    r2 <- tempfile(fileext = ".fastq")
    merged <- tempfile(fileext = ".fastq")
    raw <- tempfile(fileext = ".fastq")
    trim <- tempfile(fileext = ".fastq")
    withr::defer(
      unlink(c(r1, r2, merged, raw, trim)),
      envir = parent.frame(1L)
    )
    write_qual_fastq(pairs$fwd, r1)
    write_qual_fastq(pairs$rev, r2)
    write_qual_fastq(pairs$fwd, raw)
    write_qual_fastq(pairs$fwd, trim)
    vsearch_fastq_merge_pairs(
      r1,
      r2,
      merged,
      min_overlap = 16,
      filter_options = merged_filter_options(maxEE = 1),
      threads = 1,
      shards = 1
    )
    uc <- vsearch_cluster_unoise2(
      merged,
      min_size = 8,
      alpha = 2,
      threads = 1,
      shards = 1
    )[[1]]
    list(
      sample = label,
      raw = raw,
      trim = trim,
      merged = merged,
      uc = uc,
      merged_seq = pairs$merged
    )
  }
  s1 <- make_sample(12L, 11L, "A")
  s2 <- make_sample(10L, 22L, "B")
  seq_all <- Biostrings::DNAStringSet(unique(c(
    s1$uc$clusters$seq,
    s2$uc$clusters$seq
  )))
  names(seq_all) <- seq_along(seq_all)
  uc_list <- stats::setNames(list(s1$uc, s2$uc), c(s1$sample, s2$sample))
  dm <- make_denoise_map(uc_list, seq_all)

  scalar <- dplyr::bind_rows(
    unoise_read_map(
      s1$sample,
      s1$raw,
      s1$trim,
      s1$merged,
      s1$uc,
      dm
    ),
    unoise_read_map(
      s2$sample,
      s2$raw,
      s2$trim,
      s2$merged,
      s2$uc,
      dm
    )
  )
  vectorized <- unoise_read_map(
    sample = c(s1$sample, s2$sample),
    fq_raw = c(s1$raw, s2$raw),
    fq_trim = c(s1$trim, s2$trim),
    fq_merged = c(s1$merged, s2$merged),
    uc = list(s1$uc, s2$uc),
    denoise_map = dm
  )
  expect_equal(vectorized, scalar)
})

test_that("unoise_read_map handles empty merged FASTQ", {
  skip_if_no_vsearch()
  empty <- withr::local_tempfile(fileext = ".fastq")
  writeLines(character(), empty)
  raw <- withr::local_tempfile(fileext = ".fastq")
  trim <- withr::local_tempfile(fileext = ".fastq")
  pairs <- make_overlapping_pairs(2L, 1L, c("r1", "r2"))
  write_qual_fastq(pairs$fwd, raw)
  write_qual_fastq(pairs$fwd, trim)
  uc <- vsearch_cluster_unoise2(empty, min_size = 8, threads = 1, shards = 1)[[
    1
  ]]
  seq_all <- Biostrings::DNAStringSet(pairs$merged)
  names(seq_all) <- "1"
  dm <- make_denoise_map(uc, seq_all)
  smap <- unoise_read_map(
    sample = "empty",
    fq_raw = raw,
    fq_trim = trim,
    fq_merged = empty,
    uc = uc,
    denoise_map = dm
  )
  expect_equal(nrow(smap), 2L)
  expect_true(all(is.na(smap$seq_idx)))
  expect_true(all(bitwAnd(as.integer(smap$flags), 0x04) == 0L))
})
