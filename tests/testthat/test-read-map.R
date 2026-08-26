test_that("match_to_seq_all accepts FASTA path, DNAStringSet, and character", {
  seqs <- c("ACGTAAAATTTT", "GGGGCCCCAAAA")
  names(seqs) <- c("1", "2")
  dss <- Biostrings::DNAStringSet(seqs)
  fasta <- withr::local_tempfile(fileext = ".fasta")
  fasta_gz <- withr::local_tempfile(fileext = ".fasta.gz")
  Biostrings::writeXStringSet(dss, fasta)
  Biostrings::writeXStringSet(dss, fasta_gz, compress = TRUE)

  expect_equal(match_to_seq_all(seqs[[1]], fasta), 1L)
  expect_equal(match_to_seq_all(seqs[[2]], fasta_gz), 2L)
  expect_equal(match_to_seq_all(seqs[[1]], dss), 1L)
  expect_equal(match_to_seq_all(seqs[[2]], unname(seqs)), 2L)
  expect_true(is.na(match_to_seq_all("TTTTTTTTTTTT", dss)))
  expect_true(is.na(match_to_seq_all(NA_character_, dss)))
})

test_that("match_to_seq_all reverse-complements queries when rc is TRUE", {
  seqs <- c("ACGTAAAATTTT", "GGGGCCCCAAAA")
  names(seqs) <- c("1", "2")
  dss <- Biostrings::DNAStringSet(seqs)
  rc1 <- as.character(Biostrings::reverseComplement(dss[1]))
  expect_equal(match_to_seq_all(rc1, dss, rc = TRUE), 1L)
  expect_true(is.na(match_to_seq_all(rc1, dss, rc = FALSE)))
})

test_that("seq_idx_lookup round-trips through match_to_seq_all", {
  seqs <- c("ACGTAAAATTTT", "GGGGCCCCAAAA", "TTTTCCCCGGGG")
  names(seqs) <- as.character(seq_along(seqs))
  dss <- Biostrings::DNAStringSet(seqs)
  queries <- c(seqs[[2]], NA_character_, seqs[[1]], seqs[[2]])
  lookup <- seq_idx_lookup(queries, dss)
  expect_s3_class(lookup, "seq_idx_lookup")
  expect_equal(lookup$keys, c(seqs[[2]], seqs[[1]]))
  expect_equal(match_to_seq_all(queries, lookup), c(2L, NA, 1L, 2L))
  expect_equal(
    match_to_seq_all(queries, lookup),
    match_to_seq_all(queries, dss)
  )
})

test_that("seq_idx_lookup guards rc consistency", {
  seqs <- c("ACGTAAAATTTT", "GGGGCCCCAAAA")
  dss <- Biostrings::DNAStringSet(seqs)
  lookup <- seq_idx_lookup(seqs[[1]], dss, rc = FALSE)
  expect_error(
    match_to_seq_all(seqs[[1]], lookup, rc = TRUE),
    "does not match seq_idx_lookup"
  )
  rc_query <- as.character(Biostrings::reverseComplement(dss[1]))
  lookup_rc <- seq_idx_lookup(rc_query, dss, rc = TRUE)
  expect_equal(match_to_seq_all(rc_query, lookup_rc, rc = TRUE), 1L)
})

test_that("Biostrings reverseComplement agrees with dada2::rc on IUPAC", {
  skip_if_not_installed("dada2")
  seqs <- c("ACGTMRWSYKVHDBN", "GGGGCCCCAAAATT")
  bs <- as.character(Biostrings::reverseComplement(
    Biostrings::DNAStringSet(seqs)
  ))
  d2 <- dada2::rc(seqs)
  expect_equal(bs, d2)
})

test_that("make_denoise_map.list matches once and projects to seqtable", {
  seqs <- Biostrings::DNAStringSet(c(
    "ACGTAAAATTTT",
    "GGGGCCCCAAAA",
    "TTTTCCCCGGGG"
  ))
  names(seqs) <- as.character(seq_along(seqs))
  uc1 <- structure(
    list(
      clusters = tibble::tibble(
        clust_idx = 0:1,
        size = c(10L, 5L),
        seq = as.character(seqs[1:2])
      ),
      map = tibble::tibble(clust_idx = integer(), seq_id = character())
    ),
    class = "uc_cluster"
  )
  uc2 <- structure(
    list(
      clusters = tibble::tibble(
        clust_idx = 0L,
        size = 7L,
        seq = as.character(seqs[2])
      ),
      map = tibble::tibble(clust_idx = integer(), seq_id = character())
    ),
    class = "uc_cluster"
  )
  dm <- make_denoise_map(list(s1 = uc1, s2 = uc2), seqs)
  expect_equal(
    dm,
    tibble::tibble(
      sample = c("s1", "s1", "s2"),
      denoise_idx = c(0L, 1L, 0L),
      seq_idx = c(1L, 2L, 2L),
      nread = c(10L, 5L, 7L)
    )
  )
  out <- denoise_map_to_seqtable(dm)
  expect_equal(
    out,
    tibble::tibble(
      sample = c("s1", "s1", "s2"),
      seq_idx = c(1L, 2L, 2L),
      nread = c(10L, 5L, 7L)
    )
  )
  expect_equal(
    make_mapped_sequence_table(list(s1 = uc1, s2 = uc2), seqs),
    out
  )
  out_rc <- make_mapped_sequence_table(
    list(s1 = uc1),
    Biostrings::reverseComplement(seqs),
    rc = TRUE
  )
  expect_equal(out_rc$seq_idx, c(1L, 2L))
})

test_that("make_denoise_map.data.frame keeps pre-accept row indices", {
  merged <- tibble::tibble(
    sequence = c("AAA", "CCC", "GGG"),
    abundance = c(10L, 5L, 1L),
    accept = c(TRUE, FALSE, TRUE)
  )
  seqs <- Biostrings::DNAStringSet(c("AAA", "CCC", "GGG"))
  names(seqs) <- as.character(seq_along(seqs))
  dm <- make_denoise_map(merged, seqs)
  expect_equal(dm$denoise_idx, c(1L, 3L))
  expect_equal(dm$seq_idx, c(1L, 3L))
  expect_equal(dm$nread, c(10L, 1L))
})

test_that("make_denoise_map retains unmatched seq_idx rows", {
  uc <- structure(
    list(
      clusters = tibble::tibble(
        clust_idx = 0:1,
        size = c(10L, 5L),
        seq = c("ACGTAAAATTTT", "NNNNNNNNNNNN")
      ),
      map = tibble::tibble(clust_idx = integer(), seq_id = character())
    ),
    class = "uc_cluster"
  )
  seqs <- Biostrings::DNAStringSet("ACGTAAAATTTT")
  names(seqs) <- "1"
  dm <- make_denoise_map(uc, seqs)
  expect_equal(dm$denoise_idx, 0:1)
  expect_equal(dm$seq_idx, c(1L, NA_integer_))
  expect_equal(
    denoise_map_to_seqtable(dm),
    tibble::tibble(seq_idx = 1L, nread = 10L)
  )
})

test_that("add_lulu_to_read_map rewrites daughters and keeps prelulu_idx", {
  read_map <- tibble::tibble(
    sample = "s1",
    raw_idx = 1:4,
    seq_idx = c(1L, 2L, 1L, NA_integer_),
    flags = as.raw(c(0x07, 0x07, 0x07, 0x03))
  )
  lulu_map <- tibble::tibble(
    seq_idx = c(1L, 2L),
    lulu_idx = c(1L, 1L)
  )
  out <- add_lulu_to_read_map(read_map, lulu_map)
  expect_named(out, c("sample", "raw_idx", "seq_idx", "prelulu_idx", "flags"))
  expect_equal(out$seq_idx, c(1L, 1L, 1L, NA_integer_))
  expect_equal(out$prelulu_idx, c(1L, 2L, 1L, NA_integer_))
  expect_equal(out$flags, read_map$flags)
  expect_true(out$prelulu_idx[1] == out$seq_idx[1])
  expect_true(out$prelulu_idx[2] != out$seq_idx[2])
  expect_true(is.na(out$prelulu_idx[4]))
})

test_that("remove_tag_jumps keeps seq_idx and add_uncross joins LULU daughters", {
  seqtable <- tibble::tibble(
    sample = c("s1", "s1", "s2", "s2"),
    seq_idx = c(1L, 3L, 1L, 3L),
    nread = c(1000L, 1L, 1L, 1000L)
  )
  uncross <- remove_tag_jumps(seqtable, f = 0.01, p = 1, id_col = "seq_idx")
  expect_true("seq_idx" %in% names(uncross))
  expect_equal(uncross$seq_idx, seqtable$seq_idx)
  expect_equal(uncross$sample, seqtable$sample)
  expect_equal(uncross$is_tag_jump, c(FALSE, TRUE, TRUE, FALSE))

  read_map <- tibble::tibble(
    sample = "s1",
    raw_idx = 1:4,
    seq_idx = c(1L, 2L, NA_integer_, 3L),
    flags = as.raw(c(0x07, 0x07, 0x03, 0x07))
  )
  lulu_map <- tibble::tibble(
    seq_idx = c(1L, 2L, 3L),
    lulu_idx = c(1L, 1L, 3L)
  )
  annotated <- read_map |>
    add_lulu_to_read_map(lulu_map) |>
    add_uncross_to_read_map(seqtable, uncross)

  expect_true("prelulu_idx" %in% names(annotated))
  # parent of 1 survived uncross
  expect_equal(bitwAnd(as.integer(annotated$flags[1]), 0x08), 0x08)
  # daughter 2 joins via parent 1 and also survives
  expect_equal(annotated$seq_idx[2], 1L)
  expect_equal(annotated$prelulu_idx[2], 2L)
  expect_equal(bitwAnd(as.integer(annotated$flags[2]), 0x08), 0x08)
  # never denoised: no 0x04, no 0x08
  expect_equal(bitwAnd(as.integer(annotated$flags[3]), 0x04), 0L)
  expect_equal(bitwAnd(as.integer(annotated$flags[3]), 0x08), 0L)
  # parent 3 is a tag-jump in s1
  expect_equal(bitwAnd(as.integer(annotated$flags[4]), 0x04), 0x04)
  expect_equal(bitwAnd(as.integer(annotated$flags[4]), 0x08), 0L)
})

test_that("add_uncross_to_read_map falls back to positional keys", {
  seqtable <- tibble::tibble(
    sample = c("s1", "s2"),
    seq_idx = c(1L, 1L),
    nread = c(1000L, 1L)
  )
  uncross <- remove_tag_jumps(seqtable, f = 0.01, p = 1, id_col = "seq_idx")
  uncross_noid <- uncross[c(
    "sample",
    "nread",
    "total",
    "uncross",
    "is_tag_jump"
  )]
  read_map <- tibble::tibble(
    sample = c("s1", "s2"),
    raw_idx = 1L,
    seq_idx = 1L,
    flags = as.raw(0x07),
    prelulu_idx = 1L
  )
  out <- add_uncross_to_read_map(read_map, seqtable, uncross_noid)
  expect_equal(out$prelulu_idx, c(1L, 1L))
  expect_equal(bitwAnd(as.integer(out$flags[1]), 0x08), 0x08)
  expect_equal(bitwAnd(as.integer(out$flags[2]), 0x08), 0L)
})

test_that("with_read_map_annotate is a no-op when LULU and UNCROSS are off", {
  withr::local_options(
    optimotu.pipeline.do_lulu = FALSE,
    optimotu.pipeline.do_tag_jump = FALSE
  )
  expr <- quote(dada2_read_map(x))
  expect_identical(with_read_map_annotate(expr), expr)
})

test_that("with_read_map_annotate pipes LULU then UNCROSS without evaluating", {
  withr::local_options(
    optimotu.pipeline.do_lulu = TRUE,
    optimotu.pipeline.do_tag_jump = TRUE
  )
  got <- with_read_map_annotate(quote(dada2_read_map(x)))
  expect_true(is.call(got))
  expect_equal(
    got,
    quote(
      dada2_read_map(x) |>
        optimotu.pipeline::add_lulu_to_read_map(lulu_asv_map) |>
        optimotu.pipeline::add_uncross_to_read_map(seqtable_lulu, uncross)
    )
  )
})

test_that("with_read_map_annotate pipes LULU only", {
  withr::local_options(
    optimotu.pipeline.do_lulu = TRUE,
    optimotu.pipeline.do_tag_jump = FALSE
  )
  expect_equal(
    with_read_map_annotate(quote(dada2_read_map(x))),
    quote(
      dada2_read_map(x) |>
        optimotu.pipeline::add_lulu_to_read_map(lulu_asv_map)
    )
  )
})

test_that("with_read_map_annotate defaults UNCROSS seqtable without LULU", {
  withr::local_options(
    optimotu.pipeline.do_lulu = FALSE,
    optimotu.pipeline.do_tag_jump = TRUE
  )
  expect_equal(
    with_read_map_annotate(quote(dada2_read_map(x))),
    quote(
      dada2_read_map(x) |>
        optimotu.pipeline::add_uncross_to_read_map(seqtable_raw, uncross)
    )
  )
})

test_that("with_read_map_annotate accepts an explicit seqtable expression", {
  withr::local_options(
    optimotu.pipeline.do_lulu = TRUE,
    optimotu.pipeline.do_tag_jump = TRUE
  )
  expect_equal(
    with_read_map_annotate(
      quote(dada2_read_map(x)),
      seqtable = quote(seqtable_pre_uncross)
    ),
    quote(
      dada2_read_map(x) |>
        optimotu.pipeline::add_lulu_to_read_map(lulu_asv_map) |>
        optimotu.pipeline::add_uncross_to_read_map(
          seqtable_pre_uncross,
          uncross
        )
    )
  )
})
