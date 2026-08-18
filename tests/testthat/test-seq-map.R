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

test_that("add_lulu_to_seq_map rewrites daughters and keeps denoise_idx", {
  seqmap <- tibble::tibble(
    sample = "s1",
    raw_idx = 1:4,
    seq_idx = c(1L, 2L, 1L, NA_integer_),
    flags = as.raw(c(0x07, 0x07, 0x07, 0x03))
  )
  lulu_map <- tibble::tibble(
    seq_idx = c(1L, 2L),
    lulu_idx = c(1L, 1L)
  )
  out <- add_lulu_to_seq_map(seqmap, lulu_map)
  expect_named(out, c("sample", "raw_idx", "seq_idx", "denoise_idx", "flags"))
  expect_equal(out$seq_idx, c(1L, 1L, 1L, NA_integer_))
  expect_equal(out$denoise_idx, c(1L, 2L, 1L, NA_integer_))
  expect_equal(out$flags, seqmap$flags)
  expect_true(out$denoise_idx[1] == out$seq_idx[1])
  expect_true(out$denoise_idx[2] != out$seq_idx[2])
  expect_true(is.na(out$denoise_idx[4]))
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

  seqmap <- tibble::tibble(
    sample = "s1",
    raw_idx = 1:4,
    seq_idx = c(1L, 2L, NA_integer_, 3L),
    flags = as.raw(c(0x07, 0x07, 0x03, 0x07))
  )
  lulu_map <- tibble::tibble(
    seq_idx = c(1L, 2L, 3L),
    lulu_idx = c(1L, 1L, 3L)
  )
  annotated <- seqmap |>
    add_lulu_to_seq_map(lulu_map) |>
    add_uncross_to_seq_map(seqtable, uncross)

  expect_true("denoise_idx" %in% names(annotated))
  # parent of 1 survived uncross
  expect_equal(bitwAnd(as.integer(annotated$flags[1]), 0x08), 0x08)
  # daughter 2 joins via parent 1 and also survives
  expect_equal(annotated$seq_idx[2], 1L)
  expect_equal(annotated$denoise_idx[2], 2L)
  expect_equal(bitwAnd(as.integer(annotated$flags[2]), 0x08), 0x08)
  # never denoised: no 0x04, no 0x08
  expect_equal(bitwAnd(as.integer(annotated$flags[3]), 0x04), 0L)
  expect_equal(bitwAnd(as.integer(annotated$flags[3]), 0x08), 0L)
  # parent 3 is a tag-jump in s1
  expect_equal(bitwAnd(as.integer(annotated$flags[4]), 0x04), 0x04)
  expect_equal(bitwAnd(as.integer(annotated$flags[4]), 0x08), 0L)
})

test_that("add_uncross_to_seq_map falls back to positional keys", {
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
  seqmap <- tibble::tibble(
    sample = c("s1", "s2"),
    raw_idx = 1L,
    seq_idx = 1L,
    flags = as.raw(0x07),
    denoise_idx = 1L
  )
  out <- add_uncross_to_seq_map(seqmap, seqtable, uncross_noid)
  expect_equal(out$denoise_idx, c(1L, 1L))
  expect_equal(bitwAnd(as.integer(out$flags[1]), 0x08), 0x08)
  expect_equal(bitwAnd(as.integer(out$flags[2]), 0x08), 0L)
})
