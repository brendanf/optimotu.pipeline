test_that("fasta_is_hamming_ready accepts equal-length MSA", {
  fa <- withr::local_tempfile(fileext = ".fasta")
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c(s1 = "ACGT-", s2 = "ACG-T")),
    fa
  )
  expect_true(fasta_is_hamming_ready(fa))
})

test_that("fasta_is_hamming_ready accepts A2M with insert columns", {
  fa <- withr::local_tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">s1",
      "ACGT",
      ">s2",
      "ACggGT"
    ),
    fa
  )
  expect_true(fasta_is_hamming_ready(fa))
})

test_that("fasta_is_hamming_ready rejects variable-length unaligned FASTA", {
  fa <- withr::local_tempfile(fileext = ".fasta")
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c(s1 = "ACGT", s2 = "ACGTA")),
    fa
  )
  expect_false(fasta_is_hamming_ready(fa))
})

test_that("fasta_is_hamming_ready rejects empty files", {
  fa <- withr::local_tempfile(fileext = ".fasta")
  writeLines(character(), fa)
  expect_false(fasta_is_hamming_ready(fa))
})

test_that("fasta_is_hamming_ready accepts wrapped multiline FASTA", {
  fa <- withr::local_tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">s1",
      "ACGTACGT",
      "ACGTACGT",
      ">s2",
      "ACGTACGT",
      "ACGTACGT"
    ),
    fa
  )
  expect_true(fasta_is_hamming_ready(fa))

  fa_a2m <- withr::local_tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">s1",
      "ACGT",
      "ACGT",
      ">s2",
      "ACggGT",
      "ACGT"
    ),
    fa_a2m
  )
  expect_true(fasta_is_hamming_ready(fa_a2m))

  fa_wrap <- withr::local_tempfile(fileext = ".fasta")
  long <- paste(rep("ACGT", 20L), collapse = "")
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c(s1 = long, s2 = long)),
    fa_wrap,
    width = 20L
  )
  expect_true(fasta_is_hamming_ready(fa_wrap))
})

test_that("fasta_is_hamming_ready rejects wrapped sequences of different width", {
  fa <- withr::local_tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">s1",
      "ACGT",
      "ACGT",
      ">s2",
      "ACGT",
      "ACGTA"
    ),
    fa
  )
  expect_false(fasta_is_hamming_ready(fa))
})

test_that("fasta_is_hamming_ready streams gzip", {
  fa <- withr::local_tempfile(fileext = ".fasta.gz")
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c(s1 = "ACGT", s2 = "ACGA")),
    fa,
    compress = TRUE
  )
  expect_true(fasta_is_hamming_ready(fa))
})
