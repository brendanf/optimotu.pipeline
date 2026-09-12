test_that("seq_batch_make_chunk_files returns character(0) for empty FASTA", {
  tf <- tempfile(fileext = ".fasta")
  write_sequence(Biostrings::DNAStringSet(), tf)
  on.exit(unlink(tf), add = TRUE)
  e <- environment()
  out <- optimotu.pipeline:::seq_batch_make_chunk_files(
    seqs = tf,
    files = NULL,
    seq_idx = NULL,
    ncpu = 2L,
    local_envir = e
  )
  expect_identical(out, character(0))
})

test_that("seq_batch_make_chunk_files returns character(0) for empty DNAStringSet", {
  e <- environment()
  out <- optimotu.pipeline:::seq_batch_make_chunk_files(
    seqs = Biostrings::DNAStringSet(),
    files = NULL,
    seq_idx = NULL,
    ncpu = 2L,
    local_envir = e
  )
  expect_identical(out, character(0))
})

test_that("seq_batch_make_chunk_files preserves multi-file record order (ncpu=1)", {
  e <- environment()
  f1 <- tempfile(fileext = ".fasta")
  f2 <- tempfile(fileext = ".fasta")
  on.exit(unlink(c(f1, f2)), add = TRUE)
  s1 <- Biostrings::DNAStringSet("AAAA")
  names(s1) <- "AAAA"
  s2 <- Biostrings::DNAStringSet("TTTT")
  names(s2) <- "TTTT"
  write_sequence(s1, f1)
  write_sequence(s2, f2)
  out <- optimotu.pipeline:::seq_batch_make_chunk_files(
    seqs = c(f1, f2),
    files = NULL,
    seq_idx = NULL,
    ncpu = 1L,
    local_envir = e
  )
  expect_length(out, 1L)
  got <- Biostrings::readDNAStringSet(out[[1L]])
  expect_length(got, 2L)
  expect_equal(names(got), c("AAAA", "TTTT"))
})

test_that("vsearch_uchime_ref rejects list seq_idx with more than one partition", {
  # find_vsearch() is evaluated in the default args before the partition check.
  skip_if_no_vsearch()
  f1 <- tempfile(fileext = ".fasta")
  f2 <- tempfile(fileext = ".fasta")
  on.exit(unlink(c(f1, f2)), add = TRUE)
  d1 <- Biostrings::DNAStringSet("ACGT")
  names(d1) <- "a"
  d2 <- Biostrings::DNAStringSet("ACGT")
  names(d2) <- "b"
  write_sequence(d1, f1)
  write_sequence(d2, f2)
  expect_error(
    vsearch_uchime_ref(
      query = f1,
      ref = f2,
      seq_idx = list(1L, 1L)
    ),
    "partition"
  )
})

test_that("sintax rejects list seq_idx with more than one partition", {
  # find_vsearch() is evaluated in the default args before the partition check.
  skip_if_no_vsearch()
  f1 <- tempfile(fileext = ".fasta")
  f2 <- tempfile(fileext = ".fasta")
  on.exit(unlink(c(f1, f2)), add = TRUE)
  d1 <- Biostrings::DNAStringSet("ACGT")
  names(d1) <- "a"
  d2 <- Biostrings::DNAStringSet("ACGT")
  names(d2) <- "b"
  write_sequence(d1, f1)
  write_sequence(d2, f2)
  expect_error(
    sintax(
      query = f1,
      ref = f2,
      ncpu = 1L,
      seq_idx = list(1L, 1L)
    ),
    "partition"
  )
})

test_that("seq_batch_character reads a single FASTA path", {
  fa <- tempfile(fileext = ".fasta")
  on.exit(unlink(fa), add = TRUE)
  write_sequence(Biostrings::DNAStringSet(c(x = "ACGT", y = "TGCA")), fa)
  got <- optimotu.pipeline:::seq_batch_character(fa)
  expect_type(got, "character")
  expect_named(got, c("x", "y"))
  expect_equal(unname(got), c("ACGT", "TGCA"))
})

test_that("seq_batch_character concatenates multiple FASTA paths in order", {
  f1 <- tempfile(fileext = ".fasta")
  f2 <- tempfile(fileext = ".fasta")
  on.exit(unlink(c(f1, f2)), add = TRUE)
  write_sequence(Biostrings::DNAStringSet(c(a1 = "AAAA")), f1)
  write_sequence(Biostrings::DNAStringSet(c(b1 = "TTTT")), f2)
  got <- optimotu.pipeline:::seq_batch_character(c(f1, f2))
  expect_named(got, c("a1", "b1"))
  expect_equal(unname(got), c("AAAA", "TTTT"))
})

test_that("seq_batch_character handles DNAStringSet with seq_idx", {
  d <- Biostrings::DNAStringSet(c(s1 = "AA", s2 = "CC", s3 = "GG"))
  got <- optimotu.pipeline:::seq_batch_character(d, seq_idx = c(3L, 1L, 3L))
  expect_equal(unname(got), c("GG", "AA", "GG"))
  expect_equal(names(got), c("s3", "s1", "s3"))
})

test_that("seq_batch_character handles data.frame with seq_idx", {
  df <- tibble::tibble(
    seq_id = c("u", "v", "w"),
    seq = c("AAA", "CCC", "TTT")
  )
  got <- optimotu.pipeline:::seq_batch_character(df, seq_idx = c(2L, 2L, 1L))
  expect_named(got, c("v", "v", "u"))
  expect_equal(unname(got), c("CCC", "CCC", "AAA"))
})

test_that("seq_batch_character handles named literal character sequences", {
  got <- optimotu.pipeline:::seq_batch_character(c(a = "ACGT", b = "NNNN"))
  expect_equal(got, c(a = "ACGT", b = "NNNN"))
})

test_that("seq_batch_character subsets named literal character with seq_idx", {
  got <- optimotu.pipeline:::seq_batch_character(
    c(a = "ACGT", b = "TGCA", c = "GGGG"),
    seq_idx = c(3L, 1L, 1L)
  )
  expect_equal(got, c(c = "GGGG", a = "ACGT", a = "ACGT"))
})

test_that("seq_batch_character returns character(0) for empty DNAStringSet", {
  got <- optimotu.pipeline:::seq_batch_character(Biostrings::DNAStringSet())
  expect_identical(got, character())
})

test_that("seq_batch_character errors on unsupported input type", {
  expect_error(
    optimotu.pipeline:::seq_batch_character(list(a = "ACGT")),
    "Unsupported sequence input type"
  )
})

test_that("seq_batch_character handles fastqindexr_index without seq_idx", {
  fa <- tempfile(fileext = ".fasta")
  on.exit(unlink(fa), add = TRUE)
  write_sequence(Biostrings::DNAStringSet(c(p = "ACGT", q = "TGCA")), fa)
  idx <- fastqindexr::create_index(files = fa, type = "fasta")
  got <- optimotu.pipeline:::seq_batch_character(idx)
  expect_named(got, c("p", "q"))
  expect_equal(unname(got), c("ACGT", "TGCA"))
})

test_that("seq_batch_character handles fastqindexr_index with seq_idx", {
  fa <- tempfile(fileext = ".fasta")
  on.exit(unlink(fa), add = TRUE)
  write_sequence(
    Biostrings::DNAStringSet(c(p = "ACGT", q = "TGCA", r = "GGGG")),
    fa
  )
  idx <- fastqindexr::create_index(files = fa, type = "fasta")
  got <- optimotu.pipeline:::seq_batch_character(idx, seq_idx = c(3L, 1L, 3L))
  expect_equal(unname(got), c("GGGG", "ACGT", "GGGG"))
  expect_equal(names(got), c("r", "p", "r"))
})

test_that("seq_batch_character uses file override for fastqindexr_index", {
  fa <- tempfile(fileext = ".fasta")
  on.exit(unlink(fa), add = TRUE)
  write_sequence(Biostrings::DNAStringSet(c(z = "NNNN")), fa)
  idx <- fastqindexr::create_index(files = fa, type = "fasta")
  got <- optimotu.pipeline:::seq_batch_character(idx, files = fa)
  expect_equal(got, c(z = "NNNN"))
})
