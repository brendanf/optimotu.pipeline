make_oneline_fasta_gz <- function(seq, ids) {
  path <- tempfile(fileext = ".fasta.gz")
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(as.vector(rbind(paste0(">", ids), seq)), con = con)
  path
}

test_that("fastx wrappers accept fastqindexr index objects", {
  seq <- c("ACGT", "TGCA", "GGGG", "CCCC")
  ids <- c("s1", "s2", "s3", "s4")
  infile <- make_oneline_fasta_gz(seq, ids)
  idx <- fastqindexr::create_index(files = infile, type = "fasta")
  req <- c(2L, 1L, 2L)
  got <- fastx_gz_random_access_extract(
    infile = infile,
    index = idx,
    i = req,
    outfile = NULL
  )
  expect_s4_class(got, "DNAStringSet")
  expect_equal(unname(as.character(got)), unname(seq[req]))
  expect_equal(names(got), ids[req])
})

test_that("fastx wrappers accept fqi path indexes", {
  skip_if(!nzchar(Sys.which("fastqindex")), "fastqindex CLI not available")
  seq <- c("ACGT", "TGCA", "GGGG")
  ids <- c("p1", "p2", "p3")
  infile <- make_oneline_fasta_gz(seq, ids)
  idx <- fastx_gz_index(infile)
  req <- c(3L, 1L)
  got <- fastx_gz_random_access_extract(
    infile = infile,
    index = idx,
    i = req,
    outfile = NULL
  )
  expect_equal(as.character(got), seq[req])
  expect_equal(names(got), ids[req])
})

test_that("wrappers handle multi-file inputs with logical concatenation", {
  seq1 <- c("AAAA", "AAAC")
  ids1 <- c("a1", "a2")
  seq2 <- c("TTTT", "TTTG")
  ids2 <- c("b1", "b2")
  infile1 <- make_oneline_fasta_gz(seq1, ids1)
  infile2 <- make_oneline_fasta_gz(seq2, ids2)
  idx <- fastqindexr::create_index(files = c(infile1, infile2), type = "fasta")
  req <- c(1L, 2L, 3L, 4L, 3L)
  got <- fastx_gz_random_access_extract(
    infile = c(infile1, infile2),
    index = idx,
    i = req,
    outfile = NULL
  )
  expect_equal(
    unname(as.character(got)),
    unname(c(seq1, seq2, seq2[1]))
  )
  expect_equal(names(got), c(ids1, ids2, ids2[1]))
})

test_that("fastx_gz_extract preserves duplicate and request order", {
  seq <- c("ACGT", "TGCA", "GGGG", "CCCC", "TTAA")
  ids <- c("o1", "o2", "o3", "o4", "o5")
  infile <- make_oneline_fasta_gz(seq, ids)
  idx <- fastqindexr::create_index(files = infile, type = "fasta")
  req <- c(5L, 2L, 5L, 1L)
  out <- tempfile(fileext = ".fasta")
  expect_equal(
    fastx_gz_extract(infile = infile, index = idx, i = req, outfile = out),
    out
  )
  got <- Biostrings::readDNAStringSet(out)
  expect_equal(unname(as.character(got)), unname(seq[req]))
  expect_equal(names(got), ids[req])
})

test_that("fastx_gz_extract renumber and append follow compatibility behavior", {
  seq <- c("ACGT", "TGCA", "GGGG", "CCCC")
  ids <- c("r1", "r2", "r3", "r4")
  infile <- make_oneline_fasta_gz(seq, ids)
  idx <- fastqindexr::create_index(files = infile, type = "fasta")
  out <- tempfile(fileext = ".fasta")
  fastx_gz_extract(
    infile = infile,
    index = idx,
    i = c(1L, 2L),
    outfile = out,
    renumber = TRUE
  )
  fastx_gz_extract(
    infile = infile,
    index = idx,
    i = c(3L, 4L),
    outfile = out,
    renumber = TRUE,
    append = TRUE
  )
  got <- Biostrings::readDNAStringSet(out)
  expect_equal(names(got), c("0", "1", "0", "1"))
})

test_that("fastx_gz_random_access_extract renumber and append work", {
  seq <- c("ACGT", "TGCA", "GGGG")
  ids <- c("q1", "q2", "q3")
  infile <- make_oneline_fasta_gz(seq, ids)
  idx <- fastqindexr::create_index(files = infile, type = "fasta")
  out <- tempfile(fileext = ".fasta")
  fastx_gz_random_access_extract(
    infile = infile,
    index = idx,
    i = c(1L, 2L),
    outfile = out,
    renumber = TRUE
  )
  fastx_gz_random_access_extract(
    infile = infile,
    index = idx,
    i = 3L,
    outfile = out,
    renumber = TRUE,
    append = TRUE
  )
  got <- Biostrings::readDNAStringSet(out)
  expect_equal(names(got), c("1", "2", "1"))
})

test_that("empty extraction requests keep wrapper edge behavior", {
  seq <- c("ACGT", "TGCA")
  ids <- c("e1", "e2")
  infile <- make_oneline_fasta_gz(seq, ids)
  idx <- fastqindexr::create_index(files = infile, type = "fasta")
  out <- tempfile(fileext = ".fasta")
  expect_equal(
    fastx_gz_extract(
      infile = infile,
      index = idx,
      i = integer(),
      outfile = out
    ),
    out
  )
  expect_true(file.exists(out))
  expect_equal(file.info(out)$size, 0)
  got <- fastx_gz_random_access_extract(
    infile = infile,
    index = idx,
    i = integer(),
    outfile = NULL
  )
  expect_s4_class(got, "DNAStringSet")
  expect_length(got, 0L)
})

test_that("lulu_distmx works with fastqindexr index object input", {
  seq <- c("AAAA", "AAAT", "AATT")
  ids <- c("1", "2", "3")
  seq_file <- make_oneline_fasta_gz(seq, ids)
  seq_idx <- fastqindexr::create_index(files = seq_file, type = "fasta")
  seqtable <- tibble::tibble(seq_idx = 1L:3L, nread = c(10L, 8L, 4L))
  out <- lulu_distmx(
    seqall_file = seq_file,
    seqall_index = seq_idx,
    seqtable = seqtable,
    threshold = 1
  )
  expect_s3_class(out, "data.frame")
  expect_true(all(c("seq_idx1", "seq_idx2", "dist", "nread1") %in% names(out)))
})
