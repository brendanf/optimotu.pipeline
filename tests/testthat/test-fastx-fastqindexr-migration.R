make_oneline_fasta_gz <- function(seq, ids) {
  path <- tempfile(fileext = ".fasta.gz")
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(as.vector(rbind(paste0(">", ids), seq)), con = con)
  path
}

local_dna_hmm <- function(seed = "ACGTACGT") {
  fa <- tempfile(fileext = ".fasta")
  hmm <- tempfile(fileext = ".hmm")
  writeLines(c(">seed", seed), fa)
  processx::run(
    "hmmbuild",
    c("--dna", hmm, fa),
    error_on_status = TRUE
  )
  unlink(fa)
  withr::defer(unlink(hmm), envir = parent.frame())
  hmm
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
  expect_equal(names(got), c("0", "1", "0"))
})

test_that("fastx_gz_random_access_extract in-memory renumber is zero_based", {
  seq <- c("AA", "TT", "GG")
  ids <- c("a", "b", "c")
  infile <- make_oneline_fasta_gz(seq, ids)
  idx <- fastqindexr::create_index(files = infile, type = "fasta")
  got <- fastx_gz_random_access_extract(
    infile = infile,
    index = idx,
    i = c(3L, 1L),
    outfile = NULL,
    renumber = TRUE
  )
  expect_s4_class(got, "DNAStringSet")
  expect_equal(unname(as.character(got)), c(seq[3], seq[1]))
  expect_equal(names(got), c("0", "1"))
})

reference_fastx_gz_hash <- function(infile, index, start, n) {
  tmp <- withr::local_tempfile()
  fastqindexr::extract_sequences_to_file(
    index = index,
    seq_idx = seq(start, start + n - 1L),
    file = infile,
    outfile = tmp,
    type = "auto",
    append = FALSE,
    compress = FALSE,
    collapse_sequence_lines = FALSE,
    renumber = "none"
  )
  c(strtrim(system2("md5sum", tmp, stdout = TRUE), 32))
}

test_that("fastx_gz_hash matches md5 of extracted byte stream", {
  skip_if_not(nzchar(Sys.which("md5sum")), "md5sum not on PATH")
  seq <- c("ACGT", "TGCA", "GGGG", "CCCC", "TTAA")
  ids <- c("h1", "h2", "h3", "h4", "h5")
  infile <- make_oneline_fasta_gz(seq, ids)
  idx <- fastqindexr::create_index(files = infile, type = "fasta")
  got <- fastx_gz_hash(infile = infile, index = idx, start = 2L, n = 3L)
  expect_type(got, "character")
  expect_length(got, 1L)
  expect_match(got, "^[0-9a-f]{32}$")
  expect_equal(got, reference_fastx_gz_hash(infile, idx, 2L, 3L))
})

test_that("fastx_gz_hash hashes from the first record and single records", {
  skip_if_not(nzchar(Sys.which("md5sum")), "md5sum not on PATH")
  seq <- c("AAAA", "AAAC", "AATT")
  ids <- c("x1", "x2", "x3")
  infile <- make_oneline_fasta_gz(seq, ids)
  idx <- fastqindexr::create_index(files = infile, type = "fasta")
  expect_equal(
    fastx_gz_hash(infile = infile, index = idx, start = 1L, n = length(seq)),
    reference_fastx_gz_hash(infile, idx, 1L, length(seq))
  )
  expect_equal(
    fastx_gz_hash(infile = infile, index = idx, start = 3L, n = 1L),
    reference_fastx_gz_hash(infile, idx, 3L, 1L)
  )
})

test_that("fastx_gz_hash accepts fqi path indexes", {
  skip_if_not(nzchar(Sys.which("md5sum")), "md5sum not on PATH")
  skip_if(!nzchar(Sys.which("fastqindex")), "fastqindex CLI not available")
  seq <- c("ACGT", "TGCA", "GGGG")
  ids <- c("f1", "f2", "f3")
  infile <- make_oneline_fasta_gz(seq, ids)
  idx <- fastx_gz_index(infile)
  got <- fastx_gz_hash(infile = infile, index = idx, start = 1L, n = 2L)
  expect_equal(got, reference_fastx_gz_hash(infile, idx, 1L, 2L))
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

test_that("hmmalign indexed extraction matches legacy extract-then-align", {
  skip_if_not(nzchar(Sys.which("hmmalign")), "hmmalign not on PATH")
  skip_if_not(nzchar(Sys.which("hmmbuild")), "hmmbuild not on PATH")
  seq <- c("ACGTACGT", "AAAATTTT", "GGGGCCCC")
  ids <- as.character(seq_along(seq))
  infile <- make_oneline_fasta_gz(seq, ids)
  idx <- fastqindexr::create_index(files = infile, type = "fasta")
  req <- c(3L, 1L, 3L)
  sub <- withr::local_tempfile(fileext = ".fasta")
  fastqindexr::extract_sequences_to_file(
    index = idx,
    seq_idx = req,
    file = infile,
    outfile = sub,
    type = "fasta",
    compress = FALSE
  )
  hmm <- local_dna_hmm()
  out_legacy <- withr::local_tempfile(fileext = ".a2m")
  out_idx <- withr::local_tempfile(fileext = ".a2m")
  optimotu.pipeline::hmmalign(
    sub,
    hmm,
    out_legacy,
    outformat = "A2M",
    compress = FALSE
  )
  optimotu.pipeline::hmmalign(
    idx,
    hmm,
    out_idx,
    outformat = "A2M",
    compress = FALSE,
    files = infile,
    seq_idx = req
  )
  expect_equal(
    digest::digest(file = out_legacy, algo = "md5"),
    digest::digest(file = out_idx, algo = "md5")
  )
})

test_that("hmmalign with index uses full file when seq_idx is NULL", {
  skip_if_not(nzchar(Sys.which("hmmalign")), "hmmalign not on PATH")
  skip_if_not(nzchar(Sys.which("hmmbuild")), "hmmbuild not on PATH")
  seq <- c("ACGTACGT", "TTTTAAAA")
  ids <- as.character(seq_along(seq))
  infile <- make_oneline_fasta_gz(seq, ids)
  idx <- fastqindexr::create_index(files = infile, type = "fasta")
  sub <- withr::local_tempfile(fileext = ".fasta")
  fastqindexr::extract_sequences_to_file(
    index = idx,
    seq_idx = seq_along(seq),
    file = infile,
    outfile = sub,
    type = "fasta",
    compress = FALSE
  )
  hmm <- local_dna_hmm()
  out_full <- withr::local_tempfile(fileext = ".a2m")
  out_null <- withr::local_tempfile(fileext = ".a2m")
  optimotu.pipeline::hmmalign(
    sub,
    hmm,
    out_full,
    outformat = "A2M",
    compress = FALSE
  )
  optimotu.pipeline::hmmalign(
    idx,
    hmm,
    out_null,
    outformat = "A2M",
    compress = FALSE,
    files = infile,
    seq_idx = NULL
  )
  expect_equal(
    digest::digest(file = out_full, algo = "md5"),
    digest::digest(file = out_null, algo = "md5")
  )
})

test_that("hmmalign indexed path matches legacy for multi-file logical concat", {
  skip_if_not(nzchar(Sys.which("hmmalign")), "hmmalign not on PATH")
  skip_if_not(nzchar(Sys.which("hmmbuild")), "hmmbuild not on PATH")
  seq1 <- c("ACGTACGT", "AAAATTTT")
  seq2 <- c("GGGGGGGG", "CTCTCTCT")
  f1 <- make_oneline_fasta_gz(seq1, c("a1", "a2"))
  f2 <- make_oneline_fasta_gz(seq2, c("b1", "b2"))
  idx <- fastqindexr::create_index(files = c(f1, f2), type = "fasta")
  req <- c(4L, 1L, 4L)
  sub <- withr::local_tempfile(fileext = ".fasta")
  fastqindexr::extract_sequences_to_file(
    index = idx,
    seq_idx = req,
    file = c(f1, f2),
    outfile = sub,
    type = "fasta",
    compress = FALSE
  )
  hmm <- local_dna_hmm()
  out_legacy <- withr::local_tempfile(fileext = ".a2m")
  out_idx <- withr::local_tempfile(fileext = ".a2m")
  optimotu.pipeline::hmmalign(
    sub,
    hmm,
    out_legacy,
    outformat = "A2M",
    compress = FALSE
  )
  optimotu.pipeline::hmmalign(
    idx,
    hmm,
    out_idx,
    outformat = "A2M",
    compress = FALSE,
    files = c(f1, f2),
    seq_idx = req
  )
  expect_equal(
    digest::digest(file = out_legacy, algo = "md5"),
    digest::digest(file = out_idx, algo = "md5")
  )
})

test_that("hmmalign ncpu chunking matches single-process alignment", {
  skip_if_not(nzchar(Sys.which("hmmalign")), "hmmalign not on PATH")
  skip_if_not(nzchar(Sys.which("hmmbuild")), "hmmbuild not on PATH")
  seq <- c(
    "ACGTACGT",
    "AAAATTTT",
    "GGGGCCCC",
    "CTCTCTCT",
    "TTTTAAAA",
    "GCGCGCGC"
  )
  ids <- as.character(seq_along(seq))
  infile <- make_oneline_fasta_gz(seq, ids)
  idx <- fastqindexr::create_index(files = infile, type = "fasta")
  hmm <- local_dna_hmm()
  out_one <- withr::local_tempfile(fileext = ".a2m")
  out_many <- withr::local_tempfile(fileext = ".a2m")
  optimotu.pipeline::hmmalign(
    idx,
    hmm,
    out_one,
    outformat = "A2M",
    compress = FALSE,
    files = infile,
    seq_idx = NULL,
    ncpu = 1L
  )
  optimotu.pipeline::hmmalign(
    idx,
    hmm,
    out_many,
    outformat = "A2M",
    compress = FALSE,
    files = infile,
    seq_idx = NULL,
    ncpu = 4L
  )
  aln_one <- Biostrings::readBStringSet(out_one)
  aln_many <- Biostrings::readBStringSet(out_many)
  expect_setequal(names(aln_many), names(aln_one))
  aln_one <- aln_one[order(names(aln_one))]
  aln_many <- aln_many[order(names(aln_many))]
  expect_equal(names(aln_many), names(aln_one))
  expect_equal(as.character(aln_many), as.character(aln_one))
})

test_that("hmmalign forwards dots without error (targets tracking)", {
  skip_if_not(nzchar(Sys.which("hmmalign")), "hmmalign not on PATH")
  skip_if_not(nzchar(Sys.which("hmmbuild")), "hmmbuild not on PATH")
  seq <- c("ACGTACGT")
  infile <- make_oneline_fasta_gz(seq, "1")
  idx <- fastqindexr::create_index(files = infile, type = "fasta")
  hmm <- local_dna_hmm()
  out <- withr::local_tempfile(fileext = ".a2m")
  expect_no_error(
    optimotu.pipeline::hmmalign(
      idx,
      hmm,
      out,
      outformat = "A2M",
      compress = FALSE,
      files = infile,
      hash = "not_used"
    )
  )
})

test_that("partition_vector_equal_ncpu splits as expected", {
  expect_equal(
    optimotu.pipeline:::partition_vector_equal_ncpu(1:10, 3L),
    list(1:4, 5:7, 8:10)
  )
  expect_equal(
    optimotu.pipeline:::partition_vector_equal_ncpu(1:3, 8L),
    list(1L, 2L, 3L)
  )
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
