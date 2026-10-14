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

local_dna_cm <- function(seed = "ACGTACGTACGTACGT") {
  fa <- tempfile(fileext = ".fasta")
  cm <- tempfile(fileext = ".cm")
  writeLines(c(">seed", seed), fa)
  # cmbuild needs a Stockholm alignment; use --noss for DNA without SS.
  stk <- tempfile(fileext = ".sto")
  writeLines(
    c(
      "# STOCKHOLM 1.0",
      paste0("seed\t", seed),
      paste0("#=GC RF\t", strrep("x", nchar(seed))),
      "//"
    ),
    stk
  )
  processx::run(
    "cmbuild",
    c("--noss", "-F", cm, stk),
    error_on_status = TRUE
  )
  unlink(c(fa, stk))
  withr::defer(unlink(cm), envir = parent.frame())
  cm
}

make_oneline_fasta <- function(seq, ids) {
  path <- tempfile(fileext = ".fasta")
  writeLines(as.vector(rbind(paste0(">", ids), seq)), path)
  path
}

test_that("align_fasta_batch_to_amplicon_model HMM matches hmmalign", {
  skip_if_not(nzchar(Sys.which("hmmalign")), "hmmalign not on PATH")
  skip_if_not(nzchar(Sys.which("hmmbuild")), "hmmbuild not on PATH")
  seq <- c("ACGTACGT", "AAAATTTT", "GGGGCCCC")
  ids <- as.character(seq_along(seq))
  infile <- make_oneline_fasta(seq, ids)
  hmm <- local_dna_hmm()
  out_direct <- withr::local_tempfile(fileext = ".a2m")
  out_helper <- withr::local_tempfile(fileext = ".a2m")
  optimotu.pipeline::hmmalign(
    infile,
    hmm,
    out_direct,
    outformat = "A2M",
    compress = FALSE,
    seq_idx = c(1L, 3L)
  )
  ret <- optimotu.pipeline::align_fasta_batch_to_amplicon_model(
    seqs = infile,
    model = hmm,
    outfile = out_helper,
    model_type = "HMM",
    seq_idx = c(1L, 3L)
  )
  expect_identical(ret, out_helper)
  expect_equal(
    digest::digest(file = out_direct, algo = "md5"),
    digest::digest(file = out_helper, algo = "md5")
  )
})

test_that("align_fasta_batch_to_amplicon_model rejects HMM sfile", {
  skip_if_not(nzchar(Sys.which("hmmalign")), "hmmalign not on PATH")
  skip_if_not(nzchar(Sys.which("hmmbuild")), "hmmbuild not on PATH")
  hmm <- local_dna_hmm()
  expect_error(
    optimotu.pipeline::align_fasta_batch_to_amplicon_model(
      seqs = Biostrings::DNAStringSet(),
      model = hmm,
      outfile = tempfile(fileext = ".fasta"),
      model_type = "HMM",
      sfile = tempfile(fileext = ".dat")
    ),
    "sfile"
  )
})

test_that("align_fasta_batch_to_amplicon_model empty HMM writes empty FASTA", {
  skip_if_not(nzchar(Sys.which("hmmalign")), "hmmalign not on PATH")
  skip_if_not(nzchar(Sys.which("hmmbuild")), "hmmbuild not on PATH")
  hmm <- local_dna_hmm()
  out <- withr::local_tempfile(fileext = ".fasta")
  ret <- optimotu.pipeline::align_fasta_batch_to_amplicon_model(
    seqs = Biostrings::DNAStringSet(),
    model = hmm,
    outfile = out,
    model_type = "HMM"
  )
  expect_identical(ret, out)
  expect_equal(length(Biostrings::readDNAStringSet(out)), 0L)
})

test_that("align_fasta_batch_to_amplicon_model empty CM writes empty FASTA", {
  skip_if_not(nzchar(Sys.which("cmalign")), "cmalign not on PATH")
  skip_if_not(nzchar(Sys.which("cmbuild")), "cmbuild not on PATH")
  cm <- local_dna_cm()
  out <- withr::local_tempfile(fileext = ".fasta")
  sfile <- withr::local_tempfile(fileext = ".dat")
  ret <- optimotu.pipeline::align_fasta_batch_to_amplicon_model(
    seqs = Biostrings::DNAStringSet(),
    model = cm,
    outfile = out,
    model_type = "CM",
    sfile = sfile
  )
  expect_identical(ret, c(out, sfile))
  expect_equal(length(Biostrings::readDNAStringSet(out)), 0L)
  expect_true(file.exists(sfile))
})

test_that("align_fasta_batch_to_amplicon_model CM writes FASTA and optional sfile", {
  skip_if_not(nzchar(Sys.which("cmalign")), "cmalign not on PATH")
  skip_if_not(nzchar(Sys.which("cmbuild")), "cmbuild not on PATH")
  seed <- "ACGTACGTACGTACGT"
  cm <- local_dna_cm(seed)
  seq <- c(seed, seed)
  ids <- c("s1", "s2")
  infile <- make_oneline_fasta(seq, ids)
  out <- withr::local_tempfile(fileext = ".fasta")
  sfile <- withr::local_tempfile(fileext = ".dat")
  ret <- optimotu.pipeline::align_fasta_batch_to_amplicon_model(
    seqs = infile,
    model = cm,
    outfile = out,
    model_type = "CM",
    sfile = sfile,
    ncpu = 1L
  )
  expect_identical(ret, c(out, sfile))
  aln <- Biostrings::readDNAStringSet(out)
  expect_equal(length(aln), 2L)
  expect_equal(names(aln), ids)
  expect_true(file.size(sfile) > 0)
  scores <- optimotu.pipeline::read_sfile(sfile)
  expect_equal(nrow(scores), 2L)
})

test_that("align_fasta_batch_to_amplicon_model CM without sfile returns outfile only", {
  skip_if_not(nzchar(Sys.which("cmalign")), "cmalign not on PATH")
  skip_if_not(nzchar(Sys.which("cmbuild")), "cmbuild not on PATH")
  seed <- "ACGTACGTACGTACGT"
  cm <- local_dna_cm(seed)
  out <- withr::local_tempfile(fileext = ".fasta.gz")
  ret <- optimotu.pipeline::align_fasta_batch_to_amplicon_model(
    seqs = Biostrings::DNAStringSet(c(s1 = seed)),
    model = cm,
    outfile = out,
    model_type = "CM",
    ncpu = 1L
  )
  expect_identical(ret, out)
  expect_equal(length(Biostrings::readDNAStringSet(out)), 1L)
})

test_that("align_fasta_batch_to_amplicon_model forwards dots for targets", {
  skip_if_not(nzchar(Sys.which("hmmalign")), "hmmalign not on PATH")
  skip_if_not(nzchar(Sys.which("hmmbuild")), "hmmbuild not on PATH")
  hmm <- local_dna_hmm()
  out <- withr::local_tempfile(fileext = ".a2m")
  expect_no_error(
    optimotu.pipeline::align_fasta_batch_to_amplicon_model(
      seqs = Biostrings::DNAStringSet(c(s1 = "ACGTACGT")),
      model = hmm,
      outfile = out,
      model_type = "HMM",
      hash = "ignored-for-targets"
    )
  )
})
