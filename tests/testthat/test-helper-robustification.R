test_that("make_seq_names and name_seqs handle empty inputs", {
  expect_identical(optimotu.pipeline::make_seq_names(0L, "ASV"), character())

  empty_fasta <- withr::local_tempfile(fileext = ".fasta")
  file.create(empty_fasta)
  expect_true(file.exists(empty_fasta))
  expect_identical(
    optimotu.pipeline::name_seqs(empty_fasta, "ASV"),
    empty_fasta
  )
  expect_identical(length(Biostrings::fasta.seqlengths(empty_fasta)), 0L)

  empty_fasta_gz <- withr::local_tempfile(fileext = ".fasta.gz")
  con <- gzfile(empty_fasta_gz, open = "wt")
  close(con)
  expect_identical(
    optimotu.pipeline::name_seqs(empty_fasta_gz, "ASV"),
    empty_fasta_gz
  )
  expect_identical(length(Biostrings::fasta.seqlengths(empty_fasta_gz)), 0L)
})

test_that("nomismatch_hits_vsearch short-circuits on empty sequence file", {
  empty_fasta <- withr::local_tempfile(fileext = ".fasta.gz")
  con <- gzfile(empty_fasta, open = "wt")
  close(con)

  seqtab <- tibble::tibble(seq_idx = integer(), nread = integer())
  out <- optimotu.pipeline:::nomismatch_hits_vsearch(
    seqtab = seqtab,
    seqs = empty_fasta
  )

  expect_identical(names(out), c("query", "hit"))
  expect_true(is.integer(out$query))
  expect_true(is.integer(out$hit))
  expect_identical(nrow(out), 0L)
})

test_that("read_long_sequence_table returns typed empty tibble for zero-byte input", {
  empty_tab <- withr::local_tempfile(fileext = ".tsv")
  file.create(empty_tab)

  out <- optimotu.pipeline:::read_long_sequence_table(empty_tab)

  expect_identical(names(out), c("sample", "seq_id", "nread", "seqrun"))
  expect_true(is.character(out$sample))
  expect_true(is.character(out$seq_id))
  expect_true(is.integer(out$nread))
  expect_true(is.character(out$seqrun))
  expect_identical(nrow(out), 0L)
})

test_that("hmmalign applies seq_idx to plain FASTA paths (no index object)", {
  skip_if_not(nzchar(Sys.which("hmmalign")), "hmmalign not on PATH")
  skip_if_not(nzchar(Sys.which("hmmbuild")), "hmmbuild not on PATH")
  fa <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(">a", "ACGTACGT", ">b", "AAAATTTT", ">c", "GGGGCCCC"), fa)
  sub <- withr::local_tempfile(fileext = ".fasta")
  writeLines(
    c(">c", "GGGGCCCC", ">a", "ACGTACGT", ">c", "GGGGCCCC"),
    sub
  )
  hmm_fa <- withr::local_tempfile(fileext = ".fasta")
  hmm <- withr::local_tempfile(fileext = ".hmm")
  writeLines(c(">seed", "ACGTACGT"), hmm_fa)
  processx::run(
    "hmmbuild",
    c("--dna", hmm, hmm_fa),
    error_on_status = TRUE
  )
  out_sub <- withr::local_tempfile(fileext = ".a2m")
  out_idx <- withr::local_tempfile(fileext = ".a2m")
  optimotu.pipeline::hmmalign(
    sub,
    hmm,
    out_sub,
    outformat = "A2M",
    compress = FALSE
  )
  optimotu.pipeline::hmmalign(
    fa,
    hmm,
    out_idx,
    outformat = "A2M",
    compress = FALSE,
    seq_idx = c(3L, 1L, 3L)
  )
  expect_equal(
    digest::digest(file = out_sub, algo = "md5"),
    digest::digest(file = out_idx, algo = "md5")
  )
})

test_that("hmmalign indexed empty extraction matches legacy empty fasta", {
  skip_if_not(nzchar(Sys.which("hmmalign")), "hmmalign not on PATH")
  skip_if_not(nzchar(Sys.which("hmmbuild")), "hmmbuild not on PATH")
  infile <- withr::local_tempfile(fileext = ".fasta.gz")
  con <- gzfile(infile, open = "wt")
  writeLines(as.vector(rbind(">1", "ACGTACGT")), con = con)
  close(con)
  idx <- fastqindexr::create_index(files = infile, type = "fasta")
  empty_fa <- withr::local_tempfile(fileext = ".fasta")
  file.create(empty_fa)
  fa <- withr::local_tempfile(fileext = ".fasta")
  hmm <- withr::local_tempfile(fileext = ".hmm")
  writeLines(c(">seed", "ACGTACGT"), fa)
  processx::run(
    "hmmbuild",
    c("--dna", hmm, fa),
    error_on_status = TRUE
  )
  out_legacy <- withr::local_tempfile(fileext = ".a2m")
  out_idx <- withr::local_tempfile(fileext = ".a2m")
  optimotu.pipeline::hmmalign(
    empty_fa,
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
    seq_idx = integer()
  )
  expect_equal(
    digest::digest(file = out_legacy, algo = "md5"),
    digest::digest(file = out_idx, algo = "md5")
  )
})

test_that("hmmalign empty input writes first outfile only and returns it", {
  skip_if_not(nzchar(Sys.which("hmmalign")), "hmmalign not on PATH")
  skip_if_not(nzchar(Sys.which("hmmbuild")), "hmmbuild not on PATH")
  infile <- withr::local_tempfile(fileext = ".fasta.gz")
  con <- gzfile(infile, open = "wt")
  writeLines(as.vector(rbind(">1", "ACGTACGT")), con = con)
  close(con)
  idx <- fastqindexr::create_index(files = infile, type = "fasta")
  fa <- withr::local_tempfile(fileext = ".fasta")
  hmm <- withr::local_tempfile(fileext = ".hmm")
  writeLines(c(">seed", "ACGTACGT"), fa)
  processx::run(
    "hmmbuild",
    c("--dna", hmm, fa),
    error_on_status = TRUE
  )
  out_a <- withr::local_tempfile(fileext = ".a2m")
  out_b <- withr::local_tempfile(fileext = ".a2m")
  ret <- optimotu.pipeline::hmmalign(
    idx,
    hmm,
    c(out_a, out_b),
    outformat = "A2M",
    compress = FALSE,
    files = infile,
    seq_idx = integer(),
    ncpu = 2L
  )
  expect_identical(ret, out_a)
  expect_true(file.exists(out_a))
  expect_false(file.exists(out_b))
})

test_that("hmmalign multi-outfile length matches ncpu and returns only written paths", {
  skip_if_not(nzchar(Sys.which("hmmalign")), "hmmalign not on PATH")
  skip_if_not(nzchar(Sys.which("hmmbuild")), "hmmbuild not on PATH")
  seq <- c("ACGTACGT", "AAAATTTT", "GGGGCCCC")
  ids <- as.character(seq_along(seq))
  infile <- tempfile(fileext = ".fasta.gz")
  con <- gzfile(infile, open = "wt")
  writeLines(as.vector(rbind(paste0(">", ids), seq)), con = con)
  close(con)
  on.exit(unlink(infile), add = TRUE)
  idx <- fastqindexr::create_index(files = infile, type = "fasta")
  hmm <- tempfile(fileext = ".hmm")
  fa <- tempfile(fileext = ".fasta")
  writeLines(c(">seed", "ACGTACGT"), fa)
  processx::run("hmmbuild", c("--dna", hmm, fa), error_on_status = TRUE)
  on.exit(unlink(c(hmm, fa)), add = TRUE)
  outs <- replicate(4L, tempfile(fileext = ".a2m"))
  on.exit(unlink(outs[file.exists(outs)]), add = TRUE)
  ret <- optimotu.pipeline::hmmalign(
    idx,
    hmm,
    outs,
    outformat = "A2M",
    compress = FALSE,
    files = infile,
    seq_idx = NULL,
    ncpu = 4L
  )
  expect_length(ret, 3L)
  expect_identical(ret, outs[seq_len(3L)])
  expect_true(all(file.exists(outs[seq_len(3L)])))
  expect_false(file.exists(outs[[4L]]))
})

test_that("hmmsearch and nhmmer return typed empties for empty fasta", {
  hmmsearch_exec <- Sys.which("hmmsearch")
  nhmmer_exec <- Sys.which("nhmmer")
  if (!nzchar(hmmsearch_exec) || !nzchar(nhmmer_exec)) {
    skip("Missing hmmsearch or nhmmer executable")
  }

  empty_fasta <- withr::local_tempfile(fileext = ".fasta")
  file.create(empty_fasta)
  hmm <- withr::local_tempfile(fileext = ".hmm")
  writeLines("HMMER3/f [3.1b2 | February 2015]", hmm)

  out_hmm <- optimotu.pipeline::hmmsearch(empty_fasta, hmm)
  expect_identical(nrow(out_hmm), 0L)
  expect_identical(
    names(out_hmm),
    c(
      "seq_name",
      "seq_accno",
      "seq_length",
      "hmm_name",
      "hmm_accno",
      "hmm_length",
      "Evalue",
      "full_score",
      "full_bias",
      "hit_num",
      "total_hits",
      "c_Evalue",
      "i_Evalue",
      "hit_score",
      "hit_bias",
      "hmm_from",
      "hmm_to",
      "seq_from",
      "seq_to",
      "env_from",
      "env_to",
      "acc",
      "description"
    )
  )

  out_nhmmer <- optimotu.pipeline::nhmmer(empty_fasta, hmm, ncpu = 1L)
  expect_identical(nrow(out_nhmmer), 0L)
  expect_identical(
    names(out_nhmmer),
    c(
      "seq_name",
      "seq_accno",
      "hmm_name",
      "hmm_accno",
      "hmm_from",
      "hmm_to",
      "seq_from",
      "seq_to",
      "env_from",
      "env_to",
      "seq_len",
      "strand",
      "Evalue",
      "bit_score",
      "bias",
      "description"
    )
  )
})
