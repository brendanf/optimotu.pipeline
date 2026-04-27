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
