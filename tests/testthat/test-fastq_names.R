test_that("fastq_names works", {
  expect_equal(fastq_names(raw1_fastq), test_names1)
  expect_equal(fastq_names(raw1_fastq_gz), test_names1)
})

test_that("fastq_rename renames FASTQ read names and preserves content", {
  new_names <- sprintf("renamed_fastq_%03d", seq_along(test_names1))

  out_fastq <- tempfile(fileext = ".fastq")
  out_fastq_gz <- tempfile(fileext = ".fastq.gz")

  expect_equal(fastq_rename(raw1_fastq, new_names, out_fastq), out_fastq)
  expect_equal(fastq_rename(raw1_fastq_gz, new_names, out_fastq_gz), out_fastq_gz)

  renamed <- Biostrings::readQualityScaledDNAStringSet(out_fastq)
  renamed_gz <- Biostrings::readQualityScaledDNAStringSet(out_fastq_gz)
  original <- Biostrings::readQualityScaledDNAStringSet(raw1_fastq)

  expect_equal(names(renamed), new_names)
  expect_equal(names(renamed_gz), new_names)
  expect_equal(unname(as.character(renamed)), unname(as.character(original)))
  expect_equal(
    unname(as.character(Biostrings::quality(renamed))),
    unname(as.character(Biostrings::quality(original)))
  )
})

test_that("fasta_rename renames FASTA read names and preserves content", {
  new_names <- sprintf("renamed_fasta_%03d", seq_along(test_names1))

  out_fasta <- tempfile(fileext = ".fasta")
  out_fasta_gz <- tempfile(fileext = ".fasta.gz")

  expect_equal(fasta_rename(raw1_fasta, new_names, out_fasta), out_fasta)
  expect_equal(fasta_rename(raw1_fasta_gz, new_names, out_fasta_gz), out_fasta_gz)

  renamed <- Biostrings::readDNAStringSet(out_fasta)
  renamed_gz <- Biostrings::readDNAStringSet(out_fasta_gz)
  original <- Biostrings::readDNAStringSet(raw1_fasta)

  expect_equal(names(renamed), new_names)
  expect_equal(names(renamed_gz), new_names)
  expect_equal(unname(as.character(renamed)), unname(as.character(original)))
})

test_that("fasta_rename works with multi-line FASTA sequence records", {
  new_names <- sprintf("renamed_wrapped_%03d", seq_along(test_names1))

  wrapped_in <- tempfile(fileext = ".fasta")
  wrapped_in_gz <- tempfile(fileext = ".fasta.gz")
  Biostrings::writeXStringSet(test_dss1, wrapped_in, width = 7L)
  Biostrings::writeXStringSet(test_dss1, wrapped_in_gz, width = 7L, compress = "gzip")

  out_wrapped <- tempfile(fileext = ".fasta")
  out_wrapped_gz <- tempfile(fileext = ".fasta.gz")

  expect_equal(fasta_rename(wrapped_in, new_names, out_wrapped), out_wrapped)
  expect_equal(fasta_rename(wrapped_in_gz, new_names, out_wrapped_gz), out_wrapped_gz)

  renamed <- Biostrings::readDNAStringSet(out_wrapped)
  renamed_gz <- Biostrings::readDNAStringSet(out_wrapped_gz)
  original <- Biostrings::readDNAStringSet(wrapped_in)

  expect_equal(names(renamed), new_names)
  expect_equal(names(renamed_gz), new_names)
  expect_equal(unname(as.character(renamed)), unname(as.character(original)))
  expect_equal(unname(as.character(renamed_gz)), unname(as.character(original)))
})
