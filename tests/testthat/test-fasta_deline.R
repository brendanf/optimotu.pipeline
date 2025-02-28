test_that("fasta_deline works", {
  expect_no_error(delined_fasta <- fasta_deline(raw2_fasta, withr::local_tempfile(fileext = ".fasta")))
  expect_no_error(delined <- Biostrings::readDNAStringSet(delined_fasta))
  expect_true(all(delined == test_dss2))
  expect_equal(length(readLines(delined_fasta)), 2 * length(test_dss2))
})
