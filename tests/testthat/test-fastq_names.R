test_seq <- c("ACGACG", "TGTGA")
test_names <- c("first_sequence", "test")
names(test_seq) <- test_names
test_qual <- Biostrings::PhredQuality(c("$&89()", "#1/#9"))
test_qsdss <- Biostrings::QualityScaledDNAStringSet(test_seq, test_qual)
test_fastq <- tempfile("test", fileext = ".fastq")
test_fastq_gz <- tempfile("test", fileext = ".fastq.gz")
Biostrings::writeQualityScaledXStringSet(test_qsdss, test_fastq)
Biostrings::writeQualityScaledXStringSet(test_qsdss, test_fastq_gz, compress = "gzip")

test_that("fastq_names works", {
  expect_equal(fastq_names(test_fastq), test_names)
  expect_equal(fastq_names(test_fastq_gz), test_names)
})
