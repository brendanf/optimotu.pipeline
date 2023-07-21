test_that("fastq_names works", {
  expect_equal(fastq_names(raw1_fastq), test_names1)
  expect_equal(fastq_names(raw1_fastq_gz), test_names1)
})
