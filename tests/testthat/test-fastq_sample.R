
testthat::test_that("fastq_sample works", {
  subsamp1_fastq_gz <- withr::local_tempfile(fileext = ".fastq.gz")
  testthat::expect_equal(
    optimotu.pipeline:::fastq_sample(
      file = raw1_fastq_gz,
      numerator = 1L,
      denominator = 10L,
      output = subsamp1_fastq_gz
    ),
    subsamp1_fastq_gz
  )
  testthat::expect_no_error(
    subsamp_qsdss <- Biostrings::readQualityScaledDNAStringSet(subsamp1_fastq_gz)
  )
  testthat::expect_true(all(
    test_qsdss1[names(subsamp_qsdss)] == subsamp_qsdss
  ))
})

testthat::test_that("fastq_sample with rename works", {
  subsamp2_fastq_gz <- withr::local_tempfile(fileext = ".fastq.gz")

  testthat::expect_equal(
    optimotu.pipeline:::fastq_sample(
      file = raw1_fastq_gz,
      numerator = 1L,
      denominator = 10L,
      output = subsamp2_fastq_gz,
      rename = TRUE
    ),
    subsamp2_fastq_gz
  )

  testthat::expect_no_error(
    subsamp2_qsdss <- Biostrings::readQualityScaledDNAStringSet(subsamp2_fastq_gz)
  )

  testthat::expect_equal(
    as.character(test_qsdss1[strtoi(names(subsamp2_qsdss), 16)], use.names = FALSE),
    as.character(subsamp2_qsdss, use.names = FALSE)
  )
})

testthat::test_that("fastq_sample_multiple works", {
  subsamp3_fastq_gz <- withr::local_tempfile(fileext = ".fastq.gz")
  subsamp4_fastq_gz <- withr::local_tempfile(fileext = ".fastq.gz")
  testthat::expect_equal(
    optimotu.pipeline:::fastq_sample_multiple(
      file = raw1_fastq_gz,
      numerator = c(2L, 3L),
      denominator = 11L,
      output = c(subsamp3_fastq_gz, subsamp4_fastq_gz)
    ),
    c(subsamp3_fastq_gz, subsamp4_fastq_gz)
  )
  testthat::expect_no_error(
    subsamp3_qsdss <- Biostrings::readQualityScaledDNAStringSet(subsamp3_fastq_gz)
  )
  testthat::expect_no_error(
    subsamp4_qsdss <- Biostrings::readQualityScaledDNAStringSet(subsamp4_fastq_gz)
  )
  testthat::expect_true(all(
    test_qsdss1[names(subsamp3_qsdss)] == subsamp3_qsdss
  ))
  testthat::expect_true(all(
    test_qsdss1[names(subsamp4_qsdss)] == subsamp4_qsdss
  ))
  testthat::expect_contains(names(subsamp4_qsdss), names(subsamp3_qsdss))
})
