
testthat::test_that("fastq_sample_fraction works", {
  subsamp1_fastq_gz <- withr::local_tempfile(fileext = ".fastq.gz")
  testthat::expect_equal(
    optimotu.pipeline:::fastq_sample_fraction(
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
    optimotu.pipeline:::fastq_sample_fraction(
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

testthat::test_that("fastq_sample_number works", {
  subsamp5_fastq_gz <- withr::local_tempfile(fileext = ".fastq.gz")
  testthat::expect_equal(
    optimotu.pipeline:::fastq_sample_number(
      file = raw1_fastq_gz,
      number = 7L,
      output = subsamp5_fastq_gz
    ),
    subsamp5_fastq_gz
  )
  testthat::expect_no_error(
    subsamp5_qsdss <- Biostrings::readQualityScaledDNAStringSet(subsamp5_fastq_gz)
  )
  testthat::expect_true(all(
    test_qsdss1[names(subsamp5_qsdss)] == subsamp5_qsdss
  ))
  testthat::expect_length(subsamp5_qsdss, 7)
})

testthat::test_that("fastq_sample_number with rename works", {
  subsamp6_fastq_gz <- withr::local_tempfile(fileext = ".fastq.gz")

  testthat::expect_equal(
    optimotu.pipeline:::fastq_sample_number(
      file = raw1_fastq_gz,
      number = 61,
      output = subsamp6_fastq_gz,
      rename = TRUE
    ),
    subsamp6_fastq_gz
  )

  testthat::expect_no_error(
    subsamp6_qsdss <- Biostrings::readQualityScaledDNAStringSet(subsamp6_fastq_gz)
  )

  testthat::expect_equal(
    as.character(test_qsdss1[strtoi(names(subsamp6_qsdss), 16)], use.names = FALSE),
    as.character(subsamp6_qsdss, use.names = FALSE)
  )
  testthat::expect_length(subsamp6_qsdss, 61)
})

testthat::test_that("fastq_sample_multiple works", {
  subsamp3_fastq_gz <- withr::local_tempfile(fileext = ".fastq.gz")
  subsamp4_fastq_gz <- withr::local_tempfile(fileext = ".fastq.gz")
  testthat::expect_equal(
    optimotu.pipeline:::fastq_sample_fraction_multiple(
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

testthat::test_that("fastq_sample_number_multiple works", {
  subsamp7_fastq_gz <- withr::local_tempfile(fileext = ".fastq.gz")
  subsamp8_fastq_gz <- withr::local_tempfile(fileext = ".fastq.gz")
  testthat::expect_equal(
    optimotu.pipeline:::fastq_sample_number_multiple(
      file = raw1_fastq_gz,
      number = c(20, 50),
      output = c(subsamp7_fastq_gz, subsamp8_fastq_gz)
    ),
    c(subsamp7_fastq_gz, subsamp8_fastq_gz)
  )
  testthat::expect_no_error(
    subsamp7_qsdss <- Biostrings::readQualityScaledDNAStringSet(subsamp7_fastq_gz)
  )
  testthat::expect_no_error(
    subsamp8_qsdss <- Biostrings::readQualityScaledDNAStringSet(subsamp8_fastq_gz)
  )
  testthat::expect_true(all(
    test_qsdss1[names(subsamp7_qsdss)] == subsamp7_qsdss
  ))
  testthat::expect_true(all(
    test_qsdss1[names(subsamp8_qsdss)] == subsamp8_qsdss
  ))
  testthat::expect_contains(names(subsamp8_qsdss), names(subsamp7_qsdss))
})
