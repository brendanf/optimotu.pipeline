
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

testthat::test_that("fastq_pair_sample_fraction works", {
  subsamp9_R1_fastq_gz <- withr::local_tempfile(fileext = ".fastq.gz")
  subsamp9_R2_fastq_gz <- withr::local_tempfile(fileext = ".fastq.gz")

  testthat::expect_equal(
    fastq_pair_sample_fraction(
      file_R1 = raw1_fastq_gz,
      file_R2 = raw1_fastq_R2_gz,
      numerator = 1L,
      denominator = 10L,
      output_R1 = subsamp9_R1_fastq_gz,
      output_R2 = subsamp9_R2_fastq_gz
    ),
    c(subsamp9_R1_fastq_gz, subsamp9_R2_fastq_gz)
  )
  testthat::expect_no_error(
    subsamp9_qsdss_R1 <- Biostrings::readQualityScaledDNAStringSet(subsamp9_R1_fastq_gz)
  )
  testthat::expect_no_error(

    subsamp9_qsdss_R2 <- Biostrings::readQualityScaledDNAStringSet(subsamp9_R2_fastq_gz)
  )
  testthat::expect_true(all(
    test_qsdss1[names(subsamp9_qsdss_R1)] == subsamp9_qsdss_R1
  ))
  testthat::expect_true(all(
    test_qsdss1_R2[names(subsamp9_qsdss_R2)] == subsamp9_qsdss_R2
  ))
})

testthat::test_that("fastq_pair_sample_fraction with rename works", {
  subsamp10_fastq_gz_R1 <- withr::local_tempfile(fileext = ".fastq.gz")
  subsamp10_fastq_gz_R2 <- withr::local_tempfile(fileext = ".fastq.gz")

  testthat::expect_equal(
    optimotu.pipeline:::fastq_pair_sample_fraction(
      file_R1 = raw1_fastq_gz,
      file_R2 = raw1_fastq_R2_gz,
      numerator = 1L,
      denominator = 10L,
      output_R1 = subsamp10_fastq_gz_R1,
      output_R2 = subsamp10_fastq_gz_R2,
      rename = TRUE
    ),
    c(subsamp10_fastq_gz_R1, subsamp10_fastq_gz_R2)
  )

  testthat::expect_no_error(
    subsamp10_qsdss_R1 <- Biostrings::readQualityScaledDNAStringSet(subsamp10_fastq_gz_R1)
  )
  testthat::expect_no_error(
    subsamp10_qsdss_R2 <- Biostrings::readQualityScaledDNAStringSet(subsamp10_fastq_gz_R2)
  )

  testthat::expect_equal(
    as.character(test_qsdss1[strtoi(names(subsamp10_qsdss_R1), 16)], use.names = FALSE),
    as.character(subsamp10_qsdss_R1, use.names = FALSE)
  )
  testthat::expect_equal(
    as.character(test_qsdss1_R2[strtoi(names(subsamp10_qsdss_R2), 16)], use.names = FALSE),
    as.character(subsamp10_qsdss_R2, use.names = FALSE)
  )
  testthat::expect_identical(
    names(subsamp10_qsdss_R1),
    names(subsamp10_qsdss_R2)
  )
})

testthat::test_that("fastq_pair_sample_number works", {
  subsamp11_fastq_gz_R1 <- withr::local_tempfile(fileext = ".fastq.gz")
  subsamp11_fastq_gz_R2 <- withr::local_tempfile(fileext = ".fastq.gz")

  testthat::expect_equal(
    fastq_pair_sample_number(
      file_R1 = raw1_fastq_gz,
      file_R2 = raw1_fastq_R2_gz,
      number = 34L,
      output_R1 = subsamp11_fastq_gz_R1,
      output_R2 = subsamp11_fastq_gz_R2
    ),
    c(subsamp11_fastq_gz_R1, subsamp11_fastq_gz_R2)
  )
  testthat::expect_no_error(
    subsamp11_qsdss_R1 <- Biostrings::readQualityScaledDNAStringSet(subsamp11_fastq_gz_R1)
  )
  testthat::expect_no_error(
    subsamp11_qsdss_R2 <- Biostrings::readQualityScaledDNAStringSet(subsamp11_fastq_gz_R2)
  )
  testthat::expect_true(all(
    test_qsdss1[names(subsamp11_qsdss_R1)] == subsamp11_qsdss_R1
  ))
  testthat::expect_true(all(
    test_qsdss1_R2[names(subsamp11_qsdss_R2)] == subsamp11_qsdss_R2
  ))
  testthat::expect_identical(
    names(subsamp11_qsdss_R1),
    names(subsamp11_qsdss_R2)
  )
  testthat::expect_length(subsamp11_qsdss_R1, 34)
  testthat::expect_length(subsamp11_qsdss_R2, 34)
})

testthat::test_that("fastq_pair_sample_number with rename works", {
  subsamp12_fastq_gz_R1 <- withr::local_tempfile(fileext = ".fastq.gz")
  subsamp12_fastq_gz_R2 <- withr::local_tempfile(fileext = ".fastq.gz")

  testthat::expect_equal(
    fastq_pair_sample_number(
      file_R1 = raw1_fastq_gz,
      file_R2 = raw1_fastq_R2_gz,
      number = 99L,
      output_R1 = subsamp12_fastq_gz_R1,
      output_R2 = subsamp12_fastq_gz_R2,
      rename = TRUE
    ),
    c(subsamp12_fastq_gz_R1, subsamp12_fastq_gz_R2)
  )
  testthat::expect_no_error(
    subsamp12_qsdss_R1 <- Biostrings::readQualityScaledDNAStringSet(subsamp12_fastq_gz_R1)
  )
  testthat::expect_no_error(
    subsamp12_qsdss_R2 <- Biostrings::readQualityScaledDNAStringSet(subsamp12_fastq_gz_R2)
  )
  testthat::expect_true(all(
    test_qsdss1[strtoi(names(subsamp12_qsdss_R1), 16)] == subsamp12_qsdss_R1
  ))
  testthat::expect_true(all(
    test_qsdss1_R2[strtoi(names(subsamp12_qsdss_R2), 16)] == subsamp12_qsdss_R2
  ))
  testthat::expect_identical(
    names(subsamp12_qsdss_R1),
    names(subsamp12_qsdss_R2)
  )
  testthat::expect_length(subsamp12_qsdss_R1, 99)
  testthat::expect_length(subsamp12_qsdss_R2, 99)
})

testthat::test_that("fastq_pair_sample_number_multiple works", {
  subsamp13_fastq_gz_R1 <- withr::local_tempfile(fileext = ".fastq.gz")
  subsamp13_fastq_gz_R2 <- withr::local_tempfile(fileext = ".fastq.gz")
  subsamp14_fastq_gz_R1 <- withr::local_tempfile(fileext = ".fastq.gz")
  subsamp14_fastq_gz_R2 <- withr::local_tempfile(fileext = ".fastq.gz")

  testthat::expect_equal(
    fastq_pair_sample_number_multiple(
      file_R1 = raw1_fastq_gz,
      file_R2 = raw1_fastq_R2_gz,
      number = c(17, 81),
      output_R1 = c(subsamp13_fastq_gz_R1, subsamp14_fastq_gz_R1),
      output_R2 = c(subsamp13_fastq_gz_R2, subsamp14_fastq_gz_R2)
    ),
    list(R1 = c(subsamp13_fastq_gz_R1, subsamp14_fastq_gz_R1), R2 = c(subsamp13_fastq_gz_R2, subsamp14_fastq_gz_R2))
  )
  testthat::expect_no_error(
    subsamp13_qsdss_R1 <- Biostrings::readQualityScaledDNAStringSet(subsamp13_fastq_gz_R1)
  )
  testthat::expect_no_error(
    subsamp13_qsdss_R2 <- Biostrings::readQualityScaledDNAStringSet(subsamp13_fastq_gz_R2)
  )
  testthat::expect_no_error(
    subsamp14_qsdss_R1 <- Biostrings::readQualityScaledDNAStringSet(subsamp14_fastq_gz_R1)
  )
  testthat::expect_no_error(
    subsamp14_qsdss_R2 <- Biostrings::readQualityScaledDNAStringSet(subsamp14_fastq_gz_R2)
  )
  testthat::expect_true(all(test_qsdss1[names(subsamp13_qsdss_R1)] == subsamp13_qsdss_R1))
  testthat::expect_true(all(test_qsdss1_R2[names(subsamp13_qsdss_R2)] == subsamp13_qsdss_R2))
  testthat::expect_true(all(test_qsdss1[names(subsamp14_qsdss_R1)] == subsamp14_qsdss_R1))
  testthat::expect_true(all(test_qsdss1_R2[names(subsamp14_qsdss_R2)] == subsamp14_qsdss_R2))
  testthat::expect_identical(
    names(subsamp13_qsdss_R1),
    names(subsamp13_qsdss_R2)
  )
  testthat::expect_identical(
    names(subsamp14_qsdss_R1),
    names(subsamp14_qsdss_R2)
  )
  testthat::expect_contains(
    names(subsamp14_qsdss_R1),
    names(subsamp13_qsdss_R1)
  )
  testthat::expect_contains(
    names(subsamp14_qsdss_R2),
    names(subsamp13_qsdss_R2)
  )
  testthat::expect_length(subsamp13_qsdss_R1, 17)
  testthat::expect_length(subsamp13_qsdss_R2, 17)
  testthat::expect_length(subsamp14_qsdss_R1, 81)
  testthat::expect_length(subsamp14_qsdss_R2, 81)
})

testthat::test_that("fastq_pair_sample_number_multiple with rename works", {
  subsamp15_fastq_gz_R1 <- withr::local_tempfile(fileext = ".fastq.gz")
  subsamp15_fastq_gz_R2 <- withr::local_tempfile(fileext = ".fastq.gz")
  subsamp16_fastq_gz_R1 <- withr::local_tempfile(fileext = ".fastq.gz")
  subsamp16_fastq_gz_R2 <- withr::local_tempfile(fileext = ".fastq.gz")

  testthat::expect_equal(
    fastq_pair_sample_number_multiple(
      file_R1 = raw1_fastq_gz,
      file_R2 = raw1_fastq_R2_gz,
      number = c(77, 4),
      output_R1 = c(subsamp15_fastq_gz_R1, subsamp16_fastq_gz_R1),
      output_R2 = c(subsamp15_fastq_gz_R2, subsamp16_fastq_gz_R2),
      rename = TRUE
    ),
    list(R1 = c(subsamp15_fastq_gz_R1, subsamp16_fastq_gz_R1), R2 = c(subsamp15_fastq_gz_R2, subsamp16_fastq_gz_R2))
  )
  testthat::expect_no_error(
    subsamp15_qsdss_R1 <- Biostrings::readQualityScaledDNAStringSet(subsamp15_fastq_gz_R1)
  )
  testthat::expect_no_error(
    subsamp15_qsdss_R2 <- Biostrings::readQualityScaledDNAStringSet(subsamp15_fastq_gz_R2)
  )
  testthat::expect_no_error(
    subsamp16_qsdss_R1 <- Biostrings::readQualityScaledDNAStringSet(subsamp16_fastq_gz_R1)
  )
  testthat::expect_no_error(
    subsamp16_qsdss_R2 <- Biostrings::readQualityScaledDNAStringSet(subsamp16_fastq_gz_R2)
  )
  testthat::expect_true(all(test_qsdss1[strtoi(names(subsamp15_qsdss_R1), 16)] == subsamp15_qsdss_R1))
  testthat::expect_true(all(test_qsdss1_R2[strtoi(names(subsamp15_qsdss_R2), 16)] == subsamp15_qsdss_R2))
  testthat::expect_true(all(test_qsdss1[strtoi(names(subsamp16_qsdss_R1), 16)] == subsamp16_qsdss_R1))
  testthat::expect_true(all(test_qsdss1_R2[strtoi(names(subsamp16_qsdss_R2), 16)] == subsamp16_qsdss_R2))
  testthat::expect_identical(
    names(subsamp15_qsdss_R1),
    names(subsamp15_qsdss_R2)
  )
  testthat::expect_identical(
    names(subsamp16_qsdss_R1),
    names(subsamp16_qsdss_R2)
  )
  testthat::expect_contains(
    names(subsamp15_qsdss_R1),
    names(subsamp16_qsdss_R1)
  )
  testthat::expect_contains(
    names(subsamp15_qsdss_R2),
    names(subsamp16_qsdss_R2)
  )
  testthat::expect_length(subsamp15_qsdss_R1, 77)
  testthat::expect_length(subsamp15_qsdss_R2, 77)
  testthat::expect_length(subsamp16_qsdss_R1, 4)
  testthat::expect_length(subsamp16_qsdss_R2, 4)
})

testthat::test_that("fastq_sample works", {
  subsamp17_fastq_gz <- withr::local_tempfile(fileext = ".fastq.gz")
  shuf <- sample(100L)
  testthat::expect_equal(
    fastq_sample(
      infile = raw1_fastq_gz,
      outfile = subsamp17_fastq_gz,
      n = 10L,
      sample = shuf
    ),
    subsamp17_fastq_gz
  )
  testthat::expect_no_error(
    subsamp17_qsdss <- Biostrings::readQualityScaledDNAStringSet(subsamp17_fastq_gz)
  )
  testthat::expect_length(subsamp17_qsdss, 10)
  testthat::expect_equal(
    as.character(test_qsdss1[names(subsamp17_qsdss)], use.names = FALSE),
    as.character(subsamp17_qsdss, use.names = FALSE)
  )
  testthat::expect_identical(
    names(subsamp17_qsdss),
    names(test_qsdss1[which(shuf <= 10L)])
  )
})
