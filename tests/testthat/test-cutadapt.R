skip_if_no_cutadapt <- function() {
  tc_skip_if_no_exec("cutadapt")
}

write_test_fastq <- function(seq, id, file) {
  dss <- Biostrings::DNAStringSet(setNames(seq, id))
  q <- Biostrings::PhredQuality(
    as.integer(rep(40L, Biostrings::nchar(seq)))
  )
  Biostrings::writeQualityScaledXStringSet(
    Biostrings::QualityScaledDNAStringSet(dss, q),
    file
  )
  file
}

test_that("cutadapt option constructors expose compression_level", {
  single <- cutadapt_options()
  paired <- cutadapt_paired_options()

  expect_s3_class(single, "cutadapt_options")
  expect_s3_class(paired, "cutadapt_paired_options")
  expect_equal(single$compression_level, 6L)
  expect_equal(paired$compression_level, 6L)
  expect_true("compression_level" %in% cutadapt_option_names)
  expect_true("compression_level" %in% cutadapt_paired_option_names)
})

test_that("cutadapt option constructors validate compression_level", {
  expect_error(cutadapt_options(compression_level = 0L), "Must be >= 1")
  expect_error(
    cutadapt_options(compression_level = 10L),
    "Must be element of set"
  )
  expect_error(cutadapt_paired_options(compression_level = 0L), "Must be >= 1")
  expect_error(
    cutadapt_paired_options(compression_level = 10L),
    "Must be element of set"
  )
})

test_that("update.cutadapt_options merges new values", {
  opts <- cutadapt_options(max_err = 0.2, compression_level = 6L)
  updated <- stats::update(
    opts,
    list(compression_level = 3L, min_length = 50L)
  )
  expect_s3_class(updated, "cutadapt_options")
  expect_equal(updated$compression_level, 3L)
  expect_equal(updated$min_length, 50L)
  expect_equal(updated$max_err, 0.2)
})

test_that("update.cutadapt_paired_options merges new values", {
  opts <- cutadapt_paired_options(truncQ_R1 = 2L, compression_level = 6L)
  updated <- stats::update(
    opts,
    list(compression_level = 1L, truncQ_R2 = 5L)
  )
  expect_s3_class(updated, "cutadapt_paired_options")
  expect_equal(updated$compression_level, 1L)
  expect_equal(updated$truncQ_R2, 5L)
  expect_equal(updated$truncQ_R1, 2L)
})

test_that("update.cutadapt_* rejects inconsistent batch options", {
  expect_error(
    stats::update(
      cutadapt_options(),
      data.frame(compression_level = c(1L, 2L))
    ),
    "must be the same for all samples"
  )
  expect_error(
    stats::update(
      cutadapt_paired_options(),
      data.frame(compression_level = c(1L, 2L))
    ),
    "must be the same for all samples"
  )
})

test_that("cutadapt_filter_trim trims a single-end FASTA read", {
  cutadapt <- skip_if_no_cutadapt()
  primer <- "ACGTACGT"
  body <- "AAAATTTTGGGGCCCC"
  fa <- withr::local_tempfile(fileext = ".fasta")
  trim <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(">read1", paste0(primer, body)), fa)

  out <- cutadapt_filter_trim(
    file = fa,
    primer = primer,
    trim = trim,
    options = cutadapt_options(
      action = "trim",
      discard_untrimmed = FALSE,
      compression_level = 6L
    ),
    ncpu = 1L,
    cutadapt = cutadapt
  )

  expect_equal(out, trim)
  expect_equal(unname(as.character(Biostrings::readBStringSet(trim))), body)
})

test_that("cutadapt_paired_filter_trim trims paired FASTQ reads", {
  cutadapt <- skip_if_no_cutadapt()
  primer_R1 <- "ACGTACGT"
  primer_R2 <- "TGCATGCA"
  body_R1 <- "AAAATTTTGGGGCCCC"
  body_R2 <- "CCCCGGGGTTTTAAAA"
  r1 <- withr::local_tempfile(fileext = ".fastq")
  r2 <- withr::local_tempfile(fileext = ".fastq")
  trim_R1 <- withr::local_tempfile(fileext = ".fastq.gz")
  trim_R2 <- withr::local_tempfile(fileext = ".fastq.gz")
  write_test_fastq(paste0(primer_R1, body_R1), "read1", r1)
  write_test_fastq(paste0(primer_R2, body_R2), "read1", r2)

  out <- cutadapt_paired_filter_trim(
    file_R1 = r1,
    file_R2 = r2,
    primer_R1 = primer_R1,
    primer_R2 = primer_R2,
    trim_R1 = trim_R1,
    trim_R2 = trim_R2,
    options = cutadapt_paired_options(
      action = "trim",
      discard_untrimmed = FALSE,
      compression_level = 3L
    ),
    ncpu = 1L,
    cutadapt = cutadapt
  )

  expect_equal(out, c(trim_R1, trim_R2))
  expect_equal(
    unname(as.character(Biostrings::readQualityScaledDNAStringSet(trim_R1))),
    body_R1
  )
  expect_equal(
    unname(as.character(Biostrings::readQualityScaledDNAStringSet(trim_R2))),
    body_R2
  )
})

test_that("trim_primer returns trimmed sequences", {
  skip_if_no_cutadapt()
  primer <- "ACGTACGT"
  body <- "AAAATTTTGGGGCCCC"
  out <- trim_primer(
    c(read1 = paste0(primer, body)),
    primer = primer,
    options = cutadapt_options(
      action = "trim",
      discard_untrimmed = FALSE
    ),
    ncpu = 1L
  )

  expect_s3_class(out, "tbl_df")
  expect_equal(out$seq_id, "read1")
  expect_equal(out$seq, body)
})
