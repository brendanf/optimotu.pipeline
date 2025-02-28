test_splits <- c(1, 4, 7)

test_that("fastq_split works fastq->fastq", {
  for (n in test_splits) {
    splitfiles <- replicate(n, withr::local_tempfile(fileext=".fastq"))
    expect_no_error(
      splits <- fastq_split(raw2_fastq, splitfiles)
    )
    expect_equal(splitfiles, splits)
    for (i in 1:n) {
      expect_true(file.exists(splitfiles[i]))
      expect_true(file.size(splitfiles[i]) > 0)
      expect_true(
        all(
          test_qsdss2[seq(i, length(test_qsdss2), n)] ==
            Biostrings::readQualityScaledDNAStringSet(splitfiles[i])
        )
      )
    }
    recombine_file <- withr::local_tempfile(fileext=".fastq")
    expect_no_error(
      recombine <- fastq_combine(splitfiles, recombine_file)
    )
    expect_identical(recombine, recombine_file)
    expect_true(file.exists(recombine_file))
    expect_true(file.size(recombine_file) > 0)
    expect_true(
      all(
        test_qsdss2 ==
          Biostrings::readQualityScaledDNAStringSet(recombine_file)
      )
    )

  }
})

test_that("fastq_split and combine works fastq.gz->fastq->fastq.gz ", {
  for (n in 1:10) {
    splitfiles <- replicate(n, withr::local_tempfile(fileext=".fastq"))
    expect_no_error(
      splits <- fastq_split(raw2_fastq_gz, splitfiles)
    )
    expect_equal(splitfiles, splits)
    for (i in 1:n) {
      expect_true(file.exists(splitfiles[i]))
      expect_true(file.size(splitfiles[i]) > 0)
      expect_true(
        all(
          test_qsdss2[seq(i, length(test_qsdss2), n)] ==
            Biostrings::readQualityScaledDNAStringSet(splitfiles[i])
        )
      )
    }
    recombine_file <- withr::local_tempfile(fileext=".fastq.gz")
    expect_no_error(
      recombine <- fastq_combine(splitfiles, recombine_file, compress = TRUE)
    )
    expect_identical(recombine, recombine_file)
    expect_true(file.exists(recombine_file))
    expect_true(zlib::validate_gzip_file(recombine_file))
    expect_true(file.size(recombine_file) > 0)
    expect_true(
      all(
        test_qsdss2 ==
          Biostrings::readQualityScaledDNAStringSet(recombine_file)
      )
    )
  }
})

test_that("fastq_split and combine works fastq.gz->fastq.gz->fastq.gz", {
  for (n in test_splits) {
    splitfiles <- replicate(n, withr::local_tempfile(fileext=".fastq.gz"))
    expect_no_error(
      splits <- fastq_split(raw2_fastq_gz, splitfiles, compress=TRUE)
    )
    expect_equal(splitfiles, splits)
    for (i in 1:n) {
      expect_true(file.exists(splitfiles[i]))
      expect_true(zlib::validate_gzip_file(splitfiles[i]))
      expect_true(file.size(splitfiles[i]) > 0)
      expect_true(
        all(
          test_qsdss2[seq(i, length(test_qsdss2), n)] ==
            Biostrings::readQualityScaledDNAStringSet(splitfiles[i])
        )
      )
    }
    recombine_file <- withr::local_tempfile(fileext=".fastq.gz")
    expect_no_error(
      recombine <- fastq_combine(splitfiles, recombine_file, compress=TRUE)
    )
    expect_identical(recombine, recombine_file)
    expect_true(file.exists(recombine_file))
    expect_true(zlib::validate_gzip_file(recombine_file))
    expect_true(file.size(recombine_file) > 0)
    expect_true(
      all(
        test_qsdss2 ==
          Biostrings::readQualityScaledDNAStringSet(recombine_file)
      )
    )
  }
})

test_that("fastq_split works fastq->fastq.gz->fastq", {
  for (n in test_splits) {
    splitfiles <- replicate(n, withr::local_tempfile(fileext=".fastq.gz"))
    expect_no_error(
      splits <- fastq_split(raw2_fastq, splitfiles, compress=TRUE)
    )
    expect_equal(splitfiles, splits)
    for (i in 1:n) {
      expect_true(file.exists(splitfiles[i]))
      expect_true(zlib::validate_gzip_file(splitfiles[i]))
      expect_true(file.size(splitfiles[i]) > 0)
      expect_true(
        all(
          test_qsdss2[seq(i, length(test_qsdss2), n)] ==
            Biostrings::readQualityScaledDNAStringSet(splitfiles[i])
        )
      )
    }
    recombine_file <- withr::local_tempfile(fileext=".fastq")
    expect_no_error(
      recombine <- fastq_combine(splitfiles, recombine_file)
    )
    expect_identical(recombine, recombine_file)
    expect_true(file.exists(recombine_file))
    expect_true(file.size(recombine_file) > 0)
    expect_true(
      all(
        test_qsdss2 ==
          Biostrings::readQualityScaledDNAStringSet(recombine_file)
      )
    )
  }
})

test_that("fasta_split and combine works fasta->fasta->fasta", {
  for (n in test_splits) {
    splitfiles <- replicate(n, withr::local_tempfile(fileext=".fasta"))
    expect_no_error(
      splits <- fasta_split(raw2_fasta, splitfiles)
    )
    expect_equal(splitfiles, splits)
    for (i in 1:n) {
      expect_true(file.exists(splitfiles[i]))
      expect_true(file.size(splitfiles[i]) > 0)
      expect_true(
        all(
          test_dss2[seq(i, length(test_dss2), n)] ==
            Biostrings::readDNAStringSet(splitfiles[i])
        )
      )
    }
    recombine_file <- withr::local_tempfile(fileext=".fasta")
    expect_no_error(
      recombine <- fasta_combine(splitfiles, recombine_file)
    )
    expect_identical(recombine, recombine_file)
    expect_true(file.exists(recombine_file))
    expect_true(file.size(recombine_file) > 0)
    expect_true(
      all(
        test_dss2 ==
          Biostrings::readDNAStringSet(recombine_file)
      )
    )
  }
})

test_that("fasta_split works fasta.gz->fasta->fasta.gz", {
  for (n in test_splits) {
    splitfiles <- replicate(n, withr::local_tempfile(fileext=".fasta"))
    expect_no_error(
      splits <- fasta_split(raw2_fasta_gz, splitfiles)
    )
    expect_equal(splitfiles, splits)
    for (i in 1:n) {
      expect_true(file.exists(splitfiles[i]))
      expect_true(file.size(splitfiles[i]) > 0)
      expect_true(
        all(
          test_dss2[seq(i, length(test_dss2), n)] ==
            Biostrings::readDNAStringSet(splitfiles[i])
        )
      )
    }
    recombine_file <- withr::local_tempfile(fileext=".fasta.gz")
    expect_no_error(
      recombine <- fasta_combine(splitfiles, recombine_file, compress=TRUE)
    )
    expect_identical(recombine, recombine_file)
    expect_true(file.exists(recombine_file))
    expect_true(zlib::validate_gzip_file(recombine_file))
    expect_true(file.size(recombine_file) > 0)
    expect_true(
      all(
        test_dss2 ==
          Biostrings::readDNAStringSet(recombine_file)
      )
    )
  }
})

test_that("fasta_split works fasta.gz->fasta.gz->fasta.gz", {
  for (n in test_splits) {
    splitfiles <- replicate(n, withr::local_tempfile(fileext=".fasta.gz"))
    expect_no_error(
      splits <- fasta_split(raw2_fasta_gz, splitfiles, compress=TRUE)
    )
    expect_equal(splitfiles, splits)
    for (i in 1:n) {
      expect_true(file.exists(splitfiles[i]))
      expect_true(zlib::validate_gzip_file(splitfiles[i]))
      expect_true(file.size(splitfiles[i]) > 0)
      expect_true(
        all(
          test_dss2[seq(i, length(test_dss2), n)] ==
            Biostrings::readDNAStringSet(splitfiles[i])
        )
      )
    }
    recombine_file <- withr::local_tempfile(fileext=".fasta.gz")
    expect_no_error(
      recombine <- fasta_combine(splitfiles, recombine_file, compress=TRUE)
    )
    expect_identical(recombine, recombine_file)
    expect_true(file.exists(recombine_file))
    expect_true(zlib::validate_gzip_file(recombine_file))
    expect_true(file.size(recombine_file) > 0)
    expect_true(
      all(
        test_dss2 ==
          Biostrings::readDNAStringSet(recombine_file)
      )
    )
  }
})

test_that("fasta_split works fasta->fasta.gz->fasta", {
  for (n in test_splits) {
    splitfiles <- replicate(n, withr::local_tempfile(fileext=".fasta.gz"))
    expect_no_error(
      splits <- fasta_split(raw2_fasta, splitfiles, compress=TRUE)
    )
    expect_equal(splitfiles, splits)
    for (i in 1:n) {
      expect_true(file.exists(splitfiles[i]))
      expect_true(zlib::validate_gzip_file(splitfiles[i]))
      expect_true(file.size(splitfiles[i]) > 0)
      expect_true(
        all(
          test_dss2[seq(i, length(test_dss2), n)] ==
            Biostrings::readDNAStringSet(splitfiles[i])
        )
      )
    }
    recombine_file <- withr::local_tempfile(fileext=".fasta")
    expect_no_error(
      recombine <- fasta_combine(splitfiles, recombine_file)
    )
    expect_identical(recombine, recombine_file)
    expect_true(file.exists(recombine_file))
    expect_true(file.size(recombine_file) > 0)
    expect_true(
      all(
        test_dss2 ==
          Biostrings::readDNAStringSet(recombine_file)
      )
    )
  }
})
