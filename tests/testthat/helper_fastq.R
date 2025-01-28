# generate some test fastq sequences

test_seq <- replicate(
  100,
  sample(c("A", "C", "G", "T"), 100, replace = TRUE),
  simplify = FALSE
) |>
  vapply(paste, "", collapse = "")

test_qual <- replicate(
  100,
  sample.int(40, 100, replace = TRUE),
  simplify = FALSE
) |>
  S4Vectors::List() |>
  Biostrings::PhredQuality()

test_bins <- c(3L, 12L, 26L, 38L)

test_qual_bin <- replicate(
   100,
   sample(test_bins, 100, replace = TRUE),
   simplify = FALSE
) |>
   S4Vectors::List() |>
   Biostrings::PhredQuality()

# pretend that these 100 are an unordered subsample of a total set of 1000,
# and label them with hex-encoded numbers
test_ind <- sample.int(1000, 100)
test_names1 <- format(as.hexmode(test_ind))

# alternately, just number them
test_names2 <- sprintf("read%03d", seq_len(100))

test_dss1 <- Biostrings::DNAStringSet(
  `names<-`(test_seq, test_names1)
)

test_qsdss1 <- Biostrings::QualityScaledDNAStringSet(
  `names<-`(test_seq, test_names1),
  test_qual
)

test_dss2 <- Biostrings::DNAStringSet(
  `names<-`(test_seq, test_names2)
)

test_qsdss2 <- Biostrings::QualityScaledDNAStringSet(
  `names<-`(test_seq, test_names2),
  test_qual
)

test_qsdss3 <- Biostrings::QualityScaledDNAStringSet(
   `names<-`(test_seq, test_names1),
   test_qual_bin
)

# 8 different random selections of ~80% of the full 100
# to represent different processing stages
selections <- replicate(8, rbinom(100, 1, 0.8), simplify = FALSE)

# write the initial files as fasta, fasta.gz, fastq and fastq.gz
raw1_fasta <- tempfile(fileext = ".fasta")
Biostrings::writeXStringSet(test_dss1, raw1_fasta)
raw1_fasta_gz <- tempfile(fileext = ".fasta.gz")
Biostrings::writeXStringSet(test_dss1, raw1_fasta_gz, compress = "gzip")
raw1_fastq <- tempfile(fileext = ".fastq")
Biostrings::writeQualityScaledXStringSet(test_qsdss1, raw1_fastq)
raw1_fastq_gz <- tempfile(fileext = ".fastq.gz")
Biostrings::writeQualityScaledXStringSet(test_qsdss1, raw1_fastq_gz, compress = "gzip")
raw2_fasta <- tempfile(fileext = ".fasta")
Biostrings::writeXStringSet(test_dss2, raw2_fasta)
raw2_fasta_gz <- tempfile(fileext = ".fasta.gz")
Biostrings::writeXStringSet(test_dss2, raw2_fasta_gz, compress = "gzip")
raw2_fastq <- tempfile(fileext = ".fastq")
Biostrings::writeQualityScaledXStringSet(test_qsdss2, raw2_fastq)
raw2_fastq_gz <- tempfile(fileext = ".fastq.gz")
Biostrings::writeQualityScaledXStringSet(test_qsdss2, raw2_fastq_gz, compress = "gzip")
raw3_fastq <- tempfile(fileext = ".fastq")
Biostrings::writeQualityScaledXStringSet(test_qsdss3, raw3_fastq)
raw3_fastq_gz <- tempfile(fileext = ".fastq.gz")
Biostrings::writeQualityScaledXStringSet(test_qsdss3, raw3_fastq_gz, compress = "gzip")

# write the "stage output" files as fastq and fastq.gz
stage1_fastq <- character()
stage1_fastq_gz <- character()
stage2_fastq <- character()
stage2_fastq_gz <- character()

# also calculate the "correct" bit pattern for presence/absence of the sequences
# at each stage
test_flags <- rep(0, 100)

for (i in seq_along(selections)) {
  test_flags <- test_flags + bitwShiftL(selections[[i]], i - 1)
  stage1_fastq[i] <- tempfile(fileext = ".fastq")
  stage1_fastq_gz[i] <- tempfile(fileext = ".fastq.gz")
  stage2_fastq[i] <- tempfile(fileext = ".fastq")
  stage2_fastq_gz[i] <- tempfile(fileext = ".fastq.gz")
  Biostrings::writeQualityScaledXStringSet(
    test_qsdss1[as.logical(selections[[i]])],
    stage1_fastq[i]
  )
  Biostrings::writeQualityScaledXStringSet(
    test_qsdss1[as.logical(selections[[i]])],
    stage1_fastq_gz[i],
    compress = "gzip"
  )
  Biostrings::writeQualityScaledXStringSet(
    test_qsdss2[as.logical(selections[[i]])],
    stage2_fastq[i]
  )
  Biostrings::writeQualityScaledXStringSet(
    test_qsdss2[as.logical(selections[[i]])],
    stage2_fastq_gz[i],
    compress = "gzip"
  )
}

test_flags <- as.raw(test_flags)
