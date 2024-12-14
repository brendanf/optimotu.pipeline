test_that("fastq_qual_bins works", {
   expect_equal(fastq_qual_bins(raw1_fastq), seq.int(40))
   expect_equal(fastq_qual_bins(raw1_fastq_gz), seq.int(40))
   expect_equal(fastq_qual_bins(c(raw1_fastq, raw1_fastq_gz)), seq.int(40))
   expect_equal(fastq_qual_bins(raw3_fastq), test_bins)
   expect_equal(fastq_qual_bins(raw3_fastq_gz), test_bins)
   expect_equal(fastq_qual_bins(c(raw3_fastq, raw3_fastq_gz)), test_bins)
})
