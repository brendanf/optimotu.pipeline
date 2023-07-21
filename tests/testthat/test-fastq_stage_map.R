test_result1 <- data.frame(
  read_idx = test_ind,
  flags = test_flags
)

test_result2 <- data.frame(
  read_idx = seq_along(test_names1),
  flags = test_flags
)

test_that("fastq_stage_map compressed with numbers", {
  expect_equal(
    fastq_stage_map(raw1_fastq_gz, stage1_fastq_gz),
    test_result1
  )
})

test_that("fastq_stage_map uncompressed with numbers", {
  expect_equal(
    fastq_stage_map(raw1_fastq, stage1_fastq),
    test_result1
  )
})

test_that("fastq_stage_map compressed with names", {
  expect_equal(
    fastq_stage_map(raw2_fastq_gz, stage2_fastq_gz),
    test_result2
  )
})

test_that("fastq_stage_map uncompressed with names", {
  expect_equal(
    fastq_stage_map(raw2_fastq, stage2_fastq),
    test_result2
  )
})
