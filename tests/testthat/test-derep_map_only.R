derep1 <- dada2::derepFastq(raw1_fastq_gz)

test_that("derep_map_only works on a single derep", {
  derep2 <- rlang::duplicate(derep1)
  expect_identical(derep1, derep2)
  expect_identical(object.size(derep1), object.size(derep2))
  derep_map_only(derep2)
  expect_identical(derep1$map, derep2$map)
  expect_lt(object.size(derep2), object.size(derep1))
})


derep_multi1 <- dada2::derepFastq(stage1_fastq_gz)
test_that("derep_map_only works on a multi-file derep", {
  derep_multi2 <- rlang::duplicate(derep_multi1)
  expect_identical(derep_multi1, derep_multi2)
  expect_identical(object.size(derep_multi1), object.size(derep_multi2))
  derep_map_only(derep_multi2)
  expect_identical(derep_multi1$map, derep_multi2$map)
  expect_lt(object.size(derep_multi2), object.size(derep_multi1))
})
