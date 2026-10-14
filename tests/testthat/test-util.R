test_that("is_same_path matches relative and absolute paths", {
  fa <- withr::local_tempfile(fileext = ".fasta")
  writeLines(">s\nACGT", fa)
  expect_true(is_same_path(fa, fa))
  expect_true(is_same_path(fa, normalizePath(fa)))
  expect_false(is_same_path(fa, tempfile()))
  expect_false(is_same_path(NULL, fa))
  expect_false(is_same_path(fa, NA_character_))
})
