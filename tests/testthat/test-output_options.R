test_that("normalize_output_formats defaults and case-folds", {
  expect_equal(
    optimotu.pipeline:::normalize_output_formats(NULL),
    c("rds", "tsv")
  )
  expect_equal(
    optimotu.pipeline:::normalize_output_formats(c("RDS", "TSV", "CSV")),
    c("rds", "tsv", "csv")
  )
  expect_equal(
    optimotu.pipeline:::normalize_output_formats("qd"),
    "qdata"
  )
})

test_that("normalize_output_formats accepts parquet", {
  expect_equal(
    optimotu.pipeline:::normalize_output_formats("parquet"),
    "parquet"
  )
})

test_that("normalize_output_formats rejects unknown formats", {
  expect_error(
    optimotu.pipeline:::normalize_output_formats("parquet_invalid"),
    "Unknown output format"
  )
})

test_that("parse_output_options sets formats and wide_table", {
  withr::local_options(
    list(
      optimotu.pipeline.output_formats = NULL,
      optimotu.pipeline.wide_table = NULL
    )
  )
  optimotu.pipeline:::parse_output_options(list())
  expect_equal(optimotu.pipeline::output_formats(), c("rds", "tsv"))
  expect_false(optimotu.pipeline::do_wide_otu_table())

  withr::local_options(
    list(
      optimotu.pipeline.output_formats = NULL,
      optimotu.pipeline.wide_table = NULL
    )
  )
  optimotu.pipeline:::parse_output_options(list(
    output = list(
      formats = list("csv", "Rdata"),
      wide_table = TRUE
    )
  ))
  expect_equal(
    optimotu.pipeline::output_formats(),
    c("csv", "rdata")
  )
  expect_equal(optimotu.pipeline::output_table_formats(), "csv")
  expect_true(optimotu.pipeline::do_output_rdata())
  expect_true(optimotu.pipeline::do_wide_otu_table())
})

test_that("parse_output_options honors top-level wide_table", {
  withr::local_options(
    list(
      optimotu.pipeline.output_formats = NULL,
      optimotu.pipeline.wide_table = NULL
    )
  )
  optimotu.pipeline:::parse_output_options(list(wide_table = TRUE))
  expect_true(optimotu.pipeline::do_wide_otu_table())
})

test_that("parse_output_options dense_table alias enables wide table", {
  withr::local_options(
    list(
      optimotu.pipeline.output_formats = NULL,
      optimotu.pipeline.wide_table = NULL
    )
  )
  optimotu.pipeline:::parse_output_options(list(dense_table = TRUE))
  expect_true(optimotu.pipeline::do_wide_otu_table())
})

test_that("parse_output_options errors when wide_table and dense_table both set", {
  expect_error(
    optimotu.pipeline:::parse_output_options(list(
      wide_table = TRUE,
      dense_table = TRUE
    )),
    "Only one of 'wide_table' and 'dense_table'"
  )
})

test_that("parse_otu_table_options delegates to parse_output_options", {
  withr::local_options(
    list(
      optimotu.pipeline.output_formats = NULL,
      optimotu.pipeline.wide_table = NULL
    )
  )
  optimotu.pipeline::parse_otu_table_options(list(
    output = list(formats = "fst")
  ))
  expect_equal(optimotu.pipeline::output_formats(), "fst")
})
