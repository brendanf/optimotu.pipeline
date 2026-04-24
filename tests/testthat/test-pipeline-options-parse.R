test_that("unnest_yaml_list flattens length-1 nested lists", {
  nested <- list(list(a = 1), list(b = 2))
  out <- optimotu.pipeline:::unnest_yaml_list(nested)
  expect_equal(out, list(a = 1, b = 2))
})

test_that("parse_project_name handles default, valid, and invalid names", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())

  expect_warning(
    out <- optimotu.pipeline:::parse_project_name(list()),
    "Missing project name"
  )
  expect_equal(out, "metabarcoding_project")
  expect_equal(optimotu.pipeline::project_name(), "metabarcoding_project")

  optimotu.pipeline:::parse_project_name(list(project_name = "my_project-1"))
  expect_equal(optimotu.pipeline::project_name(), "my_project-1")

  expect_error(
    optimotu.pipeline:::parse_project_name(list(project_name = "bad name")),
    "Project name should consist"
  )
})

test_that("parse_orient sets default and explicit orientation", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())

  expect_message(
    optimotu.pipeline:::parse_orient(list(orient = NULL)),
    "No read orientation specified"
  )
  expect_equal(optimotu.pipeline::read_orientation(), "fwd")

  optimotu.pipeline:::parse_orient(list(orient = "mixed"))
  expect_equal(optimotu.pipeline::read_orientation(), "mixed")

  expect_error(
    optimotu.pipeline:::parse_orient(list(orient = "forward")),
    "Must comply to pattern"
  )
})

test_that("parse_custom_sample_table accepts FALSE and valid path", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())

  expect_null(optimotu.pipeline::custom_sample_table())
  optimotu.pipeline:::parse_custom_sample_table(
    list(custom_sample_table = FALSE)
  )
  expect_null(optimotu.pipeline::custom_sample_table())

  sample_tbl <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("sample\tseqrun", "s1\tr1"), sample_tbl)
  optimotu.pipeline:::parse_custom_sample_table(
    list(custom_sample_table = sample_tbl)
  )
  expect_equal(optimotu.pipeline::custom_sample_table(), sample_tbl)
  expect_true(optimotu.pipeline:::do_custom_sample_table())
})

test_that("parse_added_reference requires both files and sets options", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())

  fasta <- withr::local_tempfile(fileext = ".fa")
  table <- withr::local_tempfile(fileext = ".xlsx")
  writeLines(">x\nACGT", fasta)
  writeLines("placeholder", table)

  expect_error(
    optimotu.pipeline:::parse_added_reference(
      list(added_reference = list(fasta = fasta))
    ),
    "both must be given"
  )

  optimotu.pipeline:::parse_added_reference(
    list(added_reference = list(list(fasta = fasta), list(table = table)))
  )
  expect_true(optimotu.pipeline::do_added_reference())
  expect_equal(optimotu.pipeline::added_reference_fasta(), fasta)
  expect_equal(optimotu.pipeline::added_reference_table(), table)
})

test_that("parse_parallel_options handles workers precedence and bounds", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())

  expect_warning(
    optimotu.pipeline:::parse_parallel_options(
      list(
        local_threads = 3L,
        max_batchsize = 1000L,
        workers_per_seqrun = 4L,
        jobs_per_seqrun = 2L,
        min_workers = 2L,
        max_workers = 10L
      )
    ),
    "both 'workers_per_seqrun' and 'jobs_per_seqrun'"
  )

  expect_equal(getOption("optimotu_num_threads"), 3L)
  expect_equal(optimotu.pipeline::max_batchsize(), 1000L)
  expect_equal(optimotu.pipeline::workers_per_seqrun(), 4L)
  expect_equal(optimotu.pipeline::min_workers(), 2L)
  expect_equal(optimotu.pipeline::max_workers(), 10L)
})

test_that("parse_rarefy_options supports number and fraction modes", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())

  optimotu.pipeline:::parse_rarefy_options(
    list(rarefy = list(number = 100L))
  )
  expect_equal(optimotu.pipeline::rarefy_number(), 100L)
  expect_true(optimotu.pipeline::do_rarefy())

  optimotu.pipeline:::parse_rarefy_options(
    list(rarefy = list(numerator = 1L, denominator = 10L))
  )
  expect_equal(optimotu.pipeline::rarefy_numerator(), 1L)
  expect_equal(optimotu.pipeline::rarefy_denominator(), 10L)

  expect_error(
    optimotu.pipeline:::parse_rarefy_options(
      list(rarefy = list(number = 100L, numerator = 1L))
    ),
    "cannot be used in conjunction"
  )
})
