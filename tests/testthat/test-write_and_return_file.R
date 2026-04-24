test_that("ensure_directory creates parent directory and returns file invisibly", {
  file <- file.path(withr::local_tempdir(), "nested", "path", "out.txt")
  expect_false(dir.exists(dirname(file)))

  out <- optimotu.pipeline::ensure_directory(file)

  expect_true(dir.exists(dirname(file)))
  expect_equal(out, file)
})

test_that("write_and_return_file.character round-trips content and returns filename", {
  x <- c("alpha", "beta", "gamma")
  file <- withr::local_tempfile(fileext = ".txt")

  out <- optimotu.pipeline::write_and_return_file(x, file)

  expect_equal(out, file)
  expect_equal(readLines(file), x)
})

test_that("write_and_return_file.XStringSet round-trips sequences and returns filename", {
  x <- Biostrings::DNAStringSet(c(seq1 = "ACGT", seq2 = "GGTTAA"))
  file <- withr::local_tempfile(fileext = ".fa")

  out <- optimotu.pipeline::write_and_return_file(x, file)
  read_back <- Biostrings::readDNAStringSet(file)

  expect_equal(out, file)
  expect_identical(read_back, x)
})

test_that("write_and_return_file.data.frame type='rds' round-trips and returns filename", {
  x <- data.frame(a = 1:3, b = c("x", "y", "z"), stringsAsFactors = FALSE)
  file <- withr::local_tempfile(fileext = ".rds")

  out <- optimotu.pipeline::write_and_return_file(x, file, type = "rds")
  read_back <- readRDS(file)

  expect_equal(out, file)
  expect_identical(read_back, x)
})

test_that("write_and_return_file.data.frame type='tsv' round-trips and returns filename", {
  x <- data.frame(
    a = c(1.1, 2.2, 3.3),
    b = c("x", "y", "z"),
    stringsAsFactors = FALSE
  )
  file <- withr::local_tempfile(fileext = ".tsv")

  out <- optimotu.pipeline::write_and_return_file(x, file, type = "tsv")
  read_back <- readr::read_tsv(file, show_col_types = FALSE)

  expect_equal(out, file)
  expect_identical(as.data.frame(read_back, stringsAsFactors = FALSE), x)
})

test_that("write_and_return_file.default type='rds' round-trips and returns filename", {
  x <- list(a = 1L, b = c("x", "y"), c = TRUE)
  file <- withr::local_tempfile(fileext = ".rds")

  out <- optimotu.pipeline::write_and_return_file(x, file)
  read_back <- readRDS(file)

  expect_equal(out, file)
  expect_identical(read_back, x)
})

test_that("write_and_return_file.default type='qs' round-trips and returns filename", {
  testthat::skip_if_not_installed("qs")
  x <- list(a = 1L, b = c("x", "y"), c = TRUE)
  file <- withr::local_tempfile(fileext = ".qs")

  out <- optimotu.pipeline::write_and_return_file(x, file)
  read_back <- qs::qread(file)

  expect_equal(out, file)
  expect_identical(read_back, x)
})

test_that("write_and_return_file.default type='qs2' round-trips and returns filename", {
  testthat::skip_if_not_installed("qs2")
  x <- list(a = 1L, b = c("x", "y"), c = TRUE)
  file <- withr::local_tempfile(fileext = ".qs2")

  out <- optimotu.pipeline::write_and_return_file(x, file)
  read_back <- qs2::qs_read(file)

  expect_equal(out, file)
  expect_identical(read_back, x)
})

test_that("write_and_return_file.default type='qd' round-trips and returns filename", {
  testthat::skip_if_not_installed("qs2")
  x <- list(a = 1L, b = c("x", "y"), c = TRUE)
  file <- withr::local_tempfile(fileext = ".qd")

  out <- optimotu.pipeline::write_and_return_file(x, file)
  read_back <- qs2::qd_read(file)

  expect_equal(out, file)
  expect_identical(read_back, x)
})

test_that("write_and_return_file.default type='qdata' round-trips and returns filename", {
  testthat::skip_if_not_installed("qs2")
  x <- list(a = 1L, b = c("x", "y"), c = TRUE)
  file <- withr::local_tempfile(fileext = ".qdata")

  out <- optimotu.pipeline::write_and_return_file(x, file)
  read_back <- qs2::qd_read(file)

  expect_equal(out, file)
  expect_identical(read_back, x)
})

test_that("write_and_return_file.ggplot returns filename and writes readable image", {
  testthat::skip_if_not_installed("png")
  x <- ggplot2::ggplot(mtcars, ggplot2::aes(x = wt, y = mpg)) +
    ggplot2::geom_point()
  file <- withr::local_tempfile(fileext = ".png")

  out <- optimotu.pipeline::write_and_return_file(
    x,
    file,
    width = 4,
    height = 3
  )
  read_back <- png::readPNG(file)

  expect_equal(out, file)
  expect_true(file.exists(file))
  expect_true(is.array(read_back))
})
