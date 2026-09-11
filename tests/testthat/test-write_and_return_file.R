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

test_that("write_and_return_file.default type='qs' errors as deprecated", {
  x <- list(a = 1L, b = c("x", "y"), c = TRUE)
  file <- withr::local_tempfile(fileext = ".qs")

  expect_error(
    optimotu.pipeline::write_and_return_file(x, file),
    "deprecated"
  )
  expect_error(
    optimotu.pipeline::write_and_return_file(x, file, type = "qs"),
    "deprecated"
  )
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

test_that("write_and_return_file.data.frame type='csv' round-trips", {
  x <- data.frame(a = 1:2, b = c("x", "y"), stringsAsFactors = FALSE)
  file <- withr::local_tempfile(fileext = ".csv")
  out <- optimotu.pipeline::write_and_return_file(x, file, type = "csv")
  read_back <- readr::read_csv(file, show_col_types = FALSE)
  expect_equal(out, file)
  expect_equal(as.data.frame(read_back, stringsAsFactors = FALSE), x)
})

test_that("write_and_return_file.matrix writes tsv via data.frame", {
  x <- matrix(1:4, nrow = 2, dimnames = list(c("r1", "r2"), c("c1", "c2")))
  file <- withr::local_tempfile(fileext = ".tsv")
  out <- optimotu.pipeline::write_and_return_file(x, file, type = "tsv")
  read_back <- readr::read_tsv(file, show_col_types = FALSE)
  expect_equal(out, file)
  expect_equal(as.matrix(read_back), x, ignore_attr = TRUE)
})

test_that("write_and_return_file.matrix preserves matrix for rds", {
  x <- matrix(1:4, nrow = 2, dimnames = list(c("r1", "r2"), c("c1", "c2")))
  file <- withr::local_tempfile(fileext = ".rds")
  out <- optimotu.pipeline::write_and_return_file(x, file, type = "rds")
  read_back <- readRDS(file)
  expect_equal(out, file)
  expect_identical(read_back, x)
})

test_that("write_and_return_file.matrix preserves matrix for qs2", {
  testthat::skip_if_not_installed("qs2")
  x <- matrix(1:4, nrow = 2, dimnames = list(NULL, c("c1", "c2")))
  file <- withr::local_tempfile(fileext = ".qs2")
  out <- optimotu.pipeline::write_and_return_file(x, file, type = "qs2")
  read_back <- qs2::qs_read(file)
  expect_equal(out, file)
  expect_identical(read_back, x)
})

test_that("write_tabular_outputs writes multiple extensions", {
  withr::local_options(
    list(optimotu.pipeline.output_formats = c("rds", "tsv"))
  )
  x <- data.frame(a = 1L, b = "z", stringsAsFactors = FALSE)
  stem <- file.path(withr::local_tempdir(), "out")
  paths <- optimotu.pipeline::write_tabular_outputs(x, stem)
  expect_setequal(paths, c(paste0(stem, ".rds"), paste0(stem, ".tsv")))
  expect_identical(readRDS(paths[[1]]), x)
})

test_that("write_and_return_file.list saves named objects to RData", {
  var1 <- data.frame(x = 1L)
  var2 <- 42L
  file <- withr::local_tempfile(fileext = ".RData")
  out <- optimotu.pipeline::write_and_return_file(
    list(one = var1, two = var2),
    file,
    type = "Rdata"
  )
  env <- new.env()
  load(file, envir = env)
  expect_equal(out, file)
  expect_identical(env$one, var1)
  expect_equal(env$two, 42L)
})

test_that("write_and_return_file.list resolves string names", {
  var1 <- 1:3
  var2 <- letters[1:2]
  file <- withr::local_tempfile(fileext = ".RData")
  optimotu.pipeline::write_and_return_file(
    list("var1", "var2"),
    file,
    type = "rdata",
    envir = environment()
  )
  env <- new.env()
  load(file, envir = env)
  expect_equal(env$var1, var1)
  expect_equal(env$var2, var2)
})

test_that("write_and_return_file.list falls back to default for non-rdata type", {
  x <- list(a = 1L, b = 2L)
  file <- withr::local_tempfile(fileext = ".rds")
  out <- optimotu.pipeline::write_and_return_file(x, file, type = "rds")
  expect_equal(out, file)
  expect_identical(readRDS(file), x)
})

test_that("write_and_return_file.data.frame type='parquet' round-trips", {
  testthat::skip_if_not_installed("arrow")
  x <- data.frame(a = 1:2, b = c("x", "y"), stringsAsFactors = FALSE)
  file <- withr::local_tempfile(fileext = ".parquet")
  out <- optimotu.pipeline::write_and_return_file(x, file, type = "parquet")
  read_back <- as.data.frame(arrow::read_parquet(file))
  expect_equal(out, file)
  expect_equal(read_back, x)
})

test_that("write_and_return_file.data.frame type='fst' round-trips", {
  testthat::skip_if_not_installed("fst")
  x <- data.frame(a = 1:2, b = c("x", "y"), stringsAsFactors = FALSE)
  file <- withr::local_tempfile(fileext = ".fst")
  out <- optimotu.pipeline::write_and_return_file(x, file, type = "fst")
  read_back <- as.data.frame(fst::read_fst(file))
  expect_equal(out, file)
  expect_equal(read_back, x)
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
